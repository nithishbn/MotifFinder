use bio::io::fasta;
use chrono::{DateTime, Utc};
use clap::{Args, Parser, Subcommand};
use motif_finder::alignment::local_alignment_score_only;
use motif_finder::gibbs_sampler::iterate_gibbs_sampler;
use motif_finder::median_string::median_string;
use motif_finder::randomized_motif_search::iterate_randomized_motif_search;
use motif_finder::{consensus_string, Error};
use rayon::prelude::*;
use std::fs::{self, File};
use std::io::{self, Write};
use std::str;
/// Motif Finder
#[derive(Parser)]
#[command(author, version, about, long_about = None)]
struct Cli {
    #[clap(flatten)]
    global_opts: GlobalOpts,
    #[command(subcommand)]
    command: Commands,
}

#[derive(Debug, Args)]
struct GlobalOpts {
    /// file path of the genome
    input_file: String,

    /// how many entries to read
    #[arg(short = 'e', long = "entries")]
    num_entries: usize,

    /// motif length
    #[arg(short)]
    k: usize,

    /// alignment
    #[arg(short = 'a', long = "align")]
    align: bool,

    /// save motifs to file
    #[arg(short = 'o', long = "output")]
    output_file: Option<Option<String>>,
}

#[derive(Subcommand)]
enum Commands {
    #[clap(name = "gibbs", about = "Run the Gibbs Sampler algorithm")]
    GibbsSampler {
        /// number of runs
        #[arg(short = 'r', long = "runs")]
        num_runs: usize,

        /// number of iterations per run
        #[arg(short = 't', long = "iters")]
        num_iterations: usize,
    },

    #[clap(
        name = "median",
        about = "Run the Median String algorithm (Warning: this can take a long time to run for large values of k)"
    )]
    MedianString,

    #[clap(
        name = "randomized",
        about = "Run the Randomized Motif Search algorithm"
    )]
    Randomized {
        /// number of runs
        #[arg(short = 'r', long = "runs")]
        num_runs: usize,
    },
}
struct Summary {
    consensus_string: String,
    best_motif: Option<String>,
    best_motif_score: Option<isize>,
}
fn main() -> Result<(), Error> {
    let mut args = Cli::parse();
    let dt = Utc::now();
    let start_time: i64 = dt.timestamp_micros();
    println!("start at {}", dt.format("%Y-%m-%d %H:%M:%S"));
    let sequences = load_data(&args.global_opts.input_file, args.global_opts.num_entries)?;
    args.global_opts.num_entries = sequences.len();
    let GlobalOpts { k, .. } = args.global_opts;
    if k == 0 {
        return Err(Error::InvalidMotifLength);
    }
    let (file, file_path) = if let Some(save_flag) = &args.global_opts.output_file {
        let (mut file, file_path) = create_output_file(save_flag, k, start_time)?;
        match write_file_header(&mut file, &args, dt) {
            Ok(()) => {}
            Err(_err) => return Err(Error::IOError),
        }
        (Some(file), Some(file_path))
    } else {
        (None, None)
    };

    let motifs = match args.command {
        Commands::GibbsSampler {
            num_iterations,
            num_runs,
        } => run_gibbs_sampler(&sequences, k, num_runs, num_iterations),
        Commands::MedianString => run_median_string(&sequences, &args),
        Commands::Randomized { num_runs } => run_randomized_motif_search(&sequences, k, num_runs),
    }?;

    let consensus_string = generate_consensus_string(&motifs, k)?;
    println!("Consensus string: {}", consensus_string);
    let motifs_clone = motifs.clone();
    let (best_motif_score, best_motif) = if args.global_opts.align {
        let top_five = align_motifs_multi_threaded(sequences, motifs)?;
        println!("Top 5 motifs:");
        for (score, motif) in &top_five {
            println!("{}: {}", score, motif);
        }
        let (best_motif_score, best_motif) = top_five[0].clone();
        (Some(best_motif_score), Some(best_motif))
    } else {
        (None, None)
    };
    if let Some(mut file) = file {
        let summary = Summary {
            consensus_string,
            best_motif,
            best_motif_score,
        };
        match output_results_to_file(&mut file, &motifs_clone, &summary) {
            Ok(()) => {
                println!("Results saved to {}", file_path.ok_or(Error::IOError)?);
            }
            Err(_err) => {
                return Err(Error::IOError);
            }
        }
    }

    let dt_end = Utc::now();
    println!("End at {}", dt_end.format("%Y-%m-%d %H:%M:%S"));
    if let Some(duration) = dt_end.signed_duration_since(dt).num_microseconds() {
        println!("Done in {} seconds", duration as f64 / 1_000_000.0);
    }

    Ok(())
}

fn align_motifs_multi_threaded(
    sequences: Vec<String>,
    motifs: Vec<String>,
) -> Result<Vec<(isize, String)>, Error> {
    println!("Aligning motifs to sequences...");
    let mut top_five: Vec<(isize, String)> = motifs
        .par_iter()
        .map(|motif| {
            let mut highest_score = 0;
            for sequence in sequences.iter() {
                highest_score += local_alignment_score_only(sequence, motif, 1, 0, -10);
            }
            (highest_score, motif.to_owned())
        })
        .collect();

    top_five.par_sort_by(|a, b| b.0.cmp(&a.0));
    top_five.dedup();
    Ok(top_five[0..5].to_vec())
}

fn load_data(path_to_file: &str, num_entries: usize) -> Result<Vec<String>, Error> {
    println!("Loading data from '{}'...", path_to_file);
    let mut sequences = vec![];
    let file = match File::open(path_to_file) {
        Ok(file) => file,
        Err(_) => return Err(Error::FileNotFoundError),
    };
    let mut records = fasta::Reader::new(file).records();
    let mut count = 0;
    while let Some(Ok(record)) = records.next() {
        count += 1;
        if count > num_entries {
            break;
        }
        let s = match str::from_utf8(record.seq()) {
            Ok(v) => v,
            Err(_e) => return Err(Error::InvalidSequence),
        }
        .to_string()
        .to_uppercase();

        sequences.push(s);
    }
    println!("Done loading data: {} entries", sequences.len());
    Ok(sequences)
}
fn run_gibbs_sampler(
    sequences: &Vec<String>,
    k: usize,
    num_runs: usize,
    num_iterations: usize,
) -> Result<Vec<String>, Error> {
    if num_runs == 0 {
        return Err(Error::InvalidNumberOfRuns);
    }
    if num_iterations == 0 {
        return Err(Error::InvalidNumberOfIterations);
    }

    iterate_gibbs_sampler(sequences, k, sequences.len(), num_iterations, num_runs)
}

fn run_median_string(sequences: &[String], args: &Cli) -> Result<Vec<String>, Error> {
    let median_string = median_string(args.global_opts.k, sequences)?;
    println!("median string: {}", median_string);
    let vec = vec![median_string];
    Ok(vec)
}

fn run_randomized_motif_search(
    sequences: &[String],
    k: usize,
    num_runs: usize,
) -> Result<Vec<String>, Error> {
    if num_runs == 0 {
        return Err(Error::InvalidNumberOfRuns);
    }
    iterate_randomized_motif_search(sequences, k, num_runs)
}

fn write_file_header(file: &mut fs::File, args_cli: &Cli, dt: DateTime<Utc>) -> io::Result<()> {
    let Cli {
        global_opts: args,
        command: _,
    } = args_cli;
    let version = env!("CARGO_PKG_VERSION");
    writeln!(file, "MotifFinder {}", version)?;
    let command = match &args_cli.command {
        Commands::Randomized { .. } => "Randomized",
        Commands::GibbsSampler { .. } => "Gibbs Sampler",
        Commands::MedianString => "Median String",
    };
    writeln!(file, "Command: {}", command)?;
    writeln!(file, "k: {}", args.k)?;
    writeln!(file, "number of entries: {}", args.num_entries)?;
    match &args_cli.command {
        Commands::Randomized { num_runs } => {
            writeln!(file, "runs: {}", num_runs)?;
        }
        Commands::GibbsSampler {
            num_runs,
            num_iterations,
        } => {
            writeln!(file, "runs: {}", num_runs)?;
            writeln!(file, "iterations: {}", num_iterations)?;
        }
        Commands::MedianString => {}
    }
    writeln!(file, "Start time: {}", dt.format("%Y-%m-%d %H:%M:%S"))?;

    Ok(())
}

fn create_output_file(
    save_flag: &Option<String>,
    k: usize,
    timestamp: i64,
) -> Result<(File, String), Error> {
    let save_path: String = save_flag
        .clone()
        .unwrap_or_else(|| format!("MotifFinder-output-{timestamp}-{}.txt", k));
    let file = match fs::File::create(&save_path) {
        Ok(file) => file,
        Err(_err) => return Err(Error::IOError),
    };
    Ok((file, save_path))
}
fn output_results_to_file(
    file: &mut fs::File,
    motifs: &[String],
    summary: &Summary,
) -> Result<(), Error> {
    let Summary {
        consensus_string,
        best_motif_score,
        best_motif,
    } = summary;

    writeln!(file, "Consensus string: {}", consensus_string).map_err(|_| Error::IOError)?;
    if let Some(best_motif) = best_motif {
        writeln!(file, "Best motif: {}", best_motif).map_err(|_| Error::IOError)?;
    }
    if let Some(best_motif_score) = best_motif_score {
        writeln!(file, "Best motif score: {}", best_motif_score).map_err(|_| Error::IOError)?;
    }

    writeln!(
        file,
        "_________________________________________________________________________________________"
    )
    .map_err(|_| Error::IOError)?;
    write_motifs(file, motifs)?;
    Ok(())
}
fn write_motifs(file: &mut fs::File, motifs: &[String]) -> Result<(), Error> {
    for (i, motif) in motifs.iter().enumerate() {
        let motif = motif.trim();
        writeln!(file, ">motif {}", i + 1).map_err(|_| Error::IOError)?;
        writeln!(file, "{}", motif).map_err(|_| Error::IOError)?;
    }

    Ok(())
}

fn generate_consensus_string(motifs: &[String], k: usize) -> Result<String, Error> {
    if motifs.is_empty() {
        return Err(Error::NoMotifsFound);
    } else if motifs.len() == 1 {
        return Ok(motifs[0].clone());
    }
    consensus_string(motifs, k)
}

#[cfg(test)]
mod test {
    #[test]
    pub fn test_load_data() {
        let sequences = super::load_data("promoters.fasta", 4).unwrap();
        assert_eq!(sequences.len(), 4);
        let sequences = super::load_data("promoters.fasta", 3).unwrap();
        assert_eq!(sequences.len(), 3);
        let sequences = super::load_data("promoters.fasta", 2).unwrap();
        assert_eq!(sequences.len(), 2);
        let sequences = super::load_data("promoters.fasta", 1).unwrap();
        assert_eq!(sequences.len(), 1);
        let sequences = super::load_data("promoters.fasta", 0).unwrap();
        assert_eq!(sequences.len(), 0);
    }
}
