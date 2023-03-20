use chrono::{DateTime, Utc};
use clap::{Args, Parser, Subcommand};
use motif_finder::gibbs_sampler::iterate_gibbs_sampler;
use motif_finder::median_string::median_string;
use motif_finder::randomized_motif_search::iterate_randomized_motif_search;
use motif_finder::{consensus_string, Error};
use std::fs::{self, File};
use std::io::{self, BufRead, Write};
use std::path::Path;

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

fn main() -> Result<(), Error> {
    let args = Cli::parse();
    let dt = Utc::now();
    let start_time: i64 = dt.timestamp_micros();
    println!("start at {}", dt.format("%Y-%m-%d %H:%M:%S"));
    let sequences = load_data(&args.global_opts.input_file, args.global_opts.num_entries)?;
    let GlobalOpts { k, .. } = args.global_opts;
    let motifs = match args.command {
        Commands::GibbsSampler {
            num_iterations,
            num_runs,
        } => run_gibbs_sampler(&sequences, k, num_runs, num_iterations),
        Commands::MedianString => run_median_string(&sequences, &args),
        Commands::Randomized { num_runs } => run_randomized_motif_search(&sequences, k, num_runs),
    }?;
    let consensus_string = generate_consensus_string(&motifs, k)?;
    if let Some(save_flag) = &args.global_opts.output_file {
        let (mut file, file_path) =
            create_output_file(save_flag, &args, &consensus_string, dt, start_time)?;
        write_motifs(&mut file, &motifs)?;
        println!("Saved to file: {}", file_path);
    }
    // for motif in motifs {
    //     println!("{}", motif);
    // }
    println!("Consensus string: {}", consensus_string);
    let dt_end = Utc::now();
    println!("End at {}", dt_end.format("%Y-%m-%d %H:%M:%S"));
    if let Some(duration) = dt_end.signed_duration_since(dt).num_microseconds() {
        println!("Done in {} seconds", duration as f64 / 1_000_000.0);
    }

    Ok(())
}

fn read_lines<P>(filename: P) -> io::Result<io::Lines<io::BufReader<File>>>
where
    P: AsRef<Path>,
{
    let file = File::open(filename)?;
    Ok(io::BufReader::new(file).lines())
}
fn load_data(path_to_file: &str, num_entries: usize) -> Result<Vec<String>, Error> {
    println!("Loading data from '{}'...", path_to_file);
    let mut sequences = vec![];
    let mut current_sequence = String::new();
    let mut sequence_count = 0;
    if let Ok(lines) = read_lines(path_to_file) {
        for line in lines {
            if let Ok(line) = line {
                if line.starts_with('>') {
                    sequence_count += 1;
                    if sequence_count > num_entries {
                        break;
                    }
                    if !current_sequence.is_empty() {
                        sequences.push(current_sequence.clone());
                        current_sequence.clear();
                    }
                } else {
                    current_sequence.push_str(&line.to_uppercase());
                }
            } else {
                return Err(Error::IOError);
            }
        }
        if !current_sequence.is_empty() {
            sequences.push(current_sequence.clone());
        }
    } else {
        return Err(Error::FileNotFoundError);
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
    iterate_gibbs_sampler(sequences, k, sequences.len(), num_iterations, num_runs)
}

fn run_median_string(sequences: &[String], args: &Cli) -> Result<Vec<String>, Error> {
    let median_string = median_string(args.global_opts.k, sequences);
    println!("median string: {}", median_string);
    let vec = vec![median_string];
    Ok(vec)
}

fn run_randomized_motif_search(
    sequences: &[String],
    k: usize,
    num_runs: usize,
) -> Result<Vec<String>, Error> {
    iterate_randomized_motif_search(sequences, k, num_runs)
}

fn write_file_header(
    file: &mut fs::File,
    args_cli: &Cli,
    dt: DateTime<Utc>,
    consensus_string: &str,
) -> io::Result<()> {
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
    let dt = Utc::now();
    writeln!(file, "End time: {}", dt.format("%Y-%m-%d %H:%M:%S"))?;
    writeln!(file, "Consensus string: {}", consensus_string)?;
    writeln!(
        file,
        "_______________________________________________________________________________________"
    )?;
    Ok(())
}

fn create_output_file(
    save_flag: &Option<String>,
    args: &Cli,
    consensus_string: &str,
    dt: DateTime<Utc>,
    timestamp: i64,
) -> Result<(File, String), Error> {
    let save_path: String = save_flag
        .clone()
        .unwrap_or_else(|| format!("MotifFinder-output-{timestamp}-{}.txt", args.global_opts.k));
    let mut file = match fs::File::create(&save_path) {
        Ok(file) => file,
        Err(_err) => return Err(Error::IOError),
    };
    match write_file_header(&mut file, args, dt, consensus_string) {
        Ok(()) => {}
        Err(_err) => return Err(Error::IOError),
    };
    Ok((file, save_path))
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
        return Err(Error::InvalidMotifError);
    } else if motifs.len() == 1 {
        return Ok(motifs[0].clone());
    }
    consensus_string(motifs, k)
}
