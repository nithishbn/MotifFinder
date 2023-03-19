use std::fs::{self, File};
use std::io::{self, BufRead, Write};
use std::path::Path;

use chrono::{DateTime, Utc};
use clap::{Args, Parser, Subcommand};
use motif_finder::gibbs_sampler::iterate_gibbs_sampler;
use motif_finder::median_string::median_string;
use motif_finder::randomized_motif_search::iterate_randomized_motif_search;

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
    #[arg(short, long = "input")]
    input_file: String,

    /// how many entries to read
    #[arg(short = 'e', long = "entries")]
    num_entries: usize,

    /// motif length
    #[arg(short)]
    k: usize,

    /// save motifs to file
    #[arg(short, long = "output")]
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
        about = "Run the Median String algorithm (Warning: this can take a long time to run for large values of k))"
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

#[derive(Debug)]
pub enum Error {
    GenericError,
    IOError,
}

fn main() -> Result<(), Error> {
    let args = Cli::parse();
    let dt = Utc::now();
    let start_time: i64 = dt.timestamp();
    println!("start at {}", start_time);
    let value = match &args.command {
        Commands::GibbsSampler { .. } => run_gibbs_sampler(&args, dt, start_time),
        Commands::MedianString => run_median_string(&args, dt, start_time),
        Commands::Randomized { .. } => run_randomized_motif_search(&args, dt, start_time),
    };
    if let Err(_err) = value {
        return Err(Error::GenericError);
    }
    let dt = Utc::now();
    println!("done in {} seconds", dt.timestamp() - start_time);
    Ok(())
}

fn read_lines<P>(filename: P) -> io::Result<io::Lines<io::BufReader<File>>>
where
    P: AsRef<Path>,
{
    let file = File::open(filename)?;
    Ok(io::BufReader::new(file).lines())
}
fn load_data(path_to_file: &str, num_entries: usize) -> Result<Vec<String>, &str> {
    let mut genome = vec![];
    let mut chromosome = String::from("");
    let mut current = 1;
    println!("Loading data from '{}'...", path_to_file);
    if let Ok(lines) = read_lines(path_to_file) {
        for line in lines {
            if let Ok(line) = line {
                if line.starts_with(">") {
                    if chromosome.chars().count() > 0 {
                        if chromosome.chars().count() == 0 {
                            continue;
                        }
                        genome.push(chromosome);
                        chromosome = "".to_string();
                        current += 1;
                        continue;
                    } else if chromosome.chars().count() == 0 {
                        continue;
                    }
                }
                if current == num_entries + 1 {
                    if chromosome.chars().count() == 0 {
                        continue;
                    }
                    genome.push(chromosome);
                    chromosome = "".to_string();
                    break;
                }

                chromosome += line.to_uppercase().trim();
            }
        }
    } else {
        return Err("Error loading data");
    }
    if chromosome.chars().count() > 0 {
        genome.push(chromosome);
    }
    println!("Done loading data: {} entries", genome.len());
    Ok(genome)
}
fn run_gibbs_sampler(args_cli: &Cli,dt:DateTime<Utc>,start_time: i64) -> Result<(), Error> {
    let (num_runs, num_iterations) = if let Commands::GibbsSampler {
        num_runs,
        num_iterations,
    } = &args_cli.command
    {
        (*num_runs, *num_iterations)
    } else {
        return Err(Error::GenericError);
    };
    let Cli {
        global_opts: args,
        command: _,
    } = args_cli;

    println!("start at {}", start_time);
    let genome = load_data(&args.input_file, args.num_entries);
    let genome = match genome {
        Ok(genome) => genome,
        Err(_err) => return Err(Error::IOError),
    };
    let motifs = iterate_gibbs_sampler(&genome, args.k, genome.len(), num_iterations, num_runs);
    if let Some(save_flag) = &args.output_file {
        let (mut file,file_path) = create_output_file(&save_flag, args_cli, dt, start_time)?;
        match write_motifs(&mut file, &motifs) {
            Ok(()) => {}
            Err(_err) => return Err(Error::IOError),
        }
        println!("Saved to file: {}", file_path);
    }
    for motif in motifs {
        println!("{}", motif);
    }
    
    Ok(())
}

fn run_median_string(args: &Cli, dt: DateTime<Utc>, start_time: i64) -> Result<(), Error> {
    let genome = load_data(&args.global_opts.input_file, args.global_opts.num_entries);
    let genome = match genome {
        Ok(genome) => genome,
        Err(_err) => return Err(Error::IOError),
    };
    let median_string = median_string(args.global_opts.k, &genome);
    println!("median string: {}", median_string);
    if let Some(save_flag) = &args.global_opts.output_file {
        // save flag exists?
        let (mut file,file_path) = create_output_file(save_flag, args, dt, start_time)?;
        match writeln!(file, "median string: {}", median_string) {
            Ok(()) => {}
            Err(_err) => return Err(Error::IOError),
        }
        println!("Saved to file: {}", file_path);
    }

    Ok(())
}

fn run_randomized_motif_search(
    args_cli: &Cli,
    dt: DateTime<Utc>,
    start_time: i64,
) -> Result<(), Error> {
    let num_runs = if let Commands::Randomized { num_runs } = &args_cli.command {
        *num_runs
    } else {
        return Err(Error::GenericError);
    };
    let Cli {
        global_opts: args,
        command: _,
    } = args_cli;

    
    let genome = load_data(&args.input_file, args.num_entries);
    let genome = match genome {
        Ok(genome) => genome,
        Err(_err) => return Err(Error::IOError),
    };
    let motifs = iterate_randomized_motif_search(&genome, args.k, num_runs);

    if let Some(save_flag) = &args.output_file {
        let (mut file,file_path) = create_output_file(save_flag, args_cli, dt, start_time)?;
        match write_motifs(&mut file, &motifs) {
            Ok(()) => {}
            Err(_err) => return Err(Error::IOError),
        }
        println!("Saved to file: {}", file_path);
    }
    for motif in motifs {
        println!("{}", motif);
    }
    let dt = Utc::now();
    println!("done at {}", dt.timestamp() - start_time);
    Ok(())
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
    writeln!(file, "Start time: {}", dt)?;
    let dt = Utc::now();
    writeln!(file, "End time: {}", dt)?;
    writeln!(
        file,
        "_______________________________________________________________________________________"
    )?;
    Ok(())
}

fn create_output_file(
    save_flag: &Option<String>,
    args: &Cli,
    dt: DateTime<Utc>,
    timestamp: i64,
) -> Result<(File,String), Error> {
    let save_path: String = save_flag.clone()
        .unwrap_or_else(||format!("MotifFinder-output-{timestamp}-{}.txt", args.global_opts.k));
    let mut file = match fs::File::create(&save_path) {
        Ok(file) => file,
        Err(_err) => return Err(Error::IOError),
    };
    match write_file_header(&mut file, args, dt) {
        Ok(()) => {}
        Err(_err) => return Err(Error::IOError),
    };
    Ok((file,save_path))
}

fn write_motifs(file: &mut fs::File, motifs: &Vec<String>) -> Result<(), Error> {
    for (i, motif) in motifs.iter().enumerate() {
        let motif = motif.trim();
        write!(file, ">motif {}\n", i + 1).map_err(|_| Error::IOError)?;
        write!(file, "{}\n", motif).map_err(|_| Error::IOError)?;
    }

    Ok(())
}
