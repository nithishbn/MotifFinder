use std::ops::RangeInclusive;

use crate::{
    align_motifs_multi_threaded,
    alignment::align_motifs_distance,
    generate_consensus_string, load_data, run_gibbs_sampler, run_median_string,
    run_randomized_motif_search, unique_motifs,
    utils::{
        create_output_file, generate_vector_space_delimited, output_results_to_file,
        write_file_header,
    },
    Error,
};
use chrono::Utc;
use clap::{Args, Parser, Subcommand};
use clap_verbosity_flag::InfoLevel;
use rayon::prelude::*;
use tracing::{error, info, trace};
/// Motif Finder
#[derive(Debug, Parser)]
#[command(author, version, about, long_about = None)]
pub struct MotifFinder {
    #[clap(flatten)]
    global_opts: GlobalOpts,
    #[command(subcommand)]
    pub command: Commands,
    #[clap(flatten)]
    pub verbose: clap_verbosity_flag::Verbosity<InfoLevel>,
}

impl MotifFinder {
    pub fn exec(mut self) -> Result<(), Error> {
        let dt = Utc::now();
        let start_time: i64 = dt.timestamp_micros();
        println!("Welcome to MotifFinder!");
        let sequences = load_data(&self.global_opts.input_file, self.global_opts.num_entries)?;
        self.global_opts.num_entries = sequences.len();
        let GlobalOpts { k, .. } = self.global_opts;

        let (file, file_path) = if let Some(save_flag) = &self.global_opts.output_file {
            let (mut file, file_path) = create_output_file(save_flag, k, start_time)?;
            match write_file_header(
                &mut file,
                self.global_opts.k,
                self.global_opts.num_entries,
                &self.command,
                dt,
            ) {
                Ok(()) => {
                    trace!("Wrote file header to {}", file_path);
                }
                Err(_err) => return Err(Error::IOError),
            }
            (Some(file), Some(file_path))
        } else {
            (None, None)
        };
        let command_clone = (self.command).clone();
        let motifs = match self.command {
            Commands::GibbsSampler {
                num_iterations,
                num_runs,
            } => run_gibbs_sampler(&sequences, k, num_runs, num_iterations),
            Commands::MedianString => run_median_string(&sequences, k),
            Commands::Randomized { num_runs } => {
                run_randomized_motif_search(&sequences, k, num_runs)
            }
            Commands::FindMotif { motif, distance } => {
                align_motifs_distance(&sequences, &motif, distance);
                Ok(vec![motif])
            }
        }?;
        let unique_motifs: Vec<String> = unique_motifs(&motifs).into_par_iter().collect();
        let unique_motifs_string = generate_vector_space_delimited(&unique_motifs);
        println!("Unique motifs: {}", unique_motifs_string);
        let consensus_string = generate_consensus_string(&motifs, k)?;
        println!("Consensus string: {}", consensus_string);

        let (best_motif_score, best_motif) = if self.global_opts.align {
            let top_five = align_motifs_multi_threaded(&sequences, &unique_motifs)?;
            println!("Top 5 motifs:");
            for (score, motif) in &top_five {
                println!("{}: {}", score, motif);
            }
            let (best_motif_score, best_motif) = top_five[0].clone();
            align_motifs_distance(&sequences, &consensus_string, 1);
            (Some(best_motif_score), Some(best_motif))
        } else {
            (None, None)
        };
        let dt_end = if let Some(mut file) = file {
            let summary = Summary {
                consensus_string,
                best_motif,
                best_motif_score,
                unique_motifs: unique_motifs_string,
            };
            match output_results_to_file(&mut file, &motifs, &summary, command_clone) {
                Ok(dt_end) => {
                    println!("Results saved to {}", file_path.ok_or(Error::IOError)?);
                    dt_end
                }
                Err(_err) => {
                    error!(
                        "Error writing to file: {}",
                        file_path.ok_or(Error::IOError)?
                    );
                    return Err(Error::IOError);
                }
            }
        } else {
            Utc::now()
        };

        if let Some(duration) = dt_end.signed_duration_since(dt).num_microseconds() {
            info!("Done in {} seconds", duration as f64 / 1_000_000.0);
        }
        Ok(())
    }
}

#[derive(Debug, Args)]
struct GlobalOpts {
    /// file path of the genome
    input_file: String,

    /// how many entries to read
    #[arg(short = 'e', long = "entries")]
    pub num_entries: usize,

    /// motif length
    #[arg(short,value_parser=k_in_range,default_value_t=8)]
    pub k: usize,

    /// alignment
    #[arg(short = 'a', long = "align")]
    align: bool,

    /// save motifs to file
    #[arg(short = 'o', long = "output")]
    output_file: Option<Option<String>>,
}
const K_RANGE: RangeInclusive<usize> = 1..=64;

fn k_in_range(s: &str) -> Result<usize, String> {
    let k: usize = s
        .parse()
        .map_err(|_| format!("`{s}` isn't a valid k value"))?;
    if K_RANGE.contains(&k) {
        Ok(k)
    } else {
        Err(format!(
            "k not in range {}-{}",
            K_RANGE.start(),
            K_RANGE.end()
        ))
    }
}
#[derive(Subcommand, Debug, Clone)]
pub enum Commands {
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
    #[clap(name = "find_motif", about = "Find a motif in a genome")]
    FindMotif {
        motif: String,
        #[arg(short = 'd', long = "distance", default_value_t = 0)]
        distance: u8,
    },
}

pub struct Summary {
    pub consensus_string: String,
    pub unique_motifs: String,
    pub best_motif: Option<String>,
    pub best_motif_score: Option<isize>,
}
