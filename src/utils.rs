use std::{
    fmt::Display,
    fs::{self, File},
    io::{self, Write},
};

use chrono::{DateTime, Utc};

use crate::{
    command::{Commands, Summary},
    Error,
};

pub fn generate_vector_space_delimited<T: Display>(vec: &[T]) -> String {
    let mut string = "".to_string();
    for val in vec {
        string.push_str(&format!("{val} "));
    }
    let string = string.trim().to_string();
    string
}

pub fn write_file_header(
    file: &mut fs::File,
    k: usize,
    num_entries: usize,
    command: &Commands,
    dt: DateTime<Utc>,
) -> io::Result<()> {
    let version = env!("CARGO_PKG_VERSION");
    writeln!(file, "MotifFinder {}", version)?;
    let command_string = match command {
        Commands::Randomized { .. } => "Randomized",
        Commands::GibbsSampler { .. } => "Gibbs Sampler",
        Commands::MedianString => "Median String",
    };
    writeln!(file, "Command: {}", command_string)?;
    writeln!(file, "k: {}", k)?;
    writeln!(file, "number of entries: {}", num_entries)?;
    match command {
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

pub fn create_output_file(
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
pub fn output_results_to_file(
    file: &mut fs::File,
    motifs: &[String],
    summary: &Summary,
) -> Result<DateTime<Utc>, Error> {
    let Summary {
        consensus_string,
        best_motif_score,
        best_motif,
        unique_motifs,
    } = summary;
    let dt_end = Utc::now();
    writeln!(file, "End time: {}", dt_end.format("%Y-%m-%d %H:%M:%S"))
        .map_err(|_| Error::IOError)?;
    writeln!(file, "Consensus string: {}", consensus_string).map_err(|_| Error::IOError)?;

    writeln!(file, "Unique motifs: {}", unique_motifs).map_err(|_| Error::IOError)?;
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
    Ok(dt_end)
}
fn write_motifs(file: &mut fs::File, motifs: &[String]) -> Result<(), Error> {
    for (i, motif) in motifs.iter().enumerate() {
        let motif = motif.trim();
        writeln!(file, ">motif {}", i + 1).map_err(|_| Error::IOError)?;
        if i == motifs.len() - 1 {
            write!(file, "{}", motif).map_err(|_| Error::IOError)?;
        } else {
            writeln!(file, "{}", motif).map_err(|_| Error::IOError)?;
        }
    }

    Ok(())
}
