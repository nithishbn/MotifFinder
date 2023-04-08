pub mod alignment;
mod command;
mod gibbs_sampler;
mod median_string;
mod randomized_motif_search;
mod utils;

use alignment::local_alignment;
use gibbs_sampler::iterate_gibbs_sampler;
use indicatif::{MultiProgress, ParallelProgressIterator, ProgressBar, ProgressStyle};
use median_string::median_string;
use randomized_motif_search::iterate_randomized_motif_search;
use rayon::prelude::*;
use std::str;
use std::{
    collections::{HashMap, HashSet},
    fs::File,
};
use tracing::{error, info, trace};

use bio::io::fasta;
#[doc(hidden)]
pub use command::MotifFinder;

#[derive(Debug)]
pub enum Error {
    GenericError,
    IOError,
    FileNotFoundError(String),
    InvalidInputError,
    InvalidNucleotideError,
    InvalidKmerLength,
    InvalidNumberOfRuns,
    InvalidNumberOfIterations,
    InvalidMotifLength,
    NoMotifsFound,
    InvalidSequence,
    InvalidPointerError,
    InvalidNumberMotifs,
}

#[tracing::instrument(skip_all)]
fn scoring_function(motif_matrix: &[String]) -> usize {
    // given a motif matrix, generate its score by finding the highest count of nucleotide in a given position
    // and subtract that count from the total length of the column
    let mut score = 0;
    let k = motif_matrix.get(0).unwrap().chars().count();
    let motifs_length = motif_matrix.len();
    trace!(motifs_length);
    // println!("len {}",motifs_length);
    for i in 0..k {
        let mut count: HashMap<char, usize> = HashMap::new();
        for motif in motif_matrix {
            if let Some(nuc) = motif.chars().nth(i) {
                *count.entry(nuc).or_insert(0) += 1;
            } else {
                continue;
            }
        }

        let max = count.iter().max_by_key(|f| f.1).unwrap().1;
        trace!(max);
        score += motifs_length - max;
    }
    score
}

#[tracing::instrument(skip_all)]
fn generate_profile_given_motif_matrix(
    motif_matrix: &[String],
    pseudo: bool,
) -> Result<Vec<Vec<f64>>, Error> {
    // generate probabilities per column using the count matrix divided by sum of each column
    let k = motif_matrix[0].len();
    trace!(k);
    let count_matrix = generate_count_matrix(motif_matrix, k, pseudo);
    let mut profile_matrix: Vec<Vec<f64>> = vec![vec![0.0; k]; 4];
    let sum = motif_matrix.len() as f64;
    // iterating over each position
    for i in 0..k {
        // iterating over each nucleotide base
        for j in 0..4 {
            // print_vector_space_delimited(row.clone());
            // get the row associated with the nucleotide at index i
            if let Some(row) = profile_matrix.get_mut(j) {
                row[i] = (*count_matrix.get(j).unwrap().get(i).unwrap()) as f64 / sum;
            } else {
                error!("Invalid index for profile matrix");
                return Err(Error::InvalidInputError);
            }

            // divide by sum to get the percentage probability
        }
    }
    Ok(profile_matrix)
}

#[tracing::instrument(skip_all)]
fn generate_count_matrix(motif_matrix: &[String], k: usize, pseudo: bool) -> Vec<Vec<usize>> {
    // enumerate motif matrix per nucleotide per position
    let mut val = 0;
    if pseudo {
        val = 1;
    }
    let mut count_matrix: Vec<Vec<usize>> = vec![vec![val; k]; 4]; // ACGT = 4
    for i in 0..k {
        for motif in motif_matrix {
            if let Some(index) = match motif.chars().nth(i) {
                Some('A') => Some(0),
                Some('C') => Some(1),
                Some('G') => Some(2),
                Some('T') => Some(3),
                _ => None,
            } {
                if let Some(count_col) = count_matrix.get_mut(index) {
                    count_col[i] += 1;
                }
            }
        }
    }
    count_matrix
}

#[tracing::instrument(skip_all)]
fn generate_probability(kmer: &str, profile: &[Vec<f64>]) -> f64 {
    // given a kmer and a profile, generate its probability
    let mut probability = 1.0;
    for (i, nuc) in kmer.chars().enumerate() {
        let nuc_index = match nuc {
            'A' => Some(0),
            'C' => Some(1),
            'G' => Some(2),
            'T' => Some(3),
            _ => None,
        };
        if let Some(nuc_index) = nuc_index {
            // this should always be true but just in case
            let current_prob = profile.get(nuc_index).unwrap().get(i).unwrap();
            probability *= current_prob;
        }
    }
    probability
}

#[tracing::instrument(skip_all)]
fn consensus_string(motifs: &[String], k: usize) -> Result<String, Error> {
    let mut consensus = String::new();
    let count_matrix = generate_count_matrix(motifs, k, true);
    for i in 0..k {
        let mut max = 0;
        let mut max_index = 0;
        for j in 0..4 {
            let count = count_matrix
                .get(j)
                .and_then(|row| row.get(i))
                .ok_or(Error::InvalidNucleotideError)?;
            if count > &max {
                max = *count;
                max_index = j;
            }
        }
        let nuc = match max_index {
            0 => 'A',
            1 => 'C',
            2 => 'G',
            3 => 'T',
            _ => return Err(Error::InvalidNucleotideError),
        };
        consensus.push(nuc);
    }
    Ok(consensus)
}

#[tracing::instrument]
pub fn align_motifs_multi_threaded(
    sequences: &[String],
    motifs: &[String],
) -> Result<Vec<(isize, String)>, Error> {
    let motifs_len = motifs.len();
    let sequences_len = sequences.len();
    let pb = ProgressBar::new(
        motifs_len
            .try_into()
            .map_err(|_| Error::InvalidNumberMotifs)?,
    );
    let sty =
        ProgressStyle::with_template(&format!("{{prefix:.bold}}▕{{bar:.{}}}▏{{msg}} ", "9.on_0"))
            .unwrap()
            .progress_chars("█▉▊▋▌▍▎▏  ");
    pb.set_style(
        ProgressStyle::with_template(
            "[{elapsed_precise}] {spinner:.9.on_0} {bar:50.9.on_0} {pos:>2}/{len:2} {msg} ({eta})",
        )
        .unwrap(),
    );
    pb.reset_eta();

    let m = MultiProgress::new();

    let total_pb = m.add(pb);
    total_pb.println(format!(
        "Aligning {} unique motifs to {} sequences",
        motifs_len,
        sequences.len()
    ));

    let mut top_five: Vec<(isize, String)> = motifs
        .par_iter()
        .progress_with(total_pb.clone())
        .map(|motif| {
            let inner = m.add(ProgressBar::new(sequences_len.try_into().unwrap()));
            inner.set_style(sty.clone());
            inner.set_prefix(motif.to_string());
            let mut total_score = 0;
            let mut highest_score = 0;
            let mut best_motif = String::from("");
            for sequence in sequences.iter() {
                let (score, _v_align, w_align) = local_alignment(sequence, motif, 1, -10, -100)?;
                if score > highest_score {
                    highest_score = score;
                    best_motif = w_align;
                }
                total_score += score;
                inner.inc(1);
            }
            inner.finish_and_clear();
            Ok((total_score, best_motif))
        })
        .collect::<Result<Vec<(isize, String)>, Error>>()?;

    total_pb.finish_with_message("Done!");
    top_five.par_sort_by(|a, b| b.0.cmp(&a.0));
    top_five.dedup();
    top_five.truncate(5);
    Ok(top_five.to_vec())
}

#[tracing::instrument]
pub fn load_data(path_to_file: &str, num_entries: usize) -> Result<Vec<String>, Error> {
    info!("Loading data from '{}'...", path_to_file);
    let mut sequences = vec![];
    let file = match File::open(path_to_file) {
        Ok(file) => file,
        Err(_) => return Err(Error::FileNotFoundError(path_to_file.to_string())),
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
    info!("Done loading data: {} entries", sequences.len());
    Ok(sequences)
}

#[tracing::instrument(skip(sequences))]
pub fn run_gibbs_sampler(
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

#[tracing::instrument(skip(sequences))]
pub fn run_median_string(sequences: &[String], k: usize) -> Result<Vec<String>, Error> {
    let median_string = median_string(k, sequences)?;
    info!("Median string: {}", median_string);
    let vec = vec![median_string];
    Ok(vec)
}

#[tracing::instrument(skip(sequences))]
pub fn run_randomized_motif_search(
    sequences: &[String],
    k: usize,
    num_runs: usize,
) -> Result<Vec<String>, Error> {
    if num_runs == 0 {
        return Err(Error::InvalidNumberOfRuns);
    }
    iterate_randomized_motif_search(sequences, k, num_runs)
}

#[tracing::instrument(skip(motifs))]
pub fn generate_consensus_string(motifs: &[String], k: usize) -> Result<String, Error> {
    if motifs.is_empty() {
        return Err(Error::NoMotifsFound);
    } else if motifs.len() == 1 {
        return Ok(motifs[0].clone());
    }
    consensus_string(motifs, k)
}

#[tracing::instrument(skip(motifs))]
pub fn unique_motifs(motifs: &[String]) -> HashSet<String> {
    motifs.into_par_iter().cloned().collect::<HashSet<String>>()
}

#[cfg(test)]
mod test {
    use crate::align_motifs_multi_threaded;

    #[test]
    pub fn test_load_data() {
        let sequences = super::load_data("promoters.fasta", 5).unwrap();
        assert_eq!(sequences.len(), 4);
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

    #[test]
    pub fn test_entries_less_than_five() {
        let sequences = super::load_data("promoters.fasta", 4).unwrap();
        let motifs = super::run_randomized_motif_search(&sequences, 8, 20).unwrap();
        let top_five = align_motifs_multi_threaded(sequences, motifs).unwrap();
        assert!(top_five.len() <= 4);
        let sequences = super::load_data("promoters.fasta", 2).unwrap();
        assert_eq!(sequences.len(), 2);
        let motifs = super::run_randomized_motif_search(&sequences, 8, 20).unwrap();
        assert_eq!(motifs.len(), 2);
        let top_five = align_motifs_multi_threaded(sequences, motifs).unwrap();
        assert!(top_five.len() <= 2);
    }
}
