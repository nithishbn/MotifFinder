use std::collections::HashMap;

pub mod alignment;
pub mod gibbs_sampler;
pub mod median_string;
pub mod randomized_motif_search;
pub mod utils;
#[derive(Debug)]
pub enum Error {
    GenericError,
    IOError,
    FileNotFoundError,
    InvalidInputError,
    InvalidNucleotideError,
    InvalidKmerLength,
    InvalidNumberOfRuns,
    InvalidNumberOfIterations,
    InvalidMotifLength,
    NoMotifsFound,
    InvalidSequence,
    InvalidPointerError,
    InvalidNumberofMotifs,
}
pub fn scoring_function(motif_matrix: &[String]) -> usize {
    // given a motif matrix, generate its score by finding the highest count of nucleotide in a given position
    // and subtract that count from the total length of the column

    let mut score = 0;
    let k = motif_matrix.get(0).unwrap().chars().count();
    let motifs_length = motif_matrix.len();
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

        score += motifs_length - max;
    }
    score
}
pub fn generate_profile_given_motif_matrix(
    motif_matrix: &[String],
    pseudo: bool,
) -> Result<Vec<Vec<f64>>, Error> {
    // generate probabilities per column using the count matrix divided by sum of each column
    let k = motif_matrix[0].len();
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
                return Err(Error::InvalidInputError);
            }

            // divide by sum to get the percentage probability
        }
    }
    Ok(profile_matrix)
}
pub fn generate_count_matrix(motif_matrix: &[String], k: usize, pseudo: bool) -> Vec<Vec<usize>> {
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
pub fn generate_probability(kmer: &str, profile: &[Vec<f64>]) -> f64 {
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

pub fn consensus_string(motifs: &[String], k: usize) -> Result<String, Error> {
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
