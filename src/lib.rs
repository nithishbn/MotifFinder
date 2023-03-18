use std::collections::HashMap;

pub mod gibbs_sampler;
pub mod median_string;
pub mod randomized_motif_search;

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
            let nuc = motif.chars().nth(i).unwrap();
            *count.entry(nuc).or_insert(0) += 1;
        }

        let max = count.iter().max_by_key(|f| f.1).unwrap().1;

        score += motifs_length - max;
    }
    score
}
pub fn generate_profile_given_motif_matrix(motif_matrix: &[String], pseudo: bool) -> Vec<Vec<f64>> {
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
            let row = profile_matrix.get_mut(j).unwrap(); // get the row associated with the nucleotide at index i
            row[i] = (*count_matrix.get(j).unwrap().get(i).unwrap()) as f64 / sum;
            // divide by sum to get the percentage probability
        }
    }
    profile_matrix
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
            let index = match motif.chars().nth(i).unwrap() {
                'A' => Some(0),
                'C' => Some(1),
                'G' => Some(2),
                'T' => Some(3),
                _ => None,
            };
            if let Some(index) = index {
                let count_col: &mut Vec<usize> = count_matrix.get_mut(index).unwrap();
                count_col[i] += 1;
            }
        }
    }
    count_matrix
}
pub fn generate_probability(kmer: &String, profile: &[Vec<f64>]) -> f64 {
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
        if nuc_index.is_some() {
            // this should always be true but just in case
            let nuc_index = nuc_index.unwrap();
            let current_prob = profile.get(nuc_index).unwrap().get(i).unwrap();
            probability *= current_prob;
        }
    }
    probability
}
