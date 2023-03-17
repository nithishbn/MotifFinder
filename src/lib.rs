
use rand::distributions::WeightedIndex;
use rand::prelude::*;

use std::collections::HashMap;
pub mod median_string;

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




fn scoring_function(motif_matrix: &[String]) -> usize {
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
        // let mut max = 0;
        // for val in count{
        //     if val.1 > max{
        //         max = val.1;
        //     }
        // }

        let max = count.iter().max_by_key(|f| f.1).unwrap().1;

        score += motifs_length - max;
    }
    score
}


pub fn generate_count_matrix(motif_matrix: &[String], k: usize, pseudo: bool) -> Vec<Vec<usize>> {
    // enumerate motif matrix per nucleotide per position
    let mut val = 0;
    if pseudo {
        val = 1;
    }
    let mut count_matrix: Vec<Vec<usize>> = vec![vec![val; k]; 4]; // ACGT = 4
    for i in 0..k {
        // println!("hello");
        for motif in motif_matrix {
            // println!("motif {}",motif);
            let index = match motif.chars().nth(i).unwrap() {
                'A' => Some(0),
                'C' => Some(1),
                'G' => Some(2),
                'T' => Some(3),
                _ => None,
            };
            if let Some(index) = index{
                let count_col: &mut Vec<usize> = count_matrix.get_mut(index).unwrap();
                count_col[i] += 1;
            }
            
            // println!();
        }
    }
    count_matrix
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




pub fn gibbs_sampler(dna: &[String], k: usize, t: usize, n: usize) -> Vec<String> {
    // similar to randomized motif search but at every step we randomly remove one motif from the motifs list
    // we add this back in the form of the profile randomly generated kmer for that profile
    // profile_randomly_generated also adds in a level of randomness based on the profile it generates
    let mut best_motifs = vec![];

    for seq in dna {
        let dna_length = seq.chars().count();
        let start_index = thread_rng().gen_range(0..(dna_length - k + 1));
        // dbg!(dna_length);
        // dbg!(start_index+k);
        best_motifs.push(seq[start_index..start_index + k].to_string());
    }
    // println!("{} {}",best_motifs.len(),t);
    let mut best_score = scoring_function(&best_motifs);
    for _j in 0..n {
        // println!("in loop {_j}");
        let mut motifs = best_motifs.clone();
        let i = thread_rng().gen_range(0..t);
        motifs.remove(i);
        let profile = generate_profile_given_motif_matrix(&best_motifs, true);
        if let Some(motif_i) = profile_randomly_generated_kmer(&dna[i], k, &profile){
            motifs.insert(i, motif_i);
            let test_score = scoring_function(&motifs);
            if test_score < best_score {
                best_motifs = motifs;
                best_score = test_score;
            }
        }
        
    }

    best_motifs
}
pub fn profile_randomly_generated_kmer(text: &str, k: usize, profile: &[Vec<f64>]) -> Option<String> {
    // take in a profile, and for each kmer in text, generate probabilities based on the profile
    // then only output the kmer based on its probability i.e. use a weighted probability
    let n = text.chars().count();
    let mut probabilities: Vec<f64> = vec![];
    let mut kmers = vec![];
    for i in 0..n - k + 1 {
        let slice = &text[i..i + k];
        let kmer = slice.to_string();
        kmers.push(kmer.to_string());
        probabilities.push(generate_probability(&kmer, profile));
    }
    let sum: f64 = probabilities.iter().sum();
    if sum < 0.0{
        return None;
    }
    let adjusted_weights: Vec<f64> = probabilities.iter().map(|f| f / sum).collect();
    // this block of code is taken straight from the rust reference since I am not familiar with the language
    // https://docs.rs/rand/0.7.3/rand/distributions/weighted/struct.WeightedIndex.html
     // similar to random choices from python
    let mut rng = thread_rng();
    if let Ok(dist) = WeightedIndex::new(&adjusted_weights){
        return Some(kmers.get(dist.sample(&mut rng)).unwrap().to_string());
    }
    None
    
}
pub fn iterate_gibbs_sampler(dna: &[String], k: usize, t: usize, iterations: usize, runs: usize) -> Vec<String> {
    // gibbs but iterate
    println!("initializing gibbs sampler");
    let mut motifs = gibbs_sampler(dna, k, t, iterations);
    let mut best_score = scoring_function(&motifs);
    for i in 1..=runs {
        println!("Starting run {i}", );
        let check = gibbs_sampler(dna, k, t, iterations);
        let check_score = scoring_function(&check);
        if check_score < best_score {
            motifs = check;
            best_score = check_score;
        }
    }

    motifs
}
// pub fn main(){
//     let path = "tests/MotifEnumeration/dataset_865379_8.txt".to_string();
//     let (k,d,dna) = motif_enum_load_file(path);
//
// }
