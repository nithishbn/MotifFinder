use rand::distributions::WeightedIndex;
use rand::prelude::*;

use crate::Error;
use crate::{generate_probability, generate_profile_given_motif_matrix, scoring_function};
pub fn gibbs_sampler(dna: &[String], k: usize, t: usize, n: usize) -> Result<Vec<String>, Error> {
    // similar to randomized motif search but at every step we randomly remove one motif from the motifs list
    // we add this back in the form of the profile randomly generated kmer for that profile
    // profile_randomly_generated also adds in a level of randomness based on the profile it generates
    let mut best_motifs = vec![];

    for seq in dna {
        let dna_length = seq.chars().count();
        let start_index = thread_rng().gen_range(0..(dna_length - k + 1));
        // dbg!(dna_length);
        // dbg!(start_index+k);
        if k > dna_length {
            continue;
        }
        best_motifs.push(seq[start_index..start_index + k].to_string());
    }
    // println!("{} {}",best_motifs.len(),t);
    let mut best_score = scoring_function(&best_motifs);
    for _j in 0..n {
        // println!("in loop {_j}");
        let mut motifs = best_motifs.clone();
        let i = thread_rng().gen_range(0..t);
        motifs.remove(i);
        let profile = generate_profile_given_motif_matrix(&best_motifs, true)?;
        if let Some(motif_i) = profile_randomly_generated_kmer(&dna[i], k, &profile) {
            motifs.insert(i, motif_i);
            let test_score = scoring_function(&motifs);
            if test_score < best_score {
                best_motifs = motifs;
                best_score = test_score;
            }
        }
    }

    Ok(best_motifs)
}
pub fn profile_randomly_generated_kmer(
    text: &str,
    k: usize,
    profile: &[Vec<f64>],
) -> Option<String> {
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
    if sum < 0.0 {
        return None;
    }
    let adjusted_weights: Vec<f64> = probabilities.iter().map(|f| f / sum).collect();
    // this block of code is taken straight from the rust reference since I am not familiar with the language
    // https://docs.rs/rand/0.7.3/rand/distributions/weighted/struct.WeightedIndex.html
    // similar to random choices from python
    let mut rng = thread_rng();
    if let Ok(dist) = WeightedIndex::new(&adjusted_weights) {
        return Some(kmers.get(dist.sample(&mut rng)).unwrap().to_string());
    }
    None
}
pub fn iterate_gibbs_sampler(
    dna: &[String],
    k: usize,
    t: usize,
    iterations: usize,
    runs: usize,
) -> Result<Vec<String>, Error> {
    // gibbs but iterate
    println!("initializing gibbs sampler");
    let mut motifs = gibbs_sampler(dna, k, t, iterations)?;
    let mut best_score = scoring_function(&motifs);
    for i in 1..=runs {
        println!("Starting run {i}",);
        let check = gibbs_sampler(dna, k, t, iterations)?;
        let check_score = scoring_function(&check);
        if check_score < best_score {
            motifs = check;
            best_score = check_score;
        }
    }

    Ok(motifs)
}
