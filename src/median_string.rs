use std::collections::HashSet;

use crate::Error;
use indicatif::{ProgressBar, ProgressStyle};
use rayon::prelude::*;
use tracing::{trace, warn};
#[tracing::instrument(skip_all)]
fn hamming_distance(string1: &str, string2: &str) -> usize {
    trace!("Hamming distance between {} and {}", string1, string2);
    // scan linearly across both strings to find how many differences they have between each other
    let length = string1.chars().count();
    let mut distance = 0;
    let string1_vec: Vec<char> = string1.chars().collect();
    let string2_vec: Vec<char> = string2.chars().collect();
    for i in 0..length {
        if string1_vec.get(i) != string2_vec.get(i) {
            distance += 1;
        }
    }
    distance
}
#[tracing::instrument]
fn neighbors(pattern: String, d: usize) -> HashSet<String> {
    trace!("Generating neighbors of {} with distance {}", pattern, d);
    // generate all neighbors of length |pattern| by modifying at most d nucleotides
    if d == 0 {
        let pattern_set: HashSet<String> = vec![pattern].into_par_iter().collect();
        return pattern_set;
    }
    if pattern.len() == 1 {
        let base_case: HashSet<String> = vec![
            "A".to_string(),
            "C".to_string(),
            "G".to_string(),
            "T".to_string(),
        ]
        .into_par_iter()
        .collect();
        return base_case;
    }
    let mut neighborhood: HashSet<String> = HashSet::new();
    let suffix_neighbors = neighbors(pattern[1..].to_string(), d);
    for text in suffix_neighbors.iter() {
        if hamming_distance(&pattern[1..], text) < d {
            // this line is messy I apologize
            for nuc in vec![
                "A".to_string(),
                "C".to_string(),
                "G".to_string(),
                "T".to_string(),
            ] {
                neighborhood.insert(nuc + text);
            }
        } else {
            neighborhood.insert(pattern.chars().next().unwrap().to_string() + text);
        }
    }
    neighborhood
}
#[tracing::instrument(skip(dna))]
pub fn median_string(k: usize, dna: &[String]) -> Result<String, Error> {
    trace!("Finding median string of length {} in {:?}", k, dna);
    let mut distance = usize::MAX;
    let dummy_string = "A".repeat(k);
    let patterns = neighbors(dummy_string, k);
    let mut median = String::from("");
    let len = patterns.len();
    let pb = ProgressBar::new(len.try_into().map_err(|_| Error::GenericError)?);
    pb.println(format!(
        "Starting median string search for {k}-mers in {len} sequences",
    ));
    let sty = ProgressStyle::with_template(
        "[{elapsed_precise}] {spinner:.green} {bar:40.cyan/blue} {pos:>7}/{len:7} {msg} ({eta})",
    )
    .unwrap();
    pb.set_style(sty);
    pb.reset_eta();
    pb.set_message("Initializing");
    for (_i, pattern) in patterns.iter().enumerate() {
        pb.set_message(format!("Checking pattern: {pattern}"));
        pb.inc(1);
        let pattern_distance = distance_between_pattern_and_strings(pattern, dna)?;
        if distance > pattern_distance {
            distance = pattern_distance;
            median = pattern.to_string();
        }
    }
    pb.finish_with_message("Done!");

    Ok(median)
}

#[tracing::instrument(skip(dna))]
fn distance_between_pattern_and_strings(pattern: &str, dna: &[String]) -> Result<usize, Error> {
    let k = pattern.chars().count();
    let mut distance: usize = 0;
    for (_i, seq) in dna.iter().enumerate() {
        let mut hammingdist = usize::MAX;
        let seq_len = seq.chars().count();
        if k > seq_len {
            warn!("Pattern length is greater than sequence length, skipping sequence");
            continue;
        }
        for i in 0..seq_len - k + 1 {
            let kmer = &seq[i..i + k].to_string();
            let new_hamming = hamming_distance(pattern, kmer);
            if hammingdist > new_hamming {
                hammingdist = new_hamming;
            }
        }
        distance += hammingdist;
    }
    Ok(distance)
}
