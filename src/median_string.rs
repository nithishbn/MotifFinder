use rayon::prelude::*;
use std::collections::HashSet;

pub fn hamming_distance(string1: &str, string2: &str) -> usize {
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
pub fn neighbors(pattern: String, d: usize) -> HashSet<String> {
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
pub fn median_string(k: usize, dna: &[String]) -> String {
    let mut distance = usize::MAX;
    let dummy_string = "A".repeat(k);
    let patterns = neighbors(dummy_string, k);
    let mut median = String::from("");
    let len = patterns.len();
    // let mut file = fs::File::create(format!("median-string-{timestamp}-{k}-checkpoint.txt")).expect("Unable to create file");
    for (i, pattern) in patterns.iter().enumerate() {
        println!(
            "processing pattern {i} of {len} in median_string",
            i = i + 1
        );
        let pattern_distance = distance_between_pattern_and_strings(pattern, dna);
        if distance > pattern_distance {
            distance = pattern_distance;
            median = pattern.to_string();
        }
        // if i % 1000 == 0{
        //     let dt = Utc::now();
        //     write!(file, "{dt} - processing pattern {i} of {len} in median_string\n",i=i+1).expect("Unable to write data");
        // }
    }

    median
}

pub fn distance_between_pattern_and_strings(pattern: &str, dna: &[String]) -> usize {
    let k = pattern.chars().count();
    let mut distance: usize = 0;
    for (i, seq) in dna.iter().enumerate() {
        println!(
            "processing distance of seq {i} of {len}",
            i = i + 1,
            len = dna.len()
        );
        let mut hammingdist = usize::MAX;
        let seq_len = seq.chars().count();

        for i in 0..seq_len - k + 1 {
            let kmer = &seq[i..i + k].to_string();
            let new_hamming = hamming_distance(pattern, kmer);
            if hammingdist > new_hamming {
                hammingdist = new_hamming;
            }
        }
        distance += hammingdist;
    }
    distance
}
