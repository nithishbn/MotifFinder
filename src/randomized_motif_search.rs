use rand::{thread_rng, Rng};

use crate::{scoring_function, generate_profile_given_motif_matrix, generate_probability};
/*
I used "better scoring function" which is based on entropy, but I included the integer,
sum-based scoring function as well for reference
*/
pub fn randomized_motif_search(dna: &[String], k: usize) -> Vec<String> {
    let mut best_motifs = vec![];

    for seq in dna {
        let dna_length = seq.chars().count();
        let start_index = thread_rng().gen_range(0..(dna_length - k + 1));
        best_motifs.push(seq[start_index..start_index + k].to_string());
    }

    let mut best_score = scoring_function(&best_motifs);
    loop {
        let profile = generate_profile_given_motif_matrix(&best_motifs, true);
        let motifs = generate_motifs_from_profile(&profile, dna, k);
        let test_score = scoring_function(&motifs);
        if test_score < best_score {
            best_score = test_score;
            best_motifs = motifs;
        } else {
            return best_motifs;
        }
    }
}

pub fn profile_most_probable_kmer(text: &String, k: usize, profile: &[Vec<f64>]) -> String {
    // given a profile, and a DNA string, check all kmers to see which one is the most probable
    let text_len = text.chars().count();
    let mut best_probability_so_far = -1.0;
    let dummy = String::from("");
    let mut best_kmer = dummy;

    for i in 0..(text_len - k + 1) {
        let kmer = text[i..i + k].to_string();
        let kmer_prob = generate_probability(&kmer, profile);
        if kmer_prob > best_probability_so_far {
            best_kmer = kmer;
            best_probability_so_far = kmer_prob;
        }
    }

    best_kmer
}

pub fn generate_motifs_from_profile(
    profile: &Vec<Vec<f64>>,
    dna: &[String],
    k: usize,
) -> Vec<String> {
    let mut motifs: Vec<String> = vec![];
    for seq in dna {
        motifs.push(profile_most_probable_kmer(seq, k, profile));
    }
    motifs
}
pub fn iterate_randomized_motif_search(dna: &[String], k: usize, runs:usize) -> Vec<String> {
    println!("Starting randomized motif search with {} runs", runs);
    let mut motifs = randomized_motif_search(dna, k);
    let mut best_score = scoring_function(&motifs);
    for _i in 0..runs {
        println!("Run: {}", _i+1);
        let check = randomized_motif_search(dna, k);
        let check_score = scoring_function(&check);
        if check_score < best_score {
            motifs = check;
            best_score = check_score;
        }
    }

    motifs
}