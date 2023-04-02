use crate::Error;
use crate::{generate_probability, generate_profile_given_motif_matrix, scoring_function};
use indicatif::{ProgressBar, ProgressStyle};
use rand::{thread_rng, Rng};
/*
I used "better scoring function" which is based on entropy, but I included the integer,
sum-based scoring function as well for reference
*/
pub fn randomized_motif_search(dna: &[String], k: usize) -> Result<Vec<String>, Error> {
    let mut best_motifs = vec![];

    for seq in dna {
        let dna_length = seq.chars().count();
        let start_index = thread_rng().gen_range(0..(dna_length - k + 1));
        if k > dna_length {
            continue;
        }
        best_motifs.push(seq[start_index..start_index + k].to_string());
    }

    let mut best_score = scoring_function(&best_motifs);
    loop {
        let profile = generate_profile_given_motif_matrix(&best_motifs, true)?;
        let motifs = generate_motifs_from_profile(&profile, dna, k);
        let test_score = scoring_function(&motifs);
        if test_score < best_score {
            best_score = test_score;
            best_motifs = motifs;
        } else {
            return Ok(best_motifs);
        }
    }
}

pub fn profile_most_probable_kmer(text: &str, k: usize, profile: &[Vec<f64>]) -> String {
    // given a profile, and a DNA string, check all kmers to see which one is the most probable
    let text_len = text.chars().count();
    let mut best_probability_so_far = -1.0;
    let dummy = "";
    let mut best_kmer = dummy;

    for i in 0..(text_len - k + 1) {
        if k > text_len {
            continue;
        }
        let kmer = &text[i..i + k];
        let kmer_prob = generate_probability(kmer, profile);
        if kmer_prob > best_probability_so_far {
            best_kmer = kmer;
            best_probability_so_far = kmer_prob;
        }
    }

    best_kmer.to_owned()
}

pub fn generate_motifs_from_profile(profile: &[Vec<f64>], dna: &[String], k: usize) -> Vec<String> {
    let mut motifs: Vec<String> = vec![];
    for seq in dna {
        motifs.push(profile_most_probable_kmer(seq, k, profile));
    }
    motifs
}
pub fn iterate_randomized_motif_search(
    dna: &[String],
    k: usize,
    runs: usize,
) -> Result<Vec<String>, Error> {
    let pb = ProgressBar::new(runs.try_into().map_err(|_| Error::InvalidNumberOfRuns)?);
    pb.println(format!(
        "Starting randomized motif search with {} runs",
        runs
    ));
    let sty = ProgressStyle::with_template(
        "[{elapsed_precise}] {spinner:.green} {bar:40.cyan/blue} {pos:>7}/{len:7} {msg} ({eta})",
    )
    .unwrap();
    pb.set_style(sty);
    pb.reset_eta();
    pb.set_message("Initializing");
    let mut motifs = randomized_motif_search(dna, k)?;
    let mut best_score = scoring_function(&motifs);

    for _i in 1..=runs {
        pb.set_message(format!("Score so far {best_score}"));
        let check = randomized_motif_search(dna, k)?;
        let check_score = scoring_function(&check);
        pb.inc(1);
        if check_score < best_score {
            motifs = check;
            best_score = check_score;
        }
    }
    pb.finish_with_message(format!("Done! Best score: {best_score}"));
    Ok(motifs)
}
