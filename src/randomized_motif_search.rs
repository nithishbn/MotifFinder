use crate::Error;
use crate::{generate_probability, generate_profile_given_motif_matrix, scoring_function};
use indicatif::{ParallelProgressIterator, ProgressBar, ProgressStyle};
use rand::{thread_rng, Rng};
use rayon::prelude::*;
use tracing::trace;
#[tracing::instrument(skip(dna))]
fn randomized_motif_search(dna: &[String], k: usize) -> Result<Vec<String>, Error> {
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
#[tracing::instrument(skip(profile))]
fn profile_most_probable_kmer(text: &str, k: usize, profile: &[Vec<f64>]) -> String {
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

#[tracing::instrument(skip(profile, dna))]
fn generate_motifs_from_profile(profile: &[Vec<f64>], dna: &[String], k: usize) -> Vec<String> {
    let mut motifs: Vec<String> = vec![];
    for seq in dna {
        motifs.push(profile_most_probable_kmer(seq, k, profile));
    }
    motifs
}
#[tracing::instrument(skip(dna))]
pub fn iterate_randomized_motif_search(
    dna: &[String],
    k: usize,
    runs: usize,
) -> Result<Vec<String>, Error> {
    let pb = ProgressBar::new(runs.try_into().map_err(|_| Error::InvalidNumberOfRuns)?);
    trace!("Started randomized motif search");

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

    let mut result: Vec<(usize, Vec<String>)> = (1..=runs)
        .into_par_iter()
        .progress_with(pb.clone())
        .map(|_i| {
            let motifs = randomized_motif_search(dna, k)?;
            let best_score = scoring_function(&motifs);
            Ok((best_score, motifs))
        })
        .collect::<Result<Vec<(usize, Vec<String>)>, Error>>()?;
    result.par_sort_by(|a, b| a.0.cmp(&b.0));
    // dbg!(&result);
    let motifs = result[0].1.clone();
    let best_score = result[0].0;
    pb.finish_with_message(format!("Done! Best score: {best_score}"));
    Ok(motifs)
}
