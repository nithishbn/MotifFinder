use crate::Error;
use crate::{generate_probability, generate_profile_given_motif_matrix, scoring_function};
use indicatif::{ParallelProgressIterator, ProgressBar, ProgressStyle};
use rand::distributions::WeightedIndex;
use rand::prelude::*;
use rayon::prelude::*;
use tracing::{info, trace};
#[tracing::instrument(skip(dna))]
fn gibbs_sampler(dna: &[String], k: usize, t: usize, n: usize) -> Result<Vec<String>, Error> {
    // similar to randomized motif search but at every step we randomly remove one motif from the motifs list
    // we add this back in the form of the profile randomly generated kmer for that profile
    // profile_randomly_generated also adds in a level of randomness based on the profile it generates
    let mut best_motifs = vec![];

    for seq in dna {
        let dna_length = seq.chars().count();
        let start_index = thread_rng().gen_range(0..(dna_length - k + 1));
        if k > dna_length {
            continue;
        }
        best_motifs.push(seq[start_index..start_index + k].to_string());
    }
    // println!("{} {}",best_motifs.len(),t);
    let mut best_score = scoring_function(&best_motifs);
    for _j in 0..n {
        trace!("Gibbs Sampler iteration: {}", _j);
        let mut motifs = best_motifs.clone();
        let i = thread_rng().gen_range(0..t);
        trace!("Removing {}th motif", i);
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
#[tracing::instrument(skip_all)]
fn profile_randomly_generated_kmer(text: &str, k: usize, profile: &[Vec<f64>]) -> Option<String> {
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
    let sum: f64 = probabilities.par_iter().sum();
    if sum < 0.0 {
        return None;
    }
    let adjusted_weights: Vec<f64> = probabilities.par_iter().map(|f| f / sum).collect();
    // this block of code is taken straight from the rust reference since I am not familiar with the language
    // https://docs.rs/rand/0.7.3/rand/distributions/weighted/struct.WeightedIndex.html
    // similar to random choices from python
    let mut rng = thread_rng();
    if let Ok(dist) = WeightedIndex::new(&adjusted_weights) {
        return Some(kmers.get(dist.sample(&mut rng)).unwrap().to_string());
    }
    None
}
#[tracing::instrument(skip_all)]
pub fn iterate_gibbs_sampler(
    dna: &[String],
    k: usize,
    t: usize,
    iterations: usize,
    runs: usize,
) -> Result<Vec<String>, Error> {
    // gibbs but iterate
    info!("Initializing Gibbs Sampler");
    let pb = ProgressBar::new(runs.try_into().map_err(|_| Error::InvalidNumberOfRuns)?);
    let sty = ProgressStyle::with_template(
        "[{elapsed_precise}] {spinner:.green} {bar:40.cyan/blue} {pos:>7}/{len:7} {msg} ({eta})",
    )
    .unwrap();
    pb.set_style(sty);
    pb.reset_eta();
    pb.println(format!(
        "Starting Gibbs Sampler with {runs} runs and {iterations} iterations"
    ));

    let mut result: Vec<(usize, Vec<String>)> = (1..=runs)
        .into_par_iter()
        .progress_with(pb.clone())
        .map(|_i| {
            let motifs = gibbs_sampler(dna, k, t, iterations)?;
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
