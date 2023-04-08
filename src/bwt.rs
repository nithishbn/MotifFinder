#[test]
pub fn pssm() {
    use bio::pattern_matching::pssm::DNAMotif;
    use bio::pattern_matching::pssm::Motif;
    let pssm = DNAMotif::from_seqs(
        vec![
            b"AAAA".to_vec(),
            b"AATA".to_vec(),
            b"AAGA".to_vec(),
            b"AAAA".to_vec(),
        ]
        .as_ref(),
        None,
    )
    .unwrap();
    let start_pos = pssm.score(b"CCCCCAATA").unwrap().loc;
    println!("motif found at position {}", start_pos);
}
