use crate::Error;
use bio::alignment::pairwise::Aligner;
use bio::alignment::Alignment as BioAlignment;
use bio::pattern_matching::myers::Myers;
use rayon::prelude::*;
#[derive(PartialEq, Clone, Eq, Debug)]
enum Pointer {
    Down,
    Right,
    Diagonal,
    Stop,
    Empty,
}

struct Alignment {
    score: isize,
    backtrack: Vec<Vec<Pointer>>,
}

fn output_backtrack(
    backtrack: &[Vec<Pointer>],
    v: &str,
    w: &str,
    mut i: usize,
    mut j: usize,
) -> Result<(String, String), Error> {
    let mut v_alignment = "".to_string();
    let mut w_alignment = "".to_string();
    while i > 0 || j > 0 {
        match &backtrack[i][j] {
            Pointer::Down => {
                v_alignment.insert(0, v.chars().nth(i - 1).unwrap());
                w_alignment.insert(0, '-');
                i -= 1;
            }
            Pointer::Right => {
                w_alignment.insert(0, w.chars().nth(j - 1).unwrap());
                v_alignment.insert(0, '-');
                j -= 1;
            }
            Pointer::Diagonal => {
                v_alignment.insert(0, v.chars().nth(i - 1).unwrap());
                w_alignment.insert(0, w.chars().nth(j - 1).unwrap());
                i -= 1;
                j -= 1;
            }
            Pointer::Stop => {
                break;
            }
            _ => {
                return Err(Error::InvalidPointerError);
            }
        }
    }
    Ok((v_alignment, w_alignment))
}

pub fn local_alignment(
    v: &str,
    w: &str,
    match_: isize,
    mismatch: isize,
    indel: isize,
) -> Result<(isize, String, String), Error> {
    let (Alignment { score, backtrack }, row, col) =
        local_alignment_score_and_backtrack_matrix(v, w, match_, mismatch, indel)?;
    let (v_alignment, w_alignment) = output_backtrack(&backtrack, v, w, row, col)?;
    Ok((score, v_alignment, w_alignment))
}
fn max_of_matrix(matrix: &[Vec<isize>]) -> (usize, usize) {
    let mut max_so_far = isize::MIN;
    let (mut row, mut col) = (0, 0);
    for i in 0..matrix.len() {
        for j in 0..matrix[0].len() {
            let curr = matrix[i][j];
            if curr > max_so_far {
                max_so_far = curr;
                row = i;
                col = j;
            }
        }
    }
    (row, col)
}

fn local_alignment_score_and_backtrack_matrix(
    v: &str,
    w: &str,
    match_: isize,
    mismatch: isize,
    indel: isize,
) -> Result<(Alignment, usize, usize), Error> {
    let v_len = v.chars().count();
    let w_len = w.chars().count();
    let mut backtrack: Vec<Vec<Pointer>> = vec![vec![Pointer::Empty; w_len + 1]; v_len + 1];
    let mut s = vec![vec![0isize; w_len + 1]; v_len + 1];
    for i in 0..=v_len {
        s[i][0] = indel * i as isize;
        backtrack[i][0] = Pointer::Down;
    }
    for j in 0..=w_len {
        s[0][j] = indel * j as isize;
        backtrack[0][j] = Pointer::Right;
    }
    for i in 1..=v_len {
        for j in 1..=w_len {
            let v_char = v.chars().nth(i - 1);
            let w_char = w.chars().nth(j - 1);
            let matching: isize =
                if (v_char == w_char) && (v_char != Some('N') && w_char != Some('N')) {
                    match_
                } else {
                    mismatch
                };
            let temp = vec![
                s[i - 1][j] + indel,
                s[i][j - 1] + indel,
                s[i - 1][j - 1] + matching,
                0,
            ];
            s[i][j] = *temp.par_iter().max().unwrap();

            if s[i][j] == temp[0] {
                backtrack[i][j] = Pointer::Down;
            } else if s[i][j] == temp[1] {
                backtrack[i][j] = Pointer::Right;
            } else if s[i][j] == temp[2] {
                backtrack[i][j] = Pointer::Diagonal;
            } else if s[i][j] == 0 {
                backtrack[i][j] = Pointer::Stop;
            }
        }
    }
    let (row, col) = max_of_matrix(&s);
    let score = s[row][col];
    let alignment = Alignment { score, backtrack };
    Ok((alignment, row, col))
}

pub fn align_motifs_distance(sequences: &[String], consensus_string: &String) {
    let mut count = 0;
    for (i, sequence) in sequences.iter().enumerate() {
        let pattern = consensus_string.as_bytes();
        let sequence = sequence.as_bytes();
        let mut myers = Myers::<u64>::new(pattern);
        let mut aln = BioAlignment::default();
        let mut matches = myers.find_all(sequence, 2);
        println!("Sequence {}", i + 1);
        while matches.next_alignment(&mut aln) {
            println!(
                "Hit found in range: {}..{} (distance: {})",
                aln.ystart, aln.yend, aln.score
            );
            let y = if aln.ystart >= 2 {
                if aln.yend >= 2 {
                    &sequence[aln.ystart - 2..aln.yend + 2]
                } else {
                    &sequence[aln.ystart - 2..aln.yend]
                }
            } else if aln.yend >= 2 {
                &sequence[aln.ystart..aln.yend + 2]
            } else {
                &sequence[aln.ystart..aln.yend]
            };
            let x = &pattern[aln.xstart..aln.xend];
            let score = |a: u8, b: u8| if a == b { 1i32 } else { -1i32 };
            let mut aligner = Aligner::with_capacity(x.len(), y.len(), -5, -1, &score);
            let alignment = aligner.semiglobal(x, y);
            println!("{}", alignment.pretty(x.as_ref(), y.as_ref()));

            count += 1;
        }
    }
    println!("count: {}", count);
}
