use crate::Error;
use rayon::prelude::*;
#[derive(PartialEq, Clone, Eq, Debug)]
enum Pointer {
    DOWN,
    RIGHT,
    DIAGONAL,
    STOP,
    EMPTY,
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
            Pointer::DOWN => {
                v_alignment.insert(0, v.chars().nth(i - 1).unwrap());
                w_alignment.insert(0, '-');
                i -= 1;
            }
            Pointer::RIGHT => {
                w_alignment.insert(0, w.chars().nth(j - 1).unwrap());
                v_alignment.insert(0, '-');
                j -= 1;
            }
            Pointer::DIAGONAL => {
                v_alignment.insert(0, v.chars().nth(i - 1).unwrap());
                w_alignment.insert(0, w.chars().nth(j - 1).unwrap());
                i -= 1;
                j -= 1;
            }
            Pointer::STOP => {
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
    let mut backtrack: Vec<Vec<Pointer>> = vec![vec![Pointer::EMPTY; w_len + 1]; v_len + 1];
    let mut s = vec![vec![0isize; w_len + 1]; v_len + 1];
    for i in 0..=v_len {
        s[i][0] = indel * i as isize;
        backtrack[i][0] = Pointer::DOWN;
    }
    for j in 0..=w_len {
        s[0][j] = indel * j as isize;
        backtrack[0][j] = Pointer::RIGHT;
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
                backtrack[i][j] = Pointer::DOWN;
            } else if s[i][j] == temp[1] {
                backtrack[i][j] = Pointer::RIGHT;
            } else if s[i][j] == temp[2] {
                backtrack[i][j] = Pointer::DIAGONAL;
            } else if s[i][j] == 0 {
                backtrack[i][j] = Pointer::STOP;
            }
        }
    }
    let (row, col) = max_of_matrix(&s);
    let score = s[row][col];
    let alignment = Alignment { score, backtrack };
    Ok((alignment, row, col))
}
