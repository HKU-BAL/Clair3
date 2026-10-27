/*
The MIT License

Copyright (c) 2012-2015 Boston College.

Permission is hereby granted, free of charge, to any person obtaining
a copy of this software and associated documentation files (the
"Software"), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the Software, and to
permit persons to whom the Software is furnished to do so, subject to
the following conditions:

The above copyright notice and this permission notice shall be
included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

const MATCH: i32 = 4;
const MISMATCH: i32 = -6;
const GAP_OPEN: i32 = 8;
const GAP_EXTEND: i32 = 2;

#[derive(Debug, Default, Clone, PartialEq, Eq)]
pub struct Alignment {
    pub score: i32,
    pub reference_start: i32,
    pub reference_end: i32,
    pub query_start: i32,
    pub query_end: i32,
    pub cigar: String,
}

#[derive(Clone, Copy)]
struct Score {
    score: i32,
    reference_end: i32,
    query_end: i32,
}

fn base_code(base: u8) -> usize {
    match base.to_ascii_uppercase() {
        b'A' => 0,
        b'C' => 1,
        b'G' => 2,
        b'T' | b'U' => 3,
        _ => 4,
    }
}

fn substitution(query: u8, reference: u8) -> i32 {
    if base_code(query) == base_code(reference) && base_code(query) < 4 {
        MATCH
    } else {
        MISMATCH
    }
}

fn best_score(query: &[u8], reference: &[u8], reverse: bool, terminate: i32) -> Score {
    let mut previous = vec![0_i32; query.len() + 1];
    let mut current = vec![0_i32; query.len() + 1];
    let mut deletions = vec![0_i32; query.len() + 1];
    let mut maximum = 0;
    let mut reference_end = 0;
    let mut query_end = query.len() as i32 - 1;
    for index in 0..reference.len() {
        let reference_position = if reverse {
            reference.len() - index - 1
        } else {
            index
        };
        current[0] = 0;
        let mut insertion = 0;
        let mut column_maximum = 0;
        for query_position in 0..query.len() {
            let deletion = (previous[query_position + 1] - GAP_OPEN)
                .max(deletions[query_position + 1] - GAP_EXTEND)
                .max(0);
            insertion = (current[query_position] - GAP_OPEN)
                .max(insertion - GAP_EXTEND)
                .max(0);
            let score = (previous[query_position]
                + substitution(query[query_position], reference[reference_position]))
            .max(deletion)
            .max(insertion)
            .max(0);
            current[query_position + 1] = score;
            deletions[query_position + 1] = deletion;
            column_maximum = column_maximum.max(score);
        }
        if column_maximum > maximum {
            maximum = column_maximum;
            reference_end = reference_position as i32;
            query_end = (0..query.len())
                .find(|position| current[position + 1] == maximum)
                .unwrap_or(query.len() - 1) as i32;
        }
        if column_maximum == terminate {
            break;
        }
        std::mem::swap(&mut previous, &mut current);
    }
    Score {
        score: maximum,
        reference_end,
        query_end,
    }
}

#[allow(clippy::needless_range_loop)]
fn banded_cigar(query: &[u8], reference: &[u8], score: i32) -> Option<Vec<(usize, u8)>> {
    let mut band_width = query.len().abs_diff(reference.len()) + 1;
    let (width, direction) = loop {
        let width = band_width * 2 + 3;
        let direction_width = band_width * 2 + 1;
        let mut h_previous = vec![0_i32; width];
        let mut h_current = vec![0_i32; width];
        let mut e = vec![0_i32; width];
        let mut direction = vec![0_u8; query.len() * direction_width * 3];
        let mut maximum = 0;
        for query_position in 0..query.len() {
            let current_start = (query_position as isize - band_width as isize).max(0) as usize;
            let previous_start =
                (query_position as isize - 1 - band_width as isize).max(0) as usize;
            let begin = current_start;
            let end = (reference.len() - 1).min(query_position + band_width);
            let edge = (end + 1).min(width - 1);
            let mut f = 0;
            h_previous[0] = 0;
            e[0] = 0;
            h_previous[edge] = 0;
            e[edge] = 0;
            h_current[0] = 0;
            let mut last_u = 0;
            for reference_position in begin..=end {
                let u = reference_position - current_start + 1;
                let previous_u = reference_position - previous_start + 1;
                let left_u = reference_position - current_start;
                let diagonal_u = reference_position - previous_start;
                last_u = u;
                let offset =
                    query_position * direction_width * 3 + (reference_position - current_start) * 3;
                let open = if query_position == 0 {
                    -GAP_OPEN
                } else {
                    h_previous[previous_u] - GAP_OPEN
                };
                let extend = if query_position == 0 {
                    -GAP_EXTEND
                } else {
                    e[previous_u] - GAP_EXTEND
                };
                e[u] = open.max(extend);
                direction[offset] = if open > extend { 3 } else { 2 };
                let open = h_current[left_u] - GAP_OPEN;
                let extend = f - GAP_EXTEND;
                f = open.max(extend);
                direction[offset + 1] = if open > extend { 5 } else { 4 };
                let e_score = e[u].max(0);
                let f_score = f.max(0);
                let gap_score = e_score.max(f_score);
                let diagonal_score = h_previous[diagonal_u]
                    + substitution(query[query_position], reference[reference_position]);
                h_current[u] = gap_score.max(diagonal_score);
                maximum = maximum.max(h_current[u]);
                direction[offset + 2] = if gap_score <= diagonal_score {
                    1
                } else if e_score > f_score {
                    direction[offset]
                } else {
                    direction[offset + 1]
                };
            }
            h_previous[1..=last_u].copy_from_slice(&h_current[1..=last_u]);
        }
        if maximum >= score {
            break (direction_width, direction);
        }
        if band_width > query.len() + reference.len() {
            return None;
        }
        band_width *= 2;
    };
    let mut reversed = Vec::new();
    let mut query_position = query.len() - 1;
    let mut reference_position = reference.len() as isize - 1;
    let mut length = 0;
    let mut operation = b'M';
    let mut previous_operation = b'M';
    let mut state = 2_usize;
    while query_position > 0 {
        let start = (query_position as isize - (width as isize - 1) / 2).max(0) as usize;
        if reference_position < start as isize || reference_position >= reference.len() as isize {
            return None;
        }
        let direction_index =
            query_position * width * 3 + (reference_position as usize - start) * 3 + state;
        match direction[direction_index] {
            1 => {
                query_position -= 1;
                reference_position -= 1;
                state = 2;
                operation = b'M';
            }
            2 | 3 => {
                query_position -= 1;
                state = if direction[direction_index] == 2 {
                    0
                } else {
                    2
                };
                operation = b'I';
            }
            4 | 5 => {
                reference_position -= 1;
                state = if direction[direction_index] == 4 {
                    1
                } else {
                    2
                };
                operation = b'D';
            }
            _ => return None,
        }
        if operation == previous_operation {
            length += 1;
        } else {
            reversed.push((length, previous_operation));
            previous_operation = operation;
            length = 1;
        }
    }
    if operation == b'M' {
        reversed.push((length + 1, b'M'));
    } else {
        reversed.push((length, operation));
        reversed.push((1, b'M'));
    }
    reversed.reverse();
    Some(reversed)
}

fn push_cigar(cigar: &mut Vec<(usize, u8)>, length: usize, operation: u8) {
    if length == 0 {
        return;
    }
    if let Some(last) = cigar.last_mut()
        && last.1 == operation
    {
        last.0 += length;
        return;
    }
    cigar.push((length, operation));
}

pub fn smith_waterman(query: &[u8], reference: &[u8]) -> Alignment {
    if query.is_empty() || reference.is_empty() {
        return Alignment::default();
    }
    let forward = best_score(query, reference, false, -1);
    if forward.score == 0 {
        return Alignment {
            score: 0,
            reference_start: 0,
            reference_end: forward.reference_end,
            query_start: 0,
            query_end: forward.query_end,
            cigar: String::new(),
        };
    }
    let mut reverse_query = query[..=forward.query_end as usize].to_vec();
    reverse_query.reverse();
    let reverse = best_score(
        &reverse_query,
        &reference[..=forward.reference_end as usize],
        true,
        forward.score,
    );
    let reference_start = reverse.reference_end as usize;
    let query_start = forward.query_end - reverse.query_end;
    let raw = banded_cigar(
        &query[query_start as usize..=forward.query_end as usize],
        &reference[reference_start..=forward.reference_end as usize],
        forward.score,
    );
    let Some(raw) = raw else {
        return Alignment::default();
    };
    let mut cigar = Vec::new();
    push_cigar(&mut cigar, query_start as usize, b'S');
    let mut query_position = query_start as usize;
    let mut reference_position = reference_start;
    for (length, operation) in raw {
        match operation {
            b'M' => {
                for _ in 0..length {
                    let operation = if base_code(query[query_position])
                        == base_code(reference[reference_position])
                    {
                        b'='
                    } else {
                        b'X'
                    };
                    push_cigar(&mut cigar, 1, operation);
                    query_position += 1;
                    reference_position += 1;
                }
            }
            b'I' => {
                push_cigar(&mut cigar, length, operation);
                query_position += length;
            }
            b'D' => {
                push_cigar(&mut cigar, length, operation);
                reference_position += length;
            }
            _ => unreachable!(),
        }
    }
    push_cigar(
        &mut cigar,
        query.len() - forward.query_end as usize - 1,
        b'S',
    );
    let cigar = cigar
        .into_iter()
        .map(|(length, operation)| format!("{length}{}", operation as char))
        .collect();
    Alignment {
        score: forward.score,
        reference_start: reference_start as i32,
        reference_end: forward.reference_end,
        query_start,
        query_end: forward.query_end,
        cigar,
    }
}
