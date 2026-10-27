use rust_htslib::bam::{Record, record::Aux};

pub fn signal_lengths(record: &Record) -> Option<Vec<i32>> {
    let moves = match record.aux(b"mv").ok()? {
        Aux::ArrayI8(values) => values,
        _ => return None,
    };
    if moves.len() <= 1 || record.seq_len() == 0 {
        return None;
    }
    let mut signals = vec![0_i32; record.seq_len()];
    let mut base_index = -1_i64;
    for movement in moves.iter().skip(1) {
        if movement != 0 {
            base_index += 1;
            if base_index >= record.seq_len() as i64 {
                break;
            }
            signals[base_index as usize] += 1;
        } else {
            if base_index < 0 {
                continue;
            }
            if base_index >= record.seq_len() as i64 {
                break;
            }
            signals[base_index as usize] += 1;
        }
    }
    if record.flags() & 0x10 != 0 {
        signals.reverse();
    }
    Some(signals)
}

#[cfg(test)]
mod tests {
    use rust_htslib::bam::{self, Read as BamRead};

    use super::*;

    #[test]
    fn parses_valid_and_invalid_move_tables() {
        let root = std::path::Path::new(env!("CARGO_MANIFEST_DIR"));
        let mut reader =
            bam::Reader::from_path(root.join("tests/fixtures/data/toy_dwell.bam")).unwrap();
        let mut records = reader.records();
        let record = records.next().unwrap().unwrap();
        let signals = signal_lengths(&record).unwrap();
        assert_eq!(&signals[..6], &[1, 2, 3, 1, 2, 3]);
        assert!(signals.iter().all(|value| *value > 0));

        let mut reader =
            bam::Reader::from_path(root.join("tests/fixtures/data/invalid_mv.bam")).unwrap();
        let mut records = reader.records();
        let record = records.next().unwrap().unwrap();
        assert!(signal_lengths(&record).is_none());
    }
}
