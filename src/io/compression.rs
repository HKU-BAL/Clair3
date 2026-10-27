use std::{
    fs::File,
    io::{self, BufRead, BufReader},
    path::Path,
};

use anyhow::Result;
use flate2::read::MultiGzDecoder;

pub fn open_text(path: &Path) -> Result<Box<dyn BufRead>> {
    if path == Path::new("PIPE") || path == Path::new("-") {
        return Ok(Box::new(BufReader::new(io::stdin())));
    }
    let file = File::open(path)?;
    let reader = BufReader::new(file);
    match path.extension().and_then(|value| value.to_str()) {
        Some("gz") | Some("bgz") | Some("bgzf") => {
            Ok(Box::new(BufReader::new(MultiGzDecoder::new(reader))))
        }
        Some("zst") | Some("zstd") => Ok(Box::new(BufReader::new(
            zstd::stream::read::Decoder::new(reader)?,
        ))),
        _ => Ok(Box::new(reader)),
    }
}

pub fn read_text_lines(path: &Path) -> Result<impl Iterator<Item = Result<String>>> {
    let reader = open_text(path)?;
    Ok(reader.lines().map(|line| line.map_err(anyhow::Error::from)))
}

#[cfg(test)]
mod tests {
    use std::io::Write;

    use super::*;

    #[test]
    fn reads_plain_gzip_and_zstd_text() {
        let directory = tempfile::tempdir().unwrap();
        let plain = directory.path().join("input.txt");
        std::fs::write(&plain, b"one\ntwo\n").unwrap();
        let mut lines = read_text_lines(&plain).unwrap();
        assert_eq!(lines.next().unwrap().unwrap(), "one");
        assert_eq!(lines.next().unwrap().unwrap(), "two");

        let gzip = directory.path().join("input.txt.gz");
        let mut encoder = flate2::write::GzEncoder::new(
            File::create(&gzip).unwrap(),
            flate2::Compression::default(),
        );
        encoder.write_all(b"three\nfour\n").unwrap();
        encoder.finish().unwrap();
        let mut lines = read_text_lines(&gzip).unwrap();
        assert_eq!(lines.next().unwrap().unwrap(), "three");
        assert_eq!(lines.next().unwrap().unwrap(), "four");

        let zstd_path = directory.path().join("input.txt.zst");
        let mut encoder =
            zstd::stream::write::Encoder::new(File::create(&zstd_path).unwrap(), 1).unwrap();
        encoder.write_all(b"five\nsix\n").unwrap();
        encoder.finish().unwrap();
        let mut lines = read_text_lines(&zstd_path).unwrap();
        assert_eq!(lines.next().unwrap().unwrap(), "five");
        assert_eq!(lines.next().unwrap().unwrap(), "six");
    }
}
