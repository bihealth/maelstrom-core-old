pub mod library;

use std::path::Path;

use bio_types::genome::Interval;
use regex::Regex;
use rust_htslib::{bam, bam::Read};

use super::error::Error;

/// Return `bam::Format` for the given filename.
pub fn guess_bam_format(filename: &str) -> bam::Format {
    if filename.ends_with(".bam") {
        bam::Format::BAM
    } else {
        bam::Format::SAM
    }
}

/// Generate list of all contigs from BAM header.
pub fn build_chroms_bam(
    header: &bam::HeaderView,
    regex: Option<String>,
) -> Result<Vec<Interval>, Error> {
    let mut result: Vec<Interval> = Vec::new();

    let re = if let Some(regex) = regex {
        Some(Regex::new(&regex).unwrap())
    } else {
        None
    };

    for (no, name) in header.target_names().iter().enumerate() {
        let name = std::str::from_utf8(name)?.to_string();
        let len = header.target_len(no as u32).unwrap() as u64;
        let do_push = match &re {
            Some(re) => re.is_match(&name),
            None => true,
        };
        if do_push {
            result.push(Interval::new(name, 0..len));
        }
    }

    Ok(result)
}

/// Parse @RG lane into triple (id, sm).
fn parse_line_rg(line: String) -> Option<(String, String)> {
    let line_split = line.split('\t');
    let mut id: Option<String> = None;
    let mut sm: Option<String> = None;
    for s in line_split {
        let token: Vec<&str> = s.split(':').collect();
        if token.len() >= 2 {
            match token[0] {
                "ID" => {
                    id = Some(token[1].to_string());
                }
                "SM" => {
                    sm = Some(token[1].to_string());
                }
                _ => (),
            }
        }
    }

    match (id, sm) {
        (Some(id), Some(sm)) => Some((id, sm)),
        _ => None,
    }
}

/// Construct sample list from BAM header.
fn samples_from_header(header: &[u8]) -> Vec<String> {
    let text = String::from_utf8(header.to_vec()).unwrap();
    let mut samples = Vec::new();

    for line in text.lines() {
        if line.starts_with("@RG") {
            if let Some((_id, sm)) = parse_line_rg(line.to_string()) {
                if !samples.contains(&sm) {
                    samples.push(sm.clone());
                }
            }
        }
    }

    samples
}

/// Construct sample list from BAM file.
pub fn samples_from_file<P: AsRef<Path>>(path: P) -> Result<Vec<String>, Error> {
    let path = path.as_ref().to_str().unwrap();
    let reader = bam::Reader::from_path(path)?;

    Ok(samples_from_header(&Vec::from(reader.header().as_bytes())))
}

#[cfg(test)]
mod tests {
    use super::*;
    use matches::assert_matches;

    #[test]
    fn test_guess_bam_format() {
        assert_matches!(guess_bam_format("ex.sam"), bam::Format::SAM);
        assert_matches!(guess_bam_format("ex.bam"), bam::Format::BAM);
        assert_matches!(guess_bam_format("ex.xxx"), bam::Format::SAM);
    }
}
