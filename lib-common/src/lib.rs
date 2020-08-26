/// lib-common -- shared functionality
use rust_htslib::bam;

/// Return `bam::Format` for the given filename.
pub fn guess_bam_format(filename: &str) -> bam::Format {
    if filename.ends_with(".bam") {
        bam::Format::BAM
    } else {
        bam::Format::SAM
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use matches::assert_matches;

    #[test]
    fn guesss_bam_format() {
        assert_matches!(guess_bam_format("ex.sam"), bam::Format::SAM);
        assert_matches!(guess_bam_format("ex.bam"), bam::Format::BAM);
        assert_matches!(guess_bam_format("ex.xxx"), bam::Format::SAM);
    }
}
