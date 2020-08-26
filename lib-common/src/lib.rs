/// lib-common -- shared functionality
use rust_htslib::{bam, bcf};

/// Return `bam::Format` for the given filename.
pub fn guess_bam_format(filename: &str) -> bam::Format {
    if filename.ends_with(".bam") {
        bam::Format::BAM
    } else {
        bam::Format::SAM
    }
}

#[derive(Debug)]
pub struct BcfFormatInfo {
    /// Guessed format.
    pub format: bcf::Format,
    /// Guessed compression status.
    pub uncompressed: bool,
}

/// Return `bcf::Format` for the given filename.
pub fn guess_bcf_format(filename: &str) -> BcfFormatInfo {
    if filename.ends_with(".bcf") {
        return BcfFormatInfo {
            format: bcf::Format::BCF,
            uncompressed: false,
        };
    } else if filename.ends_with(".vcf.gz") {
        return BcfFormatInfo {
            format: bcf::Format::VCF,
            uncompressed: false,
        };
    } else {
        return BcfFormatInfo {
            format: bcf::Format::VCF,
            uncompressed: true,
        };
    }
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

    #[test]
    fn test_guess_bcf_format() {
        assert_matches!(guess_bcf_format("ex.vcf").format, bcf::Format::VCF);
        assert_eq!(guess_bcf_format("ex.vcf").uncompressed, true);
        assert_matches!(guess_bcf_format("ex.vcf.gz").format, bcf::Format::VCF);
        assert_eq!(guess_bcf_format("ex.vcf.gz").uncompressed, false);
        assert_matches!(guess_bcf_format("ex.bcf").format, bcf::Format::BCF);
        assert_eq!(guess_bcf_format("ex.bcf").uncompressed, false);
        assert_matches!(guess_bcf_format("ex.xxx").format, bcf::Format::VCF);
        assert_eq!(guess_bcf_format("ex.xxx").uncompressed, true);
    }
}
