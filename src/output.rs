use std::io::{self, Write};
use crate::types::TandemRepeat;

const HEADER: &str = "#seq_name\tstart\tend\tperiod\tcopies\tidentity\tindel_rate\tscore\tA_freq\tC_freq\tG_freq\tT_freq\tentropy";

/// Write the TAB-separated header line.
pub fn write_header<W: Write>(writer: &mut W) -> io::Result<()> {
    writeln!(writer, "{}", HEADER)
}

/// Write a single TandemRepeat record in TAB format.
pub fn write_record<W: Write>(writer: &mut W, tr: &TandemRepeat) -> io::Result<()> {
    writeln!(
        writer,
        "{}\t{}\t{}\t{:.0}\t{:.2}\t{:.4}\t{:.4}\t{:.1}\t{:.4}\t{:.4}\t{:.4}\t{:.4}\t{:.4}",
        tr.seq_name,
        tr.start,
        tr.end,
        tr.period,
        tr.copies,
        tr.identity,
        tr.indel_rate,
        tr.score,
        tr.composition[0],
        tr.composition[1],
        tr.composition[2],
        tr.composition[3],
        tr.entropy,
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_header() {
        let mut buf = Vec::new();
        write_header(&mut buf).unwrap();
        let output = String::from_utf8(buf).unwrap();
        assert!(output.starts_with('#'));
        assert!(output.contains("seq_name"));
        assert!(output.contains("identity"));
    }

    #[test]
    fn test_write_record() {
        let tr = TandemRepeat {
            seq_name: "chr1".to_string(),
            start: 100,
            end: 200,
            period: 4,
            copies: 25.0,
            consensus: b"ACGT".to_vec(),
            identity: 0.95,
            indel_rate: 0.02,
            score: 85.0,
            composition: [0.25, 0.25, 0.25, 0.25],
            entropy: 2.0,
        };
        let mut buf = Vec::new();
        write_record(&mut buf, &tr).unwrap();
        let output = String::from_utf8(buf).unwrap();
        assert!(output.contains("chr1"));
        assert!(output.contains("100"));
        assert!(output.contains("200"));
        assert!(output.contains("0.95"));
    }

    #[test]
    fn test_output_format() {
        let tr = TandemRepeat {
            seq_name: "test".to_string(),
            start: 42,
            end: 142,
            period: 5,
            copies: 20.0,
            consensus: b"ACGTA".to_vec(),
            identity: 0.85,
            indel_rate: 0.05,
            score: 82.5,
            composition: [0.2, 0.3, 0.3, 0.2],
            entropy: 1.97,
        };
        let mut buf = Vec::new();
        write_record(&mut buf, &tr).unwrap();
        let output = String::from_utf8(buf).unwrap();
        // Verify TAB-separated with correct field count (13 fields)
        let fields: Vec<&str> = output.trim().split('\t').collect();
        assert_eq!(fields.len(), 13, "Should have 13 tab-separated fields");
    }
}
