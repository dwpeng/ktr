use std::process::Command;

fn run_ktr(args: &[&str]) -> std::process::Output {
    Command::new(env!("CARGO_BIN_EXE_ktr"))
        .args(args)
        .output()
        .expect("Failed to run ktr")
}

fn count_output_lines(output: &str) -> usize {
    output.lines().filter(|l| !l.starts_with('#')).count()
}

#[test]
fn test_ktr_help() {
    let output = Command::new(env!("CARGO_BIN_EXE_ktr"))
        .arg("--help")
        .output()
        .expect("Failed to run ktr --help");

    assert!(output.status.success());
    let stdout = String::from_utf8_lossy(&output.stdout);
    assert!(stdout.contains("Streaming tandem repeat finder"));
    assert!(stdout.contains("--k"));
    assert!(stdout.contains("--period"));
    assert!(stdout.contains("--min-identity"));
}

#[test]
fn test_ktr_on_simple_repeat() {
    let output = run_ktr(&[
        "test_data/simple_repeat.fasta",
        "--k",
        "5",
        "--period",
        "50",
        "--min-run-length",
        "10",
        "--min-matches",
        "2",
        "--min-identity",
        "0.5",
    ]);

    assert!(
        output.status.success(),
        "ktr exited with error: {}",
        String::from_utf8_lossy(&output.stderr)
    );

    let stdout = String::from_utf8_lossy(&output.stdout);

    // Should find the ACGT tandem repeat (period ~4)
    assert!(
        stdout.contains("simple_repeat"),
        "Should output the sequence name"
    );

    // Verify the output contains a valid tab-separated record with period 4.
    // Expected format: #seq_name\tstart\tend\tperiod\tcopies\t...
    // Data line: simple_repeat\t<start>\t<end>\t4\t<copies>\t...
    let has_period_4 = stdout.lines().any(|line| {
        if line.starts_with('#') {
            return false;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        fields.len() == 13 && fields[0] == "simple_repeat" && fields[3] == "4"
    });
    assert!(
        has_period_4,
        "Expected a tab-separated output line with period=4 for simple_repeat.\nstdout:\n{}",
        stdout
    );
}

#[test]
fn test_ktr_multi_k() {
    let output = run_ktr(&[
        "test_data/simple_repeat.fasta",
        "--k",
        "3,5",
        "--period",
        "50",
        "--min-run-length",
        "10",
        "--min-matches",
        "2",
        "--min-identity",
        "0.5",
    ]);
    assert!(
        output.status.success(),
        "ktr multi-k failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    let stdout = String::from_utf8_lossy(&output.stdout);
    let has_period_4 = stdout.lines().any(|line| {
        if line.starts_with('#') {
            return false;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        fields.len() == 13 && fields[0] == "simple_repeat" && fields[3] == "4"
    });
    assert!(
        has_period_4,
        "Multi-k should detect period-4 repeat.\nstdout:\n{}",
        stdout
    );
}

#[test]
fn test_ktr_gap_tolerance() {
    // 80bp sequence: 19 perfect ACGT copies + 1 mutated copy (ACAA at position 12-15)
    let mut path = std::env::temp_dir();
    path.push("ktr_test_mutated.fasta");
    let seq = b"ACGTACGTACGTACAACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
    let mut fasta = b">mutated\n".to_vec();
    fasta.extend_from_slice(seq);
    fasta.push(b'\n');
    std::fs::write(&path, &fasta).unwrap();

    let output = run_ktr(&[
        path.to_str().unwrap(),
        "--k",
        "5",
        "--min-run-length",
        "40",
        "--min-score",
        "0",
    ]);
    assert!(
        output.status.success(),
        "ktr failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    let stdout = String::from_utf8_lossy(&output.stdout);
    let n_results = count_output_lines(&stdout);
    assert!(
        n_results >= 1,
        "Gap tolerance should detect mutated repeat, found {} results.\nstdout:\n{}",
        n_results,
        stdout
    );
    std::fs::remove_file(&path).ok();
}
