use std::process::Command;

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
    assert!(stdout.contains("--max-period"));
    assert!(stdout.contains("--min-identity"));
}

#[test]
fn test_ktr_on_simple_repeat() {
    let output = Command::new(env!("CARGO_BIN_EXE_ktr"))
        .arg("test_data/simple_repeat.fasta")
        .arg("--k")
        .arg("5")
        .arg("--max-period")
        .arg("50")
        .arg("--min-run-length")
        .arg("10")
        .arg("--min-matches")
        .arg("2")
        .arg("--min-identity")
        .arg("0.5")
        .output()
        .expect("Failed to run ktr");

    assert!(output.status.success(), "ktr exited with error: {}",
        String::from_utf8_lossy(&output.stderr));

    let stdout = String::from_utf8_lossy(&output.stdout);
    // Should find the ACGT tandem repeat (period ~4)
    assert!(stdout.contains("simple_repeat"), "Should output the sequence name");
}
