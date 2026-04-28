# ktr — Streaming Tandem Repeat Finder

ktr finds tandem repeats (TRs) in genomic sequences using a two-phase approach:
streaming kmer-based scanning followed by TRF-style wraparound dynamic programming
alignment.

## Quick Start

```bash
# Build (requires Rust nightly for SIMD)
RUSTFLAGS="-C target-cpu=native" cargo build --release

# Run on a FASTA genome
./target/release/ktr genome.fasta

# Output to file
./target/release/ktr genome.fasta -o results.tsv
```

## Output Format

13-column TAB-separated (plus 2 optional debug columns):

| Column | Field      | Description                          |
|--------|------------|--------------------------------------|
| 1      | seq_name   | Sequence name (first word of header) |
| 2      | start      | TR start (1-based)                   |
| 3      | end        | TR end (1-based, inclusive)          |
| 4      | period     | Repeat unit length (bp)              |
| 5      | copies     | Estimated copy number                |
| 6      | identity   | Average copy identity (0–1)          |
| 7      | indel_rate | Average indel rate (0–1)             |
| 8      | score      | Quality score `identity×100−indel×50`|
| 9–12   | A/C/G/T    | Base composition of consensus        |
| 13     | entropy    | Shannon entropy of base composition  |

With `--debug`, two extra columns: consensus sequence and full TR sequence.

## Algorithm

### Phase 1: Streaming Kmer Scan

The sequence is split into 500 Kbp chunks (configurable). Each chunk is scanned
independently with a sliding 5-mer (configurable) window:

1. **Kmer encoding** — Each 5-mer is encoded as a 2-bit integer (A=00, C=01,
   G=10, T=11), packed into a u64.
2. **Period voting** — When a kmer reappears within `max_period` bases, the
   distance is recorded as a vote for that period.
3. **Run tracking** — Consecutive period votes form a run. When a run ends (a
   kmer not seen within `max_period`), candidates are emitted if the vote
   concentration (`votes / effective_run_len`) exceeds `min_concentration`.
4. **Chunk overlap** — Adjacent chunks overlap by `max_period × 2` bases to
   capture boundary-spanning TRs. Overlapping candidates are merged in Phase 2.

### Phase 2: WDP Alignment

Each candidate is validated with wraparound dynamic programming:

1. **Consensus extraction** — The first `period` bases of the candidate region
   form the initial consensus (taken from the candidate start, not the context
   start with flanks — a common bug in naive implementations).
2. **Iterative alignment** — Subsequent period-length windows are aligned to the
   consensus using Smith-Waterman. Successful alignments extend the TR; the
   consensus is refined via merging.
3. **Scoring** — Average identity, indel rate, and TRF-style score are computed
   from the aligned copies.
4. **Merging** — Overlapping or adjacent TRs (gap ≤ period × 20) with similar
   periods (ratio ≤ 2.0) are merged into a single record.

## Parameters

### Detection Sensitivity

| Flag                    | Default | Description                              |
|-------------------------|---------|------------------------------------------|
| `--k`                   | 5       | Kmer length (1–32). Longer kmers reduce  |
|                         |         | false positives but miss short periods.  |
| `--max-period`          | 5000    | Maximum TR period to detect. Also sets   |
|                         |         | sliding window = max_period × 2.         |
| `--min-run-length`      | 40      | Minimum run length (bases) for Phase 1.  |
| `--min-matches`         | 3       | Minimum period votes for a candidate.    |
| `--min-concentration`   | 0.20    | Minimum vote concentration for a run     |
|                         |         | (0–1). Lower = more sensitive, higher =  |
|                         |         | more precise.                            |
| `--min-identity`        | 0.70    | Minimum copy identity for Phase 2 (0–1). |
| `--min-score`           | 0       | Minimum score for output filtering.      |
|                         |         | Score ≈ 55–100 for typical TRs.          |

### Performance

| Flag             | Default  | Description                             |
|------------------|----------|-----------------------------------------|
| `--threads`      | 0        | Thread count (0 = all available). Both  |
|                  |          | Phase 1 chunks and Phase 2 candidates   |
|                  |          | are parallelized via Rayon.             |
| `--chunk-size`   | 500000   | Chunk size for Phase 1 (0 = no chunk).  |

### Output

| Flag        | Default | Description                           |
|-------------|---------|---------------------------------------|
| `-o`        | stdout  | Output file path.                     |
| `--debug`   | off     | Append consensus and TR sequence to   |
|             |         | each output row.                      |

## Comparison with TRF

ktr is inspired by Tandem Repeats Finder (TRF) but uses a different Phase 1
approach. Key differences:

- **Phase 1** — TRF uses alignment-based detection; ktr uses kmer voting. This
  makes ktr faster for large genomes but unable to detect TRs with period <
  kmer length.
- **Parallelism** — ktr parallelizes both Phase 1 (chunk-level) and Phase 2
  (candidate-level) via Rayon.
- **Adaptive alignment** — WDP parameters are automatically scaled by period
  length: log₂(period)-based scaling softens penalties for long-period TRs.
- **Scoring** — ktr uses a simplified score (`identity×100−indel×50`). TRF's
  probabilistic scoring accounts for additional factors.

For a typical bacterial genome (E. coli), ktr detects 50–60 TRs with period ≥ 5,
covering ~80% of TRF-detected TRs (period ≥ 40). Most misses are microsatellites
(period < 20) below the kmer detection threshold.

## Building

```bash
# Native CPU optimizations (recommended)
RUSTFLAGS="-C target-cpu=native" cargo build --release

# Test
cargo test
```

The helicase crate requires AVX2 (x64) or NEON (aarch64). On systems without
these, pass `--features scalar` to helicase.

## License

MIT
