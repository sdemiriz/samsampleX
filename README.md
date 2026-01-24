# samsampleX

A C tool for depth-aware BAM file sampling. Sample reads from a source BAM file to match the depth distribution of a template BAM file.

## Features

- **map**: Extract depth of coverage from a BAM file to a BED template
- **sample**: Sample reads from a BAM to match template depth distribution
- Deterministic sampling using xxHash for reproducible results
- Support for multiple template BED files with combine modes (min, max, mean, random)
- Calculates quality metrics (Wasserstein distance, MAE)
- Collapse option for compact BED output

## Requirements

- **htslib** (for BAM I/O)
- **samtools** (for sorting/indexing, must be in PATH)
- C compiler (gcc or clang)
- make

## Installation

### 1. Install htslib

```bash
# Ubuntu/Debian
sudo apt install libhts-dev

# macOS (Homebrew)
brew install htslib

# From source
git clone https://github.com/samtools/htslib
cd htslib
make
sudo make install
```

### 2. Build samsampleX

```bash
git clone <repository>
cd samsampleX-c
make
```

### 3. Install (optional)

```bash
# System-wide (requires root)
sudo make install

# User-local install
make PREFIX=$HOME/.local install
```

## Testing

Run the full test suite:

```bash
make test
```

Or run specific test types:

```bash
# Unit tests only (region parsing, depth arrays, metrics)
make test-unit

# Integration tests only (CLI interface)
make test-integration
```

### Test Coverage

| Test File | Description |
|-----------|-------------|
| `test_parsing.c` | Region string parsing, combine mode parsing |
| `test_depth.c` | Depth array allocation, combine operations |
| `test_metrics.c` | MAE, Wasserstein distance calculations |
| `integration_tests.sh` | CLI help, subcommand arguments, error handling |

## Usage

### Map: Extract depth template

```bash
# Extract depth from template BAM to BED file
samsampleX map --template-bam template.bam --region chr1:1000-2000 --out-bed template.bed

# With collapse (merge consecutive positions with depth difference <= 5)
samsampleX map --template-bam template.bam --region chr1 --collapse 5 --out-bed template.bed
```

### Sample: Depth-aware read sampling

```bash
# Basic usage
samsampleX sample \
    --source-bam high_depth.bam \
    --template-bed template.bed \
    --region chr1:1000-2000 \
    --out-bam sampled.bam

# With multiple templates (combines using 'min' by default)
samsampleX sample \
    --source-bam source.bam \
    --template-bed template1.bed \
    --template-bed template2.bed \
    --region chr1 \
    --mode mean \
    --seed 123 \
    --out-bam sampled.bam

# Skip sorting/indexing
samsampleX sample \
    --source-bam source.bam \
    --template-bed template.bed \
    --region chr1 \
    --no-sort \
    --out-bam sampled.bam
```

## Options

### map subcommand

| Option | Description | Default |
|--------|-------------|---------|
| `--template-bam FILE` | Input BAM file (required) | - |
| `--region REGION` | Target region, samtools-style (required) | - |
| `--out-bed FILE` | Output BED file | `out.bed` |
| `--collapse INT` | Merge consecutive positions with depth diff <= INT | `0` (per-base) |

### sample subcommand

| Option | Description | Default |
|--------|-------------|---------|
| `--source-bam FILE` | Input BAM to sample from (required) | - |
| `--template-bed FILE` | Template BED file(s) (required, repeatable) | - |
| `--region REGION` | Target region, samtools-style (required) | - |
| `--out-bam FILE` | Output BAM file | `out.bam` |
| `--mode MODE` | Combine mode: min, max, mean, random | `min` |
| `--seed INT` | Random seed for reproducibility | `42` |
| `--no-sort` | Skip sorting and indexing output | false |

## Algorithm

### Mapping

1. Parse target region from BAM header
2. Compute per-position depth using read alignments
3. Write to BED4 format (chrom, start, end, depth)
4. Optionally collapse consecutive similar depths

### Sampling

1. Load template depth array from BED file(s)
2. Compute source depth array from BAM
3. Calculate ratio: `ratio[i] = min(1.0, template[i] / source[i])`
4. For each read:
   - Hash read name → fraction in [0, 1)
   - Compute mean ratio over read's positions
   - Keep read if `hash_fraction < mean_ratio`
5. Sort and index output BAM
6. Report metrics (Wasserstein distance, MAE)

## Metrics

After sampling, the tool reports:

- **Mean Template Depth**: Average depth in the template
- **Mean Output Depth**: Average depth in the output
- **MAE**: Mean Absolute Error per position
- **Wasserstein-1 Distance**: Earth Mover's Distance between depth distributions

## License

MIT License

## Acknowledgments

- Uses [xxHash](https://github.com/Cyan4973/xxHash) for fast hashing (BSD-2 license)
- Built on [htslib](https://github.com/samtools/htslib) for BAM I/O

