# samsampleX
A C-based tool for customizable BAM file downsampling. Sample reads from a source BAM file to match the depth of coverage distribution of one or more template BAM file(s).

## Features:
- Reproducable, deterministic downsampling using integer seeds.
- Aggregate depth selection using multiple BED templates via select metrics (`min`, `max`, `mean`, `random`).
- BED template compression and/or smoothing.
- Calculation of quality metrics:
    - Wasserstein distance: distribution-wide downsampling performance.
    - MAE: per-base downsampling performance.
- Depth comparison plotting for visual sampling comparisons, with the option to emit a TSV file of the same data instead.

## Installation
### Requirements

- htslib
- samtools
- C compiler (gcc/clang)
- make

### Build samsampleX
```bash
git clone https://github.com/sdemiriz/samsampleX.git
cd samsampleX
make
```

### Install
```bash
# System-wide
sudo make install
```

## Usage
### Mapping
Extract depth of coverage from template BAM to BED template with optional smoothing.
```bash
samsampleX map \
    --template-bam template.bam \
    --region chr1:1000-2000 \
    --out-bed template.bed
```
| Option | Description | Default |
|--------|-------------|---------|
| `--template-bam FILE` | Input BAM file (required) | - |
| `--region REGION` | Target region, samtools-style (required) | - |
| `--out-bed FILE` | Output BED file | `out.bed` |
| `--collapse INT` | Merge consecutive positions with depth diff <= INT | `0` (per-position) |

### Sampling
Downsample BAM based on provided BED template(s), using selected metric if multiple BEDs provided.
```bash
samsampleX sample \
    --source-bam high_depth.bam \
    --template-bed template.bed \
    --region chr1:1000-2000 \
    --out-bam sampled.bam
```

| Option | Description | Default |
|--------|-------------|---------|
| `--source-bam FILE` | Input BAM to sample from (required) | - |
| `--template-bed FILE` | Template BED file(s) (>=1 required) | - |
| `--region REGION` | Target region, samtools-style (required) | - |
| `--out-bam FILE` | Output BAM file | `out.bam` |
| `--mode MODE` | Combine mode: `min`, `max`, `mean`, `random` | `random` |
| `--seed INT` | Random seed for reproducibility | `42` |
| `--no-sort` | Skip sorting and indexing output | false |

### Plotting
Compare depth of coverage between source, template, and output BAM files. Output either as PNG plot or TSV data.

Blue is source, green is template and red is output depth.

TSV contains a column for `position`, and three for respective depths of source, template and output.
```bash
# Generate PNG plot
samsampleX plot \
    --source-bam high_depth.bam \
    --template-bam template.bam \
    --out-bam sampled.bam \
    --region chr1:1000-2000 \
    --out-png coverage_plot.png
```

| Option | Description | Default |
|--------|-------------|---------|
| `--source-bam FILE` | Source BAM file (required) | - |
| `--template-bam FILE` | Template BAM file (mutually exclusive with --template-bed) | - |
| `--template-bed FILE` | Template BED file (mutually exclusive with --template-bam) | - |
| `--out-bam FILE` | Output BAM file from sampling (required) | - |
| `--region REGION` | Target region, samtools-style (required) | - |
| `--out-png FILE` | Output PNG plot (mutually exclusive with --out-tsv) | - |
| `--out-tsv FILE` | Output TSV data (mutually exclusive with --out-png) | - |

## Testing
```bash
# All tests
make test

# Unit tests
make test-unit

# Integration tests
make test-integration
```

### Coverage
| Test File | Description |
|-----------|-------------|
| `test_parsing.c` | Region string parsing, combine mode parsing |
| `test_depth.c` | Depth array allocation, combine operations |
| `test_metrics.c` | MAE, Wasserstein distance calculations |
| `integration_tests.sh` | CLI help, subcommand arguments, error handling |

## Benchmarking
Benchmarking is carried out by a `snakemake` pipeline in the `benchmarks` directory, and is expected to be installed. The GNU `time` utility is used for gathering the runtime metrics for benchmarking. (Note: This is different from the typical `time` command.)

`samsampleX` is compared in performance against `GATK`, `samtools` and `sambamba`, all of which should be installed and available for benchmarking. Be advised that these three tools downsample uniformly against `samsampleX`'s depth-aware method, thus the comparison in performance is not entirely fair due to the difference in functionality.

Configure the benchmarking details at the top of the `Snakefile` to set the genomic region and the count for each tool to be run. Also specify the hardware requirements in terms of cores and memory allocated to each tool. Finally, execute `snakemake` workflow.

The output files `results.tsv` and `results-summary.tsv` should contain the details of the run. For clarity, both `samsampleX` steps `map` and `sample` are benchmarked separately.

`GATK`, `samtools` and `sambamba` 
## Algorithm rundown

### Mapping
1. Parse target region from source BAM header
2. Compute per-position depth of coverage for region
3. Write to BED4 format (`chrom`, `start`, `end`, `depth` columns)
4. Optionally collapse consecutive similar depths (`--collapse`)

### Sampling
1. Load template depths from BED file(s)
2. Compute source depths from BAM
3. Calculate $ratio_{depth}= depth_{template} / depth_{source}$
    - Sample all reads if $ratio_{depth} > 1.0$
4. For each read:
   - Hash read name to a fraction $f_{read}$ in [0, 1)
   - Compute read's mean ratio based on covered positions $\mu(ratio_{depth}) = \sum_{i}^{i+l_{read}} ratio_{depth}(i) / l_{read}$
   - Select and sample read if $\mu(ratio_{depth}) < f_{read}$
5. Sort and index output BAM
6. Report metrics (Wasserstein distance, MAE)

## Metrics
| Metric | Significance |
| ------ | ------------ |
| Wasserstein-1 Distance | Difference between depth distributions |
| MAE | Mean Absolute Error per position |
| Mean Template Depth | Average depth in the template |
| Mean Output Depth | Average depth in the output |


## Acknowledgments
- Uses [xxHash](https://github.com/Cyan4973/xxHash) for fast hashing (BSD-2 license)
- Uses [pbPlots](https://github.com/InductiveComputerScience/pbPlots) for plotting
- Built on [htslib](https://github.com/samtools/htslib)
