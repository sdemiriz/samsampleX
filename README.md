# samsampleX
A Python-based tool for customizable BAM file downsampling. Sample reads from a source BAM file to match the depth of coverage distribution of one or more template BAM file(s).

## Features:
- Reproducable, deterministic downsampling using integer seeds.
- Aggregate depth selection using multiple BED templates via select metrics (`min`, `max`, `mean`, `random`).
- BED template compression and/or smoothing.
- Calculation of quality metrics:
    - Wasserstein distance: distribution-wide downsampling performance.
    - Total Variation: per-base downsampling performance.
- Depth comparison plotting for visual sampling comparisons, with the option to emit a TSV file of the same data instead.

## Installation
### Requirements

- pysam
- xxHash
- numpy
- Snakemake (benchmarking only)

### Build samsampleX
```bash
git clone https://github.com/sdemiriz/samsampleX.git
cd samsampleX
pip install .
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
| `test_metrics.c` | Total Variation, Wasserstein distance calculations |
| `integration_tests.sh` | CLI help, subcommand arguments, error handling |

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
6. Report metrics (Total Variation, Wasserstein distance 1)

## Metrics
| Metric | Significance |
| ------ | ------------ |
| Wasserstein-1 Distance | Difference between depth distributions |
| Total Variation | Per-position depth difference |


## Benchmarking
Benchmarking is done by a `snakemake` workflow in the `benchmarks` directory, and thus `snakemake` should be installed beforehand (for HPC systems, also install `snakemake-executor-plugin-slurm` or other plugin compatible with your system type). 

An `Apptainer` container definition `bench.def` that contains installs for `GATK`, `samtools`, `sambamba` and `samsampleX` is included. Build this container using `apptainer build bench.sif bench.def` before running the workflow.

Configure the benchmarking parameters in `config.yaml` in the same directory: copy and rename an existing chunk with all parameters and populate the values. All input files are expected to be found in the same directory as `config.yaml`, BAM files should be indexed using `samtools index`.

```{yaml}
# config.yaml
benchmarks:             # all chunks should be children of this header
  wgs-chr21:            # arbitrary name for benchmarking instance, parameters will be children
    chr: "chr21"        # specify contig
    start: 1            # region start coordinate
    end: 46709982       # region end coordinate
    seed: 42            # random seed (base)
    n_replicates: 1     # replicate count, will affect seed 
                        # (e.g. seed=42, n=3 will use seeds 43, 44, 45)
    collapse: 0         # define smoothing during mapping step
    templates:          # specify files to use as templates in sampling
      - "template.bam"  # all files must be in the benchmarks directory

    mode: "mean"        # how to determine per-position template depths from multiple template files
    source: "source.bam" # specify file to downsample

    coefficient: 0.1    # coefficient provided to GATK, samtools, sambamba

    cpu: 1              # specify hardware resource (used by all steps)
    mem_mb: 16384
    time: "10:00"
```

When executing the workflow, navigate to the `benchmarks` directory and make sure to use the following arguments:
```
snakemake -p --use-apptainer --apptainer-args '--bind $(pwd)'
```

A directory for all intermediate files will be created for each chunk defined in `config.yaml` and the final benchmark results will be made available in the `benchmarks` directory as `benchmark-{chunk_name}.tsv`.
