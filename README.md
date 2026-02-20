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
- matplotlib
- Snakemake (benchmarking only)
- pytest (testing only)

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
| `--mode MODE` | Combine mode for multiple templates: `min`, `max`, `mean`, `random` | `random` |
| `--stat STAT` | Statistic for summarising ratio over read span: `mean`, `min`, `max`, `median` | `mean` |
| `--seed INT` | Random seed for reproducibility | `42` |
| `--no-sort` | Skip sorting and indexing output | false |
| `--no-metrics` | Skip metrics calculation after sampling | false |

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

### Stats
Compare depth distributions between two BAM files over a given region. Reports mean depth for each BAM, Total Variation distance, and normalised Wasserstein-1 distance.
```bash
samsampleX stats \
    --bam-a template.bam \
    --bam-b sampled.bam \
    --region chr1:1000-2000
```

| Option | Description | Default |
|--------|-------------|---------|
| `--bam-a FILE` | First BAM file, e.g. reference/template (required) | - |
| `--bam-b FILE` | Second BAM file, e.g. sampled output (required) | - |
| `--region REGION` | Target region, samtools-style (required) | - |

## Testing

A `pytest` test suite is available. Run using the `-v` flag for a detailed report.
```bash
pytest -v
```

## Algorithm rundown

### Mapping
1. Parse target region from source BAM header
2. Compute per-position depth of coverage for region
3. Write to BED4 format (`chrom`, `start`, `end`, `depth` columns)
4. Optionally collapse consecutive similar depths (`--collapse`)

### Sampling
1. Load template depths from BED file(s); if multiple templates are provided, combine them per-position using the selected `--mode`
2. Compute source depths from BAM
3. Calculate per-position sampling ratio: $ratio(i) = \min(1,\; depth_{template}(i) \;/\; depth_{source}(i))$
   - Positions where the template depth meets or exceeds the source depth get ratio 1.0 (keep all reads)
   - Positions with zero source depth get ratio 0.0
4. Build a cumulative sum of the ratio array for O(1) range queries
5. For each read in the source BAM:
   - Hash read name with xxHash32 to produce a deterministic fraction $f_{read} \in [0, 1)$
   - Summarise the ratio over the read's covered positions using `--stat` (default: mean via cumsum lookup)
   - Keep the read if $f_{read} < ratio_{read}$
6. Sort and index output BAM (unless `--no-sort`)
7. Report metrics: Total Variation and Wasserstein-1 distance (unless `--no-metrics`)

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
