# samsampleX
A Python-based tool for customizable BAM file downsampling. Sample reads from a source BAM file to match the depth of coverage distribution of one or more template BAM file(s).

## Features:
- Reproducable, deterministic downsampling using integer seeds.
- Uniform sampling mode: retain a fixed fraction of reads by hash, bypassing template/depth logic.
- Map depth from multiple BAM files to a single BED template using combine modes (`min`, `max`, `mean`, `median`, `random`).
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
Extract depth of coverage from one or more template BAM file(s) to a single BED template. When multiple BAMs are provided, per-position depths are combined using the selected `--mode`.
```bash
# Single BAM
samsampleX map \
    --template-bam template.bam \
    --region chr1:1000-2000 \
    --out-bed template.bed

# Multiple BAMs (combined per-position using mean)
samsampleX map \
    --template-bam a.bam b.bam c.bam \
    --region chr1:1000-2000 \
    --mode mean \
    --out-bed template.bed
```
| Option | Description | Default |
|--------|-------------|---------|
| `--template-bam FILE [FILE ...]` | Input BAM file(s) (required) | - |
| `--region REGION` | Target region, samtools-style (required) | - |
| `--out-bed FILE` | Output BED file | `out.bed` |
| `--collapse INT` | Merge consecutive positions with depth diff <= INT | `0` (per-position) |
| `--mode MODE` | Combine mode when multiple BAMs: `min`, `max`, `mean`, `median`, `random` | `mean` |
| `--seed INT` | Random seed for `--mode random` | `42` |

### Sampling
Downsample BAM based on provided BED template(s), using selected metric if multiple BEDs provided. Alternatively, use `--uniform` for position-independent uniform sampling by read-name hash.

**Depth-based sampling (template required):**
```bash
samsampleX sample \
    --source-bam high_depth.bam \
    --template-bed template.bed \
    --region chr1:1000-2000 \
    --out-bam sampled.bam
```

**Uniform sampling (no template):**
```bash
samsampleX sample \
    --source-bam high_depth.bam \
    --uniform 0.5 \
    --region chr1:1000-2000 \
    --out-bam sampled.bam
```
Retains approximately 50% of reads uniformly across the region.

| Option | Description | Default |
|--------|-------------|---------|
| `--source-bam FILE` | Input BAM to sample from (required) | - |
| `--template-bed FILE` | Template BED file(s); required unless `--uniform` is used | - |
| `--uniform FRACTION` | Uniform sampling: retain fraction of reads by hash (0–1). Bypasses template/depth logic. | - |
| `--region REGION` | Target region, samtools-style (required) | - |
| `--out-bam FILE` | Output BAM file | `out.bam` |
| `--mode MODE` | Combine mode for multiple templates: `min`, `max`, `mean`, `random` | `random` |
| `--stat STAT` | Statistic for summarising ratio over read span: `mean`, `min`, `max`, `median` | `mean` |
| `--seed INT` | Random seed for reproducibility | `42` |
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

### Mapback
Remap HLA\*LA PRG-mapped reads back to canonical chr6 coordinates. This is a preprocessing step for BAM files produced by HLA\*LA, which maps reads to a pangenome reference graph (PRG) with synthetic contig names (`PRG_1`, `PRG_2`, ...). The mapback subcommand translates these back to chr6 positions using the HLA\*LA `sequences.txt` file and known HLA gene / alt contig boundaries.

The output BAM can then be used as input to `sample` for depth-aware downsampling on chr6.

```bash
# Step 1: remap PRG reads to chr6
samsampleX mapback \
    --source-bam hlala_output.bam \
    --region chr6:28000000-34000000 \
    --genome-build GRCh38 \
    --out-bam remapped.bam

# Step 2: sample from the remapped BAM
samsampleX sample \
    --source-bam remapped.bam \
    --template-bed template.bed \
    --region chr6:28000000-34000000 \
    --out-bam sampled.bam
```

| Option | Description | Default |
|--------|-------------|---------|
| `--source-bam FILE` | HLA\*LA-remapped BAM file (required) | - |
| `--region REGION` | Target region on chr6, samtools-style (required) | - |
| `--out-bam FILE` | Output BAM file | `out.mapback.bam` |
| `--genome-build BUILD` | Reference genome build: `GRCh38` or `GRCh37` (required) | - |
| `--prg-seq FILE` | Path to HLA\*LA `sequences.txt` | `HLA-LA/graphs/PRG_MHC_GRCh38_withIMGT/sequences.txt` |

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

## Example

The following commands 

```bash
cd examples/

# Download three first three 1K Genomes 30X WGS samples from
# https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/1000G_2504_high_coverage.sequence.index
wget ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR323/ERR3239480/NA12718.final.cram -O NA12718.cram && samtools index NA12718.cram
wget ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR323/ERR3239481/NA12748.final.cram -O NA12748.cram && samtools index NA12748.cram
wget ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR323/ERR3239482/NA12775.final.cram -O NA12775.cram && samtools index NA12775.cram

# Download reference genome (GRCh38)
wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa

# Convert to BAM and index
samtools view NA12718.final.cram chr21:10000000-10010000 -b -o NA12718.bam -T GRCh38_full_analysis_set_plus_decoy_hla.fa && samtools index NA12718.bam
samtools view NA12748.final.cram chr21:10000000-10010000 -b -o NA12748.bam -T GRCh38_full_analysis_set_plus_decoy_hla.fa && samtools index NA12748.bam
samtools view NA12775.final.cram chr21:10000000-10010000 -b -o NA12775.bam -T GRCh38_full_analysis_set_plus_decoy_hla.fa && samtools index NA12775.bam

# Run samsampleX workflow
samsampleX map \
    --template-bam NA12718.bam NA12748.bam NA12775.bam \
    --region chr21:10000000-10010000 \
    --mode mean \
    --collapse 0 \
    --out-bed template.bed
# template.bed should match example-template.bed

# Source BAM+index is provided in the examples directory, created by subsetting to target region from
# https://ftp.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/NIST_HiSeq_HG002_Homogeneity-10953946/NHGRI_Illumina300X_AJtrio_novoalign_bams/HG002.GRCh38.300x.bam
samsampleX sample \
    --source-bam HG002.GRCh38.300x.chr21:10000000-10010000.bam \
    --template-bed template.bed \
    --region chr21:10000000-10010000 \
    --seed 42 \
    --out-bam sampled.bam

samtools index sampled.bam

samsampleX plot \
    --source-bam HG002.GRCh38.300x.chr21:10000000-10010000.bam \
    --template-bed template.bed \
    --out-bam sampled.bam \
    --region chr21:10000000-10010000 \
    --out-png plot.png
# plot.png should match example-plot.png
```

## Testing

A `pytest` test suite is available. Run using the `-v` flag for a detailed report.
```bash
pytest -v
```

## Algorithm rundown

### Mapping
1. Parse target region from first BAM header
2. Compute per-position depth of coverage for each BAM over the region
3. If multiple BAMs: combine depths per-position using `--mode` (min, max, mean, median, or random)
4. Write to BED4 format (`chrom`, `start`, `end`, `depth` columns)
5. Optionally collapse consecutive similar depths (`--collapse`)

### Sampling
1. **Uniform mode** (`--uniform FRACTION`): Skip template/depth logic. For each read, hash the read name with xxHash32 to get $f_{read} \in [0, 1)$; keep if $f_{read} < FRACTION$. Deterministic and position-independent.
2. **Depth-based mode**: Load template depths from BED file(s); if multiple templates are provided, combine them per-position using the selected `--mode`
3. Compute source depths from BAM
4. Calculate per-position sampling ratio: $ratio(i) = \min(1,\; depth_{template}(i) \;/\; depth_{source}(i))$
   - Positions where the template depth meets or exceeds the source depth get ratio 1.0 (keep all reads)
   - Positions with zero source depth get ratio 0.0
5. Build a cumulative sum of the ratio array for O(1) range queries
6. For each read in the source BAM:
   - Hash read name with xxHash32 to produce a deterministic fraction $f_{read} \in [0, 1)$
   - Summarise the ratio over the read's covered positions using `--stat` (default: mean via cumsum lookup)
   - Keep the read if $f_{read} < ratio_{read}$
7. Report metrics (depth-based mode only): Total Variation and Wasserstein-1 distance (unless `--no-metrics`)

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
