/*
 * sample.c - Sample subcommand for depth-aware BAM sampling
 *
 * This module samples reads from a source BAM file based on depth
 * distribution from one or more template BED files. Reads are selected
 * using a deterministic hash-based approach that allows for reproducible
 * sampling.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <getopt.h>

#include <htslib/sam.h>
#include <htslib/hts.h>

#include "sample.h"
#include "bed.h"
#include "depth.h"
#include "metrics.h"
#include "xxhash.h"
#include "samsampleX.h"

/*
 * Print sample subcommand usage.
 */
void sample_usage(void) {
    fprintf(stderr, "Usage: %s sample [options]\n\n", SAMSAMPLEX_NAME);
    fprintf(stderr, "Sample reads from a BAM file to match depth distribution from template BED(s).\n\n");
    fprintf(stderr, "Required options:\n");
    fprintf(stderr, "  --source-bam FILE     Input BAM file to sample reads from\n");
    fprintf(stderr, "  --template-bed FILE   Template BED file(s) with depth values (can be repeated)\n");
    fprintf(stderr, "  --region REGION       Target region (samtools-style)\n\n");
    fprintf(stderr, "Optional:\n");
    fprintf(stderr, "  --out-bam FILE        Output BAM file [default: %s]\n", DEFAULT_OUT_BAM);
    fprintf(stderr, "  --mode MODE           How to combine multiple templates: min, max, mean, random\n");
    fprintf(stderr, "                        [default: min]\n");
    fprintf(stderr, "  --seed INT            Random seed for deterministic sampling [default: %d]\n", DEFAULT_SEED);
    fprintf(stderr, "  --no-sort             Do not sort and index output BAM file\n");
    fprintf(stderr, "  -h, --help            Show this help message\n");
}

/*
 * Compute cumulative sum array from ratio array.
 * Returns array of size (length + 1) where cumsum[0] = 0.
 */
static double *compute_cumsum(const double *ratios, size_t length) {
    double *cumsum;
    SAFE_MALLOC(cumsum, (length + 1) * sizeof(double));
    
    cumsum[0] = 0.0;
    for (size_t i = 0; i < length; i++) {
        cumsum[i + 1] = cumsum[i] + ratios[i];
    }
    
    return cumsum;
}

/*
 * Get mean ratio for a read spanning [read_start, read_end) using cumulative sum.
 */
static double get_mean_ratio(const double *cumsum, int64_t region_start,
                             int64_t region_end, int64_t read_start, int64_t read_end) {
    /* Clip read coordinates to region */
    int64_t clip_start = MAX(read_start, region_start);
    int64_t clip_end = MIN(read_end, region_end);
    
    if (clip_start >= clip_end) {
        return 0.0;
    }
    
    /* Convert to array indices */
    size_t i1 = (size_t)(clip_start - region_start);
    size_t i2 = (size_t)(clip_end - region_start);
    
    double sum = cumsum[i2] - cumsum[i1];
    return sum / (double)(clip_end - clip_start);
}

/*
 * Run the sample subcommand.
 */
int sample_run(sample_args_t *args) {
    int ret = 0;
    
    /* Validate required arguments */
    if (!args->source_bam) {
        fprintf(stderr, "Error: --source-bam is required\n");
        sample_usage();
        return 1;
    }
    if (args->n_template_beds == 0) {
        fprintf(stderr, "Error: At least one --template-bed is required\n");
        sample_usage();
        return 1;
    }
    if (!args->region) {
        fprintf(stderr, "Error: --region is required\n");
        sample_usage();
        return 1;
    }
    
    /* Set defaults */
    if (!args->out_bam) {
        args->out_bam = DEFAULT_OUT_BAM;
    }
    
    fprintf(stderr, "[sample] Source BAM: %s\n", args->source_bam);
    fprintf(stderr, "[sample] Template BEDs: %d file(s)\n", args->n_template_beds);
    for (int i = 0; i < args->n_template_beds; i++) {
        fprintf(stderr, "[sample]   %d: %s\n", i + 1, args->template_beds[i]);
    }
    fprintf(stderr, "[sample] Region: %s\n", args->region);
    fprintf(stderr, "[sample] Mode: %s\n", combine_mode_to_string(args->mode));
    fprintf(stderr, "[sample] Seed: %u\n", args->seed);
    fprintf(stderr, "[sample] Output BAM: %s\n", args->out_bam);
    
    /* Parse region */
    region_t *region = region_parse(args->region);
    if (!region || !region->contig) {
        fprintf(stderr, "Error: Invalid region format: %s\n", args->region);
        return 1;
    }
    
    /* Open source BAM to get header and determine region bounds */
    htsFile *source_fp = hts_open(args->source_bam, "r");
    if (!source_fp) {
        fprintf(stderr, "Error: Cannot open source BAM: %s\n", args->source_bam);
        region_free(region);
        return 1;
    }
    
    sam_hdr_t *header = sam_hdr_read(source_fp);
    if (!header) {
        fprintf(stderr, "Error: Cannot read BAM header\n");
        hts_close(source_fp);
        region_free(region);
        return 1;
    }
    
    /* Load index */
    hts_idx_t *idx = sam_index_load(source_fp, args->source_bam);
    if (!idx) {
        fprintf(stderr, "Error: Cannot load BAM index. Make sure .bai file exists.\n");
        sam_hdr_destroy(header);
        hts_close(source_fp);
        region_free(region);
        return 1;
    }
    
    /* Get contig ID and length */
    int tid = sam_hdr_name2tid(header, region->contig);
    if (tid < 0) {
        fprintf(stderr, "Error: Contig '%s' not found in BAM\n", region->contig);
        hts_idx_destroy(idx);
        sam_hdr_destroy(header);
        hts_close(source_fp);
        region_free(region);
        return 1;
    }
    
    /* Set region bounds */
    if (region->start < 0) region->start = 0;
    if (region->end < 0) region->end = sam_hdr_tid2len(header, tid);
    
    fprintf(stderr, "[sample] Parsed region: %s:%ld-%ld\n",
            region->contig, (long)region->start, (long)region->end);
    
    /* Load template depth arrays from BED files */
    fprintf(stderr, "[sample] Loading template BED file(s)...\n");
    depth_array_t *template_arrays[MAX_TEMPLATE_BEDS];
    for (int i = 0; i < args->n_template_beds; i++) {
        template_arrays[i] = bed_read_depths(args->template_beds[i], region->contig,
                                             region->start, region->end);
        if (!template_arrays[i]) {
            fprintf(stderr, "Error: Failed to load BED file: %s\n", args->template_beds[i]);
            for (int j = 0; j < i; j++) {
                depth_array_free(template_arrays[j]);
            }
            hts_idx_destroy(idx);
            sam_hdr_destroy(header);
            hts_close(source_fp);
            region_free(region);
            return 1;
        }
    }
    
    /* Combine template depths */
    depth_array_t *template_depth;
    if (args->n_template_beds == 1) {
        template_depth = template_arrays[0];
    } else {
        fprintf(stderr, "[sample] Combining %d templates using '%s' mode...\n",
                args->n_template_beds, combine_mode_to_string(args->mode));
        template_depth = bed_combine_depths(template_arrays, args->n_template_beds,
                                            args->mode, args->seed);
        /* Free individual arrays */
        for (int i = 0; i < args->n_template_beds; i++) {
            depth_array_free(template_arrays[i]);
        }
        if (!template_depth) {
            fprintf(stderr, "Error: Failed to combine template depths\n");
            hts_idx_destroy(idx);
            sam_hdr_destroy(header);
            hts_close(source_fp);
            region_free(region);
            return 1;
        }
    }
    
    /* Compute source depth */
    fprintf(stderr, "[sample] Computing source depth array...\n");
    depth_array_t *source_depth = depth_from_bam(args->source_bam, region->contig,
                                                  region->start, region->end);
    if (!source_depth) {
        fprintf(stderr, "Error: Failed to compute source depth\n");
        depth_array_free(template_depth);
        hts_idx_destroy(idx);
        sam_hdr_destroy(header);
        hts_close(source_fp);
        region_free(region);
        return 1;
    }
    
    /* Compute ratio array: ratio[i] = min(1.0, template[i] / source[i]) */
    fprintf(stderr, "[sample] Computing sampling ratios...\n");
    double *ratios;
    SAFE_MALLOC(ratios, source_depth->length * sizeof(double));
    
    for (size_t i = 0; i < source_depth->length; i++) {
        if (source_depth->depths[i] == 0) {
            ratios[i] = 0.0;  /* No reads to sample */
        } else {
            double r = (double)template_depth->depths[i] / (double)source_depth->depths[i];
            ratios[i] = (r > 1.0) ? 1.0 : r;  /* Cap at 1.0 */
        }
    }
    
    /* Compute cumulative sum for efficient range queries */
    double *cumsum = compute_cumsum(ratios, source_depth->length);
    free(ratios);
    
    /* Open output BAM for writing */
    htsFile *out_fp = hts_open(args->out_bam, "wb");
    if (!out_fp) {
        fprintf(stderr, "Error: Cannot open output BAM: %s\n", args->out_bam);
        free(cumsum);
        depth_array_free(source_depth);
        depth_array_free(template_depth);
        hts_idx_destroy(idx);
        sam_hdr_destroy(header);
        hts_close(source_fp);
        region_free(region);
        return 1;
    }
    
    /* Write header */
    if (sam_hdr_write(out_fp, header) < 0) {
        fprintf(stderr, "Error: Failed to write BAM header\n");
        hts_close(out_fp);
        free(cumsum);
        depth_array_free(source_depth);
        depth_array_free(template_depth);
        hts_idx_destroy(idx);
        sam_hdr_destroy(header);
        hts_close(source_fp);
        region_free(region);
        return 1;
    }
    
    /* Iterate through reads and sample */
    fprintf(stderr, "[sample] Sampling reads...\n");
    
    hts_itr_t *iter = sam_itr_queryi(idx, tid, region->start, region->end);
    if (!iter) {
        fprintf(stderr, "Error: Cannot create iterator for region\n");
        ret = 1;
        goto cleanup;
    }
    
    bam1_t *b = bam_init1();
    uint64_t total_reads = 0;
    uint64_t kept_reads = 0;
    
    while (sam_itr_next(source_fp, iter, b) >= 0) {
        /* Skip unmapped reads */
        if (b->core.flag & BAM_FUNMAP) {
            continue;
        }
        
        total_reads++;
        
        /* Get read name and compute hash fraction */
        const char *qname = bam_get_qname(b);
        double hash_fraction = XXH32_fraction(qname, args->seed);
        
        /* Get mean ratio over read's covered positions */
        int64_t read_start = b->core.pos;
        int64_t read_end = bam_endpos(b);
        double mean_ratio = get_mean_ratio(cumsum, region->start, region->end,
                                           read_start, read_end);
        
        /* Keep read if hash_fraction < mean_ratio */
        if (hash_fraction < mean_ratio) {
            if (sam_write1(out_fp, header, b) < 0) {
                fprintf(stderr, "Error: Failed to write read\n");
                ret = 1;
                break;
            }
            kept_reads++;
        }
        
        /* Progress report every 1M reads */
        if (total_reads % 1000000 == 0) {
            fprintf(stderr, "[sample]   Processed %lu reads, kept %lu (%.1f%%)\n",
                    (unsigned long)total_reads, (unsigned long)kept_reads,
                    100.0 * kept_reads / total_reads);
        }
    }
    
    fprintf(stderr, "[sample] Processed %lu reads, kept %lu (%.1f%%)\n",
            (unsigned long)total_reads, (unsigned long)kept_reads,
            total_reads > 0 ? 100.0 * kept_reads / total_reads : 0.0);
    
    bam_destroy1(b);
    hts_itr_destroy(iter);
    
    /* Close files before sorting */
    hts_close(out_fp);
    out_fp = NULL;
    
    /* Sort and index output BAM */
    if (!args->no_sort && ret == 0) {
        fprintf(stderr, "[sample] Sorting output BAM...\n");
        
        /* Create temp file for sorted output */
        char tmp_path[1024];
        snprintf(tmp_path, sizeof(tmp_path), "%s.tmp.bam", args->out_bam);
        
        /* Use samtools sort via htslib */
        char sort_cmd[4096];
        snprintf(sort_cmd, sizeof(sort_cmd), 
                 "samtools sort -o '%s' '%s' && mv '%s' '%s' && samtools index '%s'",
                 tmp_path, args->out_bam, tmp_path, args->out_bam, args->out_bam);
        
        int sort_ret = system(sort_cmd);
        if (sort_ret != 0) {
            fprintf(stderr, "Warning: Failed to sort/index output BAM (exit code %d)\n", sort_ret);
            fprintf(stderr, "         The output BAM may be unsorted.\n");
        } else {
            fprintf(stderr, "[sample] Sorting and indexing complete.\n");
        }
    }
    
    /* Calculate and print metrics */
    if (ret == 0) {
        fprintf(stderr, "[sample] Computing metrics...\n");
        
        /* Get output depth */
        depth_array_t *output_depth = depth_from_bam(args->out_bam, region->contig,
                                                      region->start, region->end);
        if (output_depth) {
            metrics_result_t metrics;
            if (metrics_calculate(template_depth, output_depth, &metrics) == 0) {
                metrics_print(&metrics);
            }
            depth_array_free(output_depth);
        } else {
            fprintf(stderr, "Warning: Could not compute metrics\n");
        }
    }
    
cleanup:
    free(cumsum);
    depth_array_free(source_depth);
    depth_array_free(template_depth);
    hts_idx_destroy(idx);
    sam_hdr_destroy(header);
    hts_close(source_fp);
    if (out_fp) hts_close(out_fp);
    region_free(region);
    
    if (ret == 0) {
        fprintf(stderr, "[sample] Done. Output written to: %s\n", args->out_bam);
    }
    
    return ret;
}

