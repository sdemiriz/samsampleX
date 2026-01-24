/*
 * map.c - Map subcommand for extracting depth from BAM to BED
 *
 * This module extracts depth of coverage from a template BAM file
 * and writes it to a BED file. The output can be "collapsed" to
 * merge consecutive positions with similar depths.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>

#include <htslib/sam.h>
#include <htslib/hts.h>

#include "map.h"
#include "bed.h"
#include "depth.h"
#include "samsampleX.h"

/*
 * Print map subcommand usage.
 */
void map_usage(void) {
    fprintf(stderr, "Usage: %s map [options]\n\n", SAMSAMPLEX_NAME);
    fprintf(stderr, "Extract depth of coverage from a BAM file and write to BED format.\n\n");
    fprintf(stderr, "Required options:\n");
    fprintf(stderr, "  --template-bam FILE   Input BAM file to extract depth from\n");
    fprintf(stderr, "  --region REGION       Target region (samtools-style: chr1, chr1:1000-2000)\n\n");
    fprintf(stderr, "Optional:\n");
    fprintf(stderr, "  --out-bed FILE        Output BED file [default: %s]\n", DEFAULT_OUT_BED);
    fprintf(stderr, "  --collapse INT        Merge consecutive positions with depth difference <= INT\n");
    fprintf(stderr, "                        Use 0 for per-base output [default: %d]\n", DEFAULT_COLLAPSE);
    fprintf(stderr, "  -h, --help            Show this help message\n");
}

/*
 * Write depth array to BED file with optional collapsing.
 */
static int write_bed_output(FILE *fp, depth_array_t *arr, int collapse) {
    if (!arr || arr->length == 0) {
        return 0;
    }
    
    if (collapse == 0) {
        /* Per-base output: one line per position */
        for (size_t i = 0; i < arr->length; i++) {
            int64_t pos = arr->start + (int64_t)i;
            bed_write_entry(fp, arr->contig, pos, pos + 1, arr->depths[i]);
        }
    } else {
        /* Collapsed output: merge consecutive positions with similar depth */
        int64_t interval_start = arr->start;
        int32_t interval_depth = arr->depths[0];
        
        for (size_t i = 1; i < arr->length; i++) {
            int32_t current_depth = arr->depths[i];
            
            /* Check if depth changed by more than collapse threshold */
            int32_t diff = current_depth - interval_depth;
            if (diff < 0) diff = -diff;  /* abs */
            
            if (diff > collapse) {
                /* Write previous interval */
                int64_t interval_end = arr->start + (int64_t)i;
                bed_write_entry(fp, arr->contig, interval_start, interval_end, interval_depth);
                
                /* Start new interval */
                interval_start = interval_end;
                interval_depth = current_depth;
            }
        }
        
        /* Write final interval */
        bed_write_entry(fp, arr->contig, interval_start, arr->end, interval_depth);
    }
    
    return 0;
}

/*
 * Run the map subcommand.
 */
int map_run(map_args_t *args) {
    /* Validate required arguments */
    if (!args->template_bam) {
        fprintf(stderr, "Error: --template-bam is required\n");
        map_usage();
        return 1;
    }
    if (!args->region) {
        fprintf(stderr, "Error: --region is required\n");
        map_usage();
        return 1;
    }
    if (args->collapse < 0) {
        fprintf(stderr, "Error: --collapse must be non-negative\n");
        return 1;
    }
    
    /* Set default output if not specified */
    if (!args->out_bed) {
        args->out_bed = DEFAULT_OUT_BED;
    }
    
    fprintf(stderr, "[map] Template BAM: %s\n", args->template_bam);
    fprintf(stderr, "[map] Region: %s\n", args->region);
    fprintf(stderr, "[map] Collapse: %d\n", args->collapse);
    fprintf(stderr, "[map] Output BED: %s\n", args->out_bed);
    
    /* Parse region */
    region_t *region = region_parse(args->region);
    if (!region || !region->contig) {
        fprintf(stderr, "Error: Invalid region format: %s\n", args->region);
        return 1;
    }
    
    /* If start/end not specified, get contig length from BAM */
    if (region->start < 0 || region->end < 0) {
        htsFile *fp = hts_open(args->template_bam, "r");
        if (!fp) {
            fprintf(stderr, "Error: Cannot open BAM file: %s\n", args->template_bam);
            region_free(region);
            return 1;
        }
        
        sam_hdr_t *header = sam_hdr_read(fp);
        if (!header) {
            fprintf(stderr, "Error: Cannot read BAM header\n");
            hts_close(fp);
            region_free(region);
            return 1;
        }
        
        int64_t contig_len = get_contig_length(header, region->contig);
        if (contig_len < 0) {
            fprintf(stderr, "Error: Contig '%s' not found in BAM\n", region->contig);
            sam_hdr_destroy(header);
            hts_close(fp);
            region_free(region);
            return 1;
        }
        
        if (region->start < 0) region->start = 0;
        if (region->end < 0) region->end = contig_len;
        
        sam_hdr_destroy(header);
        hts_close(fp);
    }
    
    fprintf(stderr, "[map] Parsed region: %s:%ld-%ld\n", 
            region->contig, (long)region->start, (long)region->end);
    
    /* Compute depth array from BAM */
    fprintf(stderr, "[map] Computing depth array (this may take a while)...\n");
    depth_array_t *depth = depth_from_bam(args->template_bam, region->contig,
                                          region->start, region->end);
    if (!depth) {
        fprintf(stderr, "Error: Failed to compute depth array\n");
        region_free(region);
        return 1;
    }
    
    fprintf(stderr, "[map] Computed depth for %zu positions\n", depth->length);
    
    /* Open output file */
    FILE *out_fp = fopen(args->out_bed, "w");
    if (!out_fp) {
        fprintf(stderr, "Error: Cannot open output file: %s\n", args->out_bed);
        depth_array_free(depth);
        region_free(region);
        return 1;
    }
    
    /* Write BED output */
    fprintf(stderr, "[map] Writing BED file%s...\n", 
            args->collapse > 0 ? " (collapsed)" : "");
    int ret = write_bed_output(out_fp, depth, args->collapse);
    
    fclose(out_fp);
    depth_array_free(depth);
    region_free(region);
    
    if (ret == 0) {
        fprintf(stderr, "[map] Done. Output written to: %s\n", args->out_bed);
    }
    
    return ret;
}

