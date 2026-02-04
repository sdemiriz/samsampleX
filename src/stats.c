/*
 * stats.c - Statistics subcommand for comparing two BAM files
 *
 * Computes depth distributions from two BAM files and calculates
 * comparison metrics (Wasserstein distance, Total Variation).
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <htslib/sam.h>
#include <htslib/hts.h>

#include "stats.h"
#include "depth.h"
#include "bed.h"
#include "metrics.h"
#include "samsampleX.h"

/*
 * Print usage for 'stats' subcommand.
 */
void stats_usage(void) {
    fprintf(stderr, "Usage: %s stats [options]\n\n", SAMSAMPLEX_NAME);
    fprintf(stderr, "Compare depth distributions between two BAM files.\n\n");
    fprintf(stderr, "Required:\n");
    fprintf(stderr, "  -a, --bam-a FILE      First BAM file (e.g., template/reference)\n");
    fprintf(stderr, "  -b, --bam-b FILE      Second BAM file (e.g., sampled output)\n");
    fprintf(stderr, "  -r, --region REGION   Region to compare (e.g., chr1:1000-2000)\n\n");
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "  -h, --help            Show this help message\n\n");
    fprintf(stderr, "Output metrics:\n");
    fprintf(stderr, "  - Mean depth for each BAM\n");
    fprintf(stderr, "  - Total Variation (TV) distance\n");
    fprintf(stderr, "  - Wasserstein-1 distance (normalized)\n");
}

/*
 * Run the 'stats' subcommand.
 */
int stats_run(const stats_args_t *args) {
    /* Validate required arguments */
    if (!args->bam_a) {
        fprintf(stderr, "Error: --bam-a is required\n");
        stats_usage();
        return 1;
    }
    if (!args->bam_b) {
        fprintf(stderr, "Error: --bam-b is required\n");
        stats_usage();
        return 1;
    }
    if (!args->region) {
        fprintf(stderr, "Error: --region is required\n");
        stats_usage();
        return 1;
    }
    
    /* Parse region */
    region_t *region = region_parse(args->region);
    if (!region || !region->contig) {
        fprintf(stderr, "Error: Invalid region format: %s\n", args->region);
        return 1;
    }
    
    fprintf(stderr, "Computing depth for BAM A: %s\n", args->bam_a);
    fprintf(stderr, "Computing depth for BAM B: %s\n", args->bam_b);
    fprintf(stderr, "Region: %s:%ld-%ld\n", region->contig, 
            region->start + 1, region->end);  /* Print 1-based for user */
    
    /* Compute depth for first BAM */
    depth_array_t *depth_a = depth_from_bam(args->bam_a, region->contig,
                                             region->start, region->end);
    if (!depth_a) {
        fprintf(stderr, "Error: Failed to compute depth for BAM A\n");
        region_free(region);
        return 1;
    }
    
    /* Compute depth for second BAM */
    depth_array_t *depth_b = depth_from_bam(args->bam_b, region->contig,
                                             region->start, region->end);
    if (!depth_b) {
        fprintf(stderr, "Error: Failed to compute depth for BAM B\n");
        depth_array_free(depth_a);
        region_free(region);
        return 1;
    }
    
    /* Calculate and print metrics */
    metrics_result_t result;
    if (metrics_calculate(depth_a, depth_b, &result) != 0) {
        fprintf(stderr, "Error: Failed to calculate metrics\n");
        depth_array_free(depth_a);
        depth_array_free(depth_b);
        region_free(region);
        return 1;
    }
    
    /* Print with labels indicating which BAM is which */
    fprintf(stderr, "\n========== Comparison Metrics ==========\n");
    fprintf(stderr, "BAM A (reference):       %s\n", args->bam_a);
    fprintf(stderr, "BAM B (comparison):      %s\n", args->bam_b);
    fprintf(stderr, "Region:                  %s:%ld-%ld\n", region->contig,
            region->start + 1, region->end);
    fprintf(stderr, "Positions compared:      %zu\n", depth_a->length);
    fprintf(stderr, "-----------------------------------------\n");
    fprintf(stderr, "Mean Depth (BAM A):      %.2f\n", result.mean_template);
    fprintf(stderr, "Mean Depth (BAM B):      %.2f\n", result.mean_output);
    fprintf(stderr, "Total Variation:         %.4f\n", result.tv);
    fprintf(stderr, "Wasserstein-1 (norm):    %.6f\n", result.wasserstein);
    fprintf(stderr, "=========================================\n");
    
    /* Cleanup */
    depth_array_free(depth_a);
    depth_array_free(depth_b);
    region_free(region);
    
    return 0;
}
