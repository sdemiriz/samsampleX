/*
 * stats.h - Statistics subcommand for comparing two BAM files
 */

#ifndef STATS_H
#define STATS_H

#include "samsampleX.h"

/*
 * Arguments for the 'stats' subcommand.
 */
typedef struct {
    const char *bam_a;      /* First BAM file (e.g., template) */
    const char *bam_b;      /* Second BAM file (e.g., output from another tool) */
    const char *region;     /* Region to compare (samtools-style) */
} stats_args_t;

/*
 * Print usage for 'stats' subcommand.
 */
void stats_usage(void);

/*
 * Run the 'stats' subcommand.
 *
 * Computes depth arrays for both BAM files and calculates comparison metrics.
 *
 * @args: Parsed command-line arguments
 * @return: 0 on success, non-zero on error
 */
int stats_run(const stats_args_t *args);

#endif /* STATS_H */
