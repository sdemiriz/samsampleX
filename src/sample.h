/*
 * sample.h - Sample subcommand for depth-aware BAM sampling
 */

#ifndef SAMPLE_H
#define SAMPLE_H

#include "samsampleX.h"

/*
 * Maximum number of template BED files
 */
#define MAX_TEMPLATE_BEDS 64

/*
 * Sample command arguments
 */
typedef struct {
    const char *source_bam;                     /* Input BAM file to sample from */
    const char *template_beds[MAX_TEMPLATE_BEDS]; /* Template BED files */
    int n_template_beds;                        /* Number of template BED files */
    const char *region;                         /* Target region */
    const char *out_bam;                        /* Output BAM file */
    combine_mode_t mode;                        /* How to combine multiple templates */
    uint32_t seed;                              /* Random seed */
    int no_sort;                                /* Skip sorting/indexing if true */
    int no_metrics;                             /* Skip metrics calculation if true */
} sample_args_t;

/*
 * Run the sample subcommand.
 *
 * @args: Command arguments
 * @return: 0 on success, non-zero on error
 */
int sample_run(sample_args_t *args);

/*
 * Print sample subcommand usage.
 */
void sample_usage(void);

#endif /* SAMPLE_H */

