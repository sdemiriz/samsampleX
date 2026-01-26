/*
 * plot.h - Plot subcommand for visualizing depth of coverage
 */

#ifndef PLOT_H
#define PLOT_H

#include "samsampleX.h"

/*
 * Default output filenames
 */
#define DEFAULT_OUT_PNG "plot.png"
#define DEFAULT_OUT_TSV "depths.tsv"

/*
 * Default plot dimensions (pixels)
 */
#define DEFAULT_PLOT_WIDTH 1200.0
#define DEFAULT_PLOT_HEIGHT 600.0

/*
 * Plot command arguments
 */
typedef struct {
    const char *source_bam;     /* Source BAM file (required) */
    const char *template_bam;   /* Template BAM file (mutually exclusive with template_bed) */
    const char *template_bed;   /* Template BED file (mutually exclusive with template_bam) */
    const char *out_bam;        /* Output BAM file (required) */
    const char *region;         /* Target region (required) */
    const char *out_png;        /* Output PNG file (mutually exclusive with out_tsv) */
    const char *out_tsv;        /* Output TSV file (mutually exclusive with out_png) */
} plot_args_t;

/*
 * Run the plot subcommand.
 *
 * @args: Command arguments
 * @return: 0 on success, non-zero on error
 */
int plot_run(plot_args_t *args);

/*
 * Print plot subcommand usage.
 */
void plot_usage(void);

#endif /* PLOT_H */

