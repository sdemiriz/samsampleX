/*
 * map.h - Map subcommand for extracting depth from BAM to BED
 */

#ifndef MAP_H
#define MAP_H

/*
 * Map command arguments
 */
typedef struct {
    const char *template_bam;   /* Input BAM file */
    const char *region;         /* Target region (samtools-style) */
    const char *out_bed;        /* Output BED file */
    int collapse;               /* Collapse tolerance (0 = no collapse) */
} map_args_t;

/*
 * Run the map subcommand.
 *
 * @args: Command arguments
 * @return: 0 on success, non-zero on error
 */
int map_run(map_args_t *args);

/*
 * Print map subcommand usage.
 */
void map_usage(void);

#endif /* MAP_H */

