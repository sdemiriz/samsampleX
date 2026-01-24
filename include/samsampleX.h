/*
 * samsampleX.h - Shared types and constants for samsampleX
 *
 * A tool for depth-aware BAM file sampling.
 */

#ifndef SAMSAMPLEX_H
#define SAMSAMPLEX_H

#include <stdint.h>
#include <stddef.h>
#include <stdbool.h>

/* Version info */
#define SAMSAMPLEX_VERSION "0.1.0"
#define SAMSAMPLEX_NAME "samsampleX"

/* Default values */
#define DEFAULT_COLLAPSE 0
#define DEFAULT_SEED 42
#define DEFAULT_OUT_BED "out.bed"
#define DEFAULT_OUT_BAM "out.bam"
#define DEFAULT_MODE MODE_MIN

/* Combine modes for multiple BED templates */
typedef enum {
    MODE_MIN = 0,
    MODE_MAX,
    MODE_MEAN,
    MODE_RANDOM
} combine_mode_t;

/* Region specification (parsed from samtools-style string) */
typedef struct {
    char *contig;       /* Chromosome/contig name */
    int64_t start;      /* 0-based start position (-1 if whole contig) */
    int64_t end;        /* 0-based end position (exclusive, -1 if whole contig) */
} region_t;

/* Depth array for a genomic region */
typedef struct {
    int32_t *depths;    /* Array of depth values */
    size_t length;      /* Number of positions */
    char *contig;       /* Contig name */
    int64_t start;      /* Start position */
    int64_t end;        /* End position */
} depth_array_t;

/* BED entry for depth template */
typedef struct {
    char *chrom;
    int64_t start;
    int64_t end;
    int32_t depth;
} bed_entry_t;

/* Utility macros */
#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define MAX(a, b) ((a) > (b) ? (a) : (b))

/* Memory allocation with error checking */
#define SAFE_MALLOC(ptr, size) do { \
    (ptr) = malloc(size); \
    if (!(ptr)) { \
        fprintf(stderr, "Error: Memory allocation failed\n"); \
        exit(EXIT_FAILURE); \
    } \
} while(0)

#define SAFE_CALLOC(ptr, count, size) do { \
    (ptr) = calloc(count, size); \
    if (!(ptr)) { \
        fprintf(stderr, "Error: Memory allocation failed\n"); \
        exit(EXIT_FAILURE); \
    } \
} while(0)

#define SAFE_REALLOC(ptr, size) do { \
    void *_tmp = realloc(ptr, size); \
    if (!_tmp) { \
        fprintf(stderr, "Error: Memory reallocation failed\n"); \
        exit(EXIT_FAILURE); \
    } \
    (ptr) = _tmp; \
} while(0)

/* Function prototypes for region parsing */
region_t *region_parse(const char *region_str);
void region_free(region_t *region);

/* Convert combine mode string to enum */
combine_mode_t parse_combine_mode(const char *mode_str);
const char *combine_mode_to_string(combine_mode_t mode);

#endif /* SAMSAMPLEX_H */

