/*
 * bed.h - BED file reading and writing utilities
 */

#ifndef BED_H
#define BED_H

#include <stdio.h>
#include "samsampleX.h"

/*
 * Write a single BED entry to file
 */
void bed_write_entry(FILE *fp, const char *chrom, int64_t start, int64_t end, int32_t depth);

/*
 * Read depth values from a BED file into a depth array for a specific region.
 * The returned array covers positions [region_start, region_end).
 * Positions not covered by BED entries are set to 0.
 *
 * @bed_path: Path to BED file
 * @contig: Target contig name
 * @region_start: Start of region (0-based)
 * @region_end: End of region (exclusive)
 * @return: Allocated depth array (caller must free with depth_array_free)
 */
depth_array_t *bed_read_depths(const char *bed_path, const char *contig,
                               int64_t region_start, int64_t region_end);

/*
 * Combine multiple depth arrays according to a mode.
 *
 * @arrays: Array of depth_array_t pointers
 * @n_arrays: Number of arrays
 * @mode: Combine mode (min, max, mean, random)
 * @seed: Random seed (used only for MODE_RANDOM)
 * @return: Combined depth array (caller must free)
 */
depth_array_t *bed_combine_depths(depth_array_t **arrays, size_t n_arrays,
                                  combine_mode_t mode, uint32_t seed);

/*
 * Free a depth array
 */
void depth_array_free(depth_array_t *arr);

/*
 * Allocate a depth array
 */
depth_array_t *depth_array_alloc(const char *contig, int64_t start, int64_t end);

#endif /* BED_H */

