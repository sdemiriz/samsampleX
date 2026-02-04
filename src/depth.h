/*
 * depth.h - Depth array computation from BAM files using htslib
 */

#ifndef DEPTH_H
#define DEPTH_H

#include <htslib/sam.h>
#include <htslib/hts.h>
#include "samsampleX.h"

/*
 * Compute depth of coverage for a region from a BAM file.
 *
 * @bam_path: Path to indexed BAM file
 * @contig: Target contig name
 * @start: Start position (0-based)
 * @end: End position (exclusive)
 * @return: Depth array (caller must free with depth_array_free)
 */
depth_array_t *depth_from_bam(const char *bam_path, const char *contig,
                              int64_t start, int64_t end);

/*
 * Get contig length from BAM header.
 *
 * @header: BAM header
 * @contig: Contig name
 * @return: Contig length, or -1 if not found
 */
int64_t get_contig_length(sam_hdr_t *header, const char *contig);

/*
 * Parse a region string (samtools-style) into components.
 * Supports formats: "chr1", "chr1:1000", "chr1:1000-2000"
 *
 * @region_str: Region string
 * @return: Parsed region (caller must free with region_free)
 */
region_t *region_parse(const char *region_str);

/*
 * Free a region struct.
 */
void region_free(region_t *region);

/*
 * Parse combine mode string to enum.
 */
combine_mode_t parse_combine_mode(const char *mode_str);

/*
 * Convert combine mode enum to string.
 */
const char *combine_mode_to_string(combine_mode_t mode);

/*
 * Resolve contig name against BAM header, handling chr prefix mismatch.
 * Tries: exact match, then with "chr" removed, then with "chr" added.
 *
 * @header: BAM header to search
 * @contig: User-provided contig name
 * @return: Matched contig name from header (caller must free), or NULL if not found
 */
char *resolve_contig_name(sam_hdr_t *header, const char *contig);

/*
 * Check if two contig names match, considering chr prefix variations.
 * E.g., "chr21" matches "21", "chrX" matches "X".
 *
 * @name1: First contig name
 * @name2: Second contig name
 * @return: 1 if names match, 0 otherwise
 */
int contig_names_match(const char *name1, const char *name2);

#endif /* DEPTH_H */

