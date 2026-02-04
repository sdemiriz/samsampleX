/*
 * depth.c - Depth array computation from BAM files using htslib
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include <htslib/sam.h>
#include <htslib/hts.h>

#include "depth.h"
#include "bed.h"
#include "samsampleX.h"

/*
 * Check if two contig names match, considering chr prefix variations.
 * Returns 1 (true) if names match with or without "chr" prefix.
 * E.g., "chr21" matches "21", "chrX" matches "X".
 */
int contig_names_match(const char *name1, const char *name2) {
    if (!name1 || !name2) return 0;
    
    /* Exact match */
    if (strcmp(name1, name2) == 0) return 1;
    
    /* Check if one has "chr" prefix and the other doesn't */
    const char *n1_base = name1;
    const char *n2_base = name2;
    
    if (strncmp(name1, "chr", 3) == 0) {
        n1_base = name1 + 3;
    }
    if (strncmp(name2, "chr", 3) == 0) {
        n2_base = name2 + 3;
    }
    
    return strcmp(n1_base, n2_base) == 0;
}

/*
 * Resolve contig name against BAM header, handling chr prefix mismatch.
 * Tries: exact match, then with "chr" removed, then with "chr" added.
 * Returns the matching contig name from header (caller must free), or NULL.
 */
char *resolve_contig_name(sam_hdr_t *header, const char *contig) {
    if (!header || !contig) return NULL;
    
    /* Try exact match first */
    int tid = sam_hdr_name2tid(header, contig);
    if (tid >= 0) {
        return strdup(contig);
    }
    
    /* Try without "chr" prefix if it has one */
    if (strncmp(contig, "chr", 3) == 0) {
        const char *no_chr = contig + 3;
        tid = sam_hdr_name2tid(header, no_chr);
        if (tid >= 0) {
            fprintf(stderr, "Note: Using contig '%s' (matched from '%s')\n", no_chr, contig);
            return strdup(no_chr);
        }
    }
    
    /* Try with "chr" prefix if it doesn't have one */
    if (strncmp(contig, "chr", 3) != 0) {
        char with_chr[256];
        snprintf(with_chr, sizeof(with_chr), "chr%s", contig);
        tid = sam_hdr_name2tid(header, with_chr);
        if (tid >= 0) {
            fprintf(stderr, "Note: Using contig '%s' (matched from '%s')\n", with_chr, contig);
            return strdup(with_chr);
        }
    }
    
    /* No match found */
    return NULL;
}

/*
 * Compute depth of coverage for a region from a BAM file.
 */
depth_array_t *depth_from_bam(const char *bam_path, const char *contig,
                              int64_t start, int64_t end) {
    /* Open BAM file */
    htsFile *fp = hts_open(bam_path, "r");
    if (!fp) {
        fprintf(stderr, "Error: Cannot open BAM file: %s\n", bam_path);
        return NULL;
    }
    
    /* Read header */
    sam_hdr_t *header = sam_hdr_read(fp);
    if (!header) {
        fprintf(stderr, "Error: Cannot read BAM header: %s\n", bam_path);
        hts_close(fp);
        return NULL;
    }
    
    /* Load index */
    hts_idx_t *idx = sam_index_load(fp, bam_path);
    if (!idx) {
        fprintf(stderr, "Error: Cannot load BAM index for: %s\n", bam_path);
        fprintf(stderr, "       Make sure the BAM file is indexed (.bai file exists)\n");
        sam_hdr_destroy(header);
        hts_close(fp);
        return NULL;
    }
    
    /* Resolve contig name (handles chr prefix mismatch) */
    char *resolved_contig = resolve_contig_name(header, contig);
    if (!resolved_contig) {
        fprintf(stderr, "Error: Contig '%s' not found in BAM header\n", contig);
        hts_idx_destroy(idx);
        sam_hdr_destroy(header);
        hts_close(fp);
        return NULL;
    }
    
    /* Get contig ID using resolved name */
    int tid = sam_hdr_name2tid(header, resolved_contig);
    
    /* Get contig length if end is not specified */
    if (end < 0 || end > sam_hdr_tid2len(header, tid)) {
        end = sam_hdr_tid2len(header, tid);
    }
    if (start < 0) {
        start = 0;
    }
    
    /* Allocate depth array with resolved contig name */
    depth_array_t *arr = depth_array_alloc(resolved_contig, start, end);
    free(resolved_contig);
    
    /* Set up iterator for region */
    hts_itr_t *iter = sam_itr_queryi(idx, tid, start, end);
    if (!iter) {
        fprintf(stderr, "Error: Cannot create iterator for region\n");
        depth_array_free(arr);
        hts_idx_destroy(idx);
        sam_hdr_destroy(header);
        hts_close(fp);
        return NULL;
    }
    
    /* Set up pileup */
    bam_plp_t plp = bam_plp_init(NULL, NULL);
    bam_plp_set_maxcnt(plp, INT_MAX);  /* No limit on depth */
    
    /* Read all alignments in region and compute depth using pileup */
    bam1_t *b = bam_init1();
    int ret;
    
    /* We'll use a simpler approach: iterate through reads and increment depth */
    while ((ret = sam_itr_next(fp, iter, b)) >= 0) {
        /* Skip unmapped, secondary, QC failed, and duplicate reads */
        if (b->core.flag & (BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP)) {
            continue;
        }
        
        /* Get read alignment positions */
        int64_t read_start = b->core.pos;
        int64_t read_end = bam_endpos(b);
        
        /* Clip to region */
        int64_t overlap_start = MAX(read_start, start);
        int64_t overlap_end = MIN(read_end, end);
        
        if (overlap_start >= overlap_end) {
            continue;
        }
        
        /* Increment depth for each covered position */
        /* Note: This is a simplified approach. For more accurate depth
         * that accounts for CIGAR operations (insertions, deletions, etc.),
         * we would need to walk through the CIGAR string. */
        for (int64_t pos = overlap_start; pos < overlap_end; pos++) {
            size_t idx_pos = (size_t)(pos - start);
            if (idx_pos < arr->length) {
                arr->depths[idx_pos]++;
            }
        }
    }
    
    /* Cleanup */
    bam_destroy1(b);
    bam_plp_destroy(plp);
    hts_itr_destroy(iter);
    hts_idx_destroy(idx);
    sam_hdr_destroy(header);
    hts_close(fp);
    
    return arr;
}

/*
 * Get contig length from BAM header.
 * Handles chr prefix mismatch automatically.
 */
int64_t get_contig_length(sam_hdr_t *header, const char *contig) {
    char *resolved = resolve_contig_name(header, contig);
    if (!resolved) {
        return -1;
    }
    int tid = sam_hdr_name2tid(header, resolved);
    free(resolved);
    if (tid < 0) {
        return -1;
    }
    return sam_hdr_tid2len(header, tid);
}

/*
 * Parse a region string (samtools-style) into components.
 */
region_t *region_parse(const char *region_str) {
    region_t *region;
    SAFE_MALLOC(region, sizeof(region_t));
    
    /* Make a copy to work with */
    char *str = strdup(region_str);
    if (!str) {
        free(region);
        return NULL;
    }
    
    /* Initialize defaults */
    region->start = -1;
    region->end = -1;
    region->contig = NULL;
    
    /* Find colon separator */
    char *colon = strchr(str, ':');
    if (colon) {
        *colon = '\0';
        region->contig = strdup(str);
        
        /* Parse start-end */
        char *dash = strchr(colon + 1, '-');
        if (dash) {
            *dash = '\0';
            region->start = atol(colon + 1);
            region->end = atol(dash + 1);
        } else {
            /* Only start position provided */
            region->start = atol(colon + 1);
            region->end = -1;  /* Will be set later */
        }
        
        /* Convert from 1-based to 0-based coordinates */
        if (region->start > 0) {
            region->start--;
        }
        /* end stays as-is since it's exclusive in 0-based */
    } else {
        /* Just contig name, no coordinates */
        region->contig = strdup(str);
    }
    
    free(str);
    return region;
}

/*
 * Free a region struct.
 */
void region_free(region_t *region) {
    if (region) {
        free(region->contig);
        free(region);
    }
}

/*
 * Parse combine mode string to enum.
 */
combine_mode_t parse_combine_mode(const char *mode_str) {
    if (strcmp(mode_str, "min") == 0) return MODE_MIN;
    if (strcmp(mode_str, "max") == 0) return MODE_MAX;
    if (strcmp(mode_str, "mean") == 0) return MODE_MEAN;
    if (strcmp(mode_str, "random") == 0) return MODE_RANDOM;
    
    fprintf(stderr, "Warning: Unknown mode '%s', using 'min'\n", mode_str);
    return MODE_MIN;
}

/*
 * Convert combine mode enum to string.
 */
const char *combine_mode_to_string(combine_mode_t mode) {
    switch (mode) {
        case MODE_MIN: return "min";
        case MODE_MAX: return "max";
        case MODE_MEAN: return "mean";
        case MODE_RANDOM: return "random";
        default: return "unknown";
    }
}

