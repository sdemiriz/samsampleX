/*
 * bed.c - BED file reading and writing utilities
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "bed.h"
#include "samsampleX.h"

/*
 * Write a single BED entry to file
 */
void bed_write_entry(FILE *fp, const char *chrom, int64_t start, int64_t end, int32_t depth) {
    fprintf(fp, "%s\t%ld\t%ld\t%d\n", chrom, (long)start, (long)end, depth);
}

/*
 * Allocate a depth array
 */
depth_array_t *depth_array_alloc(const char *contig, int64_t start, int64_t end) {
    depth_array_t *arr;
    SAFE_MALLOC(arr, sizeof(depth_array_t));
    
    arr->length = (size_t)(end - start);
    arr->start = start;
    arr->end = end;
    arr->contig = strdup(contig);
    
    SAFE_CALLOC(arr->depths, arr->length, sizeof(int32_t));
    
    return arr;
}

/*
 * Free a depth array
 */
void depth_array_free(depth_array_t *arr) {
    if (arr) {
        free(arr->depths);
        free(arr->contig);
        free(arr);
    }
}

/*
 * Read depth values from a BED file into a depth array for a specific region.
 */
depth_array_t *bed_read_depths(const char *bed_path, const char *contig,
                               int64_t region_start, int64_t region_end) {
    FILE *fp = fopen(bed_path, "r");
    if (!fp) {
        fprintf(stderr, "Error: Cannot open BED file: %s\n", bed_path);
        return NULL;
    }
    
    depth_array_t *arr = depth_array_alloc(contig, region_start, region_end);
    
    char line[4096];
    while (fgets(line, sizeof(line), fp)) {
        /* Skip comment and empty lines */
        if (line[0] == '#' || line[0] == '\n' || line[0] == '\0') {
            continue;
        }
        
        /* Parse BED4 format: chrom, start, end, depth */
        char chrom[256];
        int64_t bed_start, bed_end;
        int32_t depth;
        
        if (sscanf(line, "%255s %ld %ld %d", chrom, &bed_start, &bed_end, &depth) != 4) {
            /* Try with tab delimiter explicitly */
            char *token = strtok(line, "\t");
            if (!token) continue;
            strncpy(chrom, token, sizeof(chrom) - 1);
            
            token = strtok(NULL, "\t");
            if (!token) continue;
            bed_start = atol(token);
            
            token = strtok(NULL, "\t");
            if (!token) continue;
            bed_end = atol(token);
            
            token = strtok(NULL, "\t\n");
            if (!token) continue;
            depth = atoi(token);
        }
        
        /* Skip if different contig */
        if (strcmp(chrom, contig) != 0) {
            continue;
        }
        
        /* Calculate overlap with target region */
        int64_t overlap_start = MAX(bed_start, region_start);
        int64_t overlap_end = MIN(bed_end, region_end);
        
        /* Skip if no overlap */
        if (overlap_start >= overlap_end) {
            continue;
        }
        
        /* Fill depth values in the array */
        size_t i_start = (size_t)(overlap_start - region_start);
        size_t i_end = (size_t)(overlap_end - region_start);
        
        for (size_t i = i_start; i < i_end; i++) {
            arr->depths[i] = depth;
        }
    }
    
    fclose(fp);
    return arr;
}

/*
 * Combine multiple depth arrays according to a mode.
 */
depth_array_t *bed_combine_depths(depth_array_t **arrays, size_t n_arrays,
                                  combine_mode_t mode, uint32_t seed) {
    if (n_arrays == 0 || !arrays || !arrays[0]) {
        return NULL;
    }
    
    /* If only one array, return a copy */
    if (n_arrays == 1) {
        depth_array_t *result = depth_array_alloc(arrays[0]->contig, 
                                                   arrays[0]->start, 
                                                   arrays[0]->end);
        memcpy(result->depths, arrays[0]->depths, result->length * sizeof(int32_t));
        return result;
    }
    
    /* Create result array with same dimensions as first input */
    depth_array_t *result = depth_array_alloc(arrays[0]->contig,
                                               arrays[0]->start,
                                               arrays[0]->end);
    
    /* Seed random number generator for MODE_RANDOM */
    if (mode == MODE_RANDOM) {
        srand(seed);
    }
    
    /* Combine depths position by position */
    for (size_t i = 0; i < result->length; i++) {
        int32_t min_val = arrays[0]->depths[i];
        int32_t max_val = arrays[0]->depths[i];
        int64_t sum = arrays[0]->depths[i];
        
        for (size_t j = 1; j < n_arrays; j++) {
            int32_t val = arrays[j]->depths[i];
            if (val < min_val) min_val = val;
            if (val > max_val) max_val = val;
            sum += val;
        }
        
        switch (mode) {
            case MODE_MIN:
                result->depths[i] = min_val;
                break;
            case MODE_MAX:
                result->depths[i] = max_val;
                break;
            case MODE_MEAN:
                result->depths[i] = (int32_t)(sum / (int64_t)n_arrays);
                break;
            case MODE_RANDOM:
                if (min_val == max_val) {
                    result->depths[i] = min_val;
                } else {
                    result->depths[i] = min_val + (rand() % (max_val - min_val + 1));
                }
                break;
        }
    }
    
    return result;
}

