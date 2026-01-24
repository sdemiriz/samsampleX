/*
 * main.c - Entry point for samsampleX
 *
 * samsampleX is a tool for depth-aware BAM file sampling.
 * It supports multiple subcommands:
 *   - map: Extract depth from BAM to BED template
 *   - sample: Sample reads to match template depth distribution
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>

#include "samsampleX.h"
#include "map.h"
#include "sample.h"
#include "depth.h"

/*
 * Print main usage message.
 */
static void print_usage(void) {
    fprintf(stderr, "samsampleX v%s - Depth-aware BAM file sampling\n\n", SAMSAMPLEX_VERSION);
    fprintf(stderr, "Usage: %s <command> [options]\n\n", SAMSAMPLEX_NAME);
    fprintf(stderr, "Commands:\n");
    fprintf(stderr, "  map      Extract depth of coverage from BAM to BED template\n");
    fprintf(stderr, "  sample   Sample reads from BAM to match template depth distribution\n\n");
    fprintf(stderr, "Use '%s <command> --help' for command-specific help.\n", SAMSAMPLEX_NAME);
}

/*
 * Parse arguments for the 'map' subcommand.
 */
static int parse_map_args(int argc, char *argv[], map_args_t *args) {
    static struct option long_options[] = {
        {"template-bam", required_argument, 0, 't'},
        {"region",       required_argument, 0, 'r'},
        {"out-bed",      required_argument, 0, 'o'},
        {"collapse",     required_argument, 0, 'c'},
        {"help",         no_argument,       0, 'h'},
        {0, 0, 0, 0}
    };
    
    /* Set defaults */
    memset(args, 0, sizeof(*args));
    args->collapse = DEFAULT_COLLAPSE;
    args->out_bed = DEFAULT_OUT_BED;
    
    int opt;
    int option_index = 0;
    
    /* Reset getopt for subcommand parsing */
    optind = 1;
    
    while ((opt = getopt_long(argc, argv, "t:r:o:c:h", long_options, &option_index)) != -1) {
        switch (opt) {
            case 't':
                args->template_bam = optarg;
                break;
            case 'r':
                args->region = optarg;
                break;
            case 'o':
                args->out_bed = optarg;
                break;
            case 'c':
                args->collapse = atoi(optarg);
                break;
            case 'h':
                map_usage();
                exit(0);
            default:
                map_usage();
                return -1;
        }
    }
    
    return 0;
}

/*
 * Parse arguments for the 'sample' subcommand.
 * 
 * Note: --template-bed accepts multiple files: --template-bed file1.bed file2.bed ...
 * Files are collected until the next option (starting with -) or end of arguments.
 */
static int parse_sample_args(int argc, char *argv[], sample_args_t *args) {
    static struct option long_options[] = {
        {"source-bam",   required_argument, 0, 's'},
        {"template-bed", required_argument, 0, 't'},
        {"region",       required_argument, 0, 'r'},
        {"out-bam",      required_argument, 0, 'o'},
        {"mode",         required_argument, 0, 'm'},
        {"seed",         required_argument, 0, 'S'},
        {"no-sort",      no_argument,       0, 'n'},
        {"no-metrics",   no_argument,       0, 'M'},
        {"help",         no_argument,       0, 'h'},
        {0, 0, 0, 0}
    };
    
    /* Set defaults */
    memset(args, 0, sizeof(*args));
    args->mode = DEFAULT_MODE;
    args->seed = DEFAULT_SEED;
    args->out_bam = DEFAULT_OUT_BAM;
    args->no_sort = 0;
    args->no_metrics = 0;
    
    int opt;
    int option_index = 0;
    
    /* Reset getopt for subcommand parsing */
    optind = 1;
    
    while ((opt = getopt_long(argc, argv, "s:t:r:o:m:S:nMh", long_options, &option_index)) != -1) {
        switch (opt) {
            case 's':
                args->source_bam = optarg;
                break;
            case 't':
                /* First BED file from optarg */
                if (args->n_template_beds > 0) {
                    fprintf(stderr, "Error: --template-bed can only be specified once\n");
                    fprintf(stderr, "       Use: --template-bed file1.bed file2.bed ...\n");
                    return -1;
                }
                if (args->n_template_beds >= MAX_TEMPLATE_BEDS) {
                    fprintf(stderr, "Error: Too many template BED files (max %d)\n", MAX_TEMPLATE_BEDS);
                    return -1;
                }
                args->template_beds[args->n_template_beds++] = optarg;
                
                /* Collect additional BED files until next option or end */
                while (optind < argc && argv[optind][0] != '-') {
                    if (args->n_template_beds >= MAX_TEMPLATE_BEDS) {
                        fprintf(stderr, "Error: Too many template BED files (max %d)\n", MAX_TEMPLATE_BEDS);
                        return -1;
                    }
                    args->template_beds[args->n_template_beds++] = argv[optind];
                    optind++;
                }
                break;
            case 'r':
                args->region = optarg;
                break;
            case 'o':
                args->out_bam = optarg;
                break;
            case 'm':
                args->mode = parse_combine_mode(optarg);
                break;
            case 'S':
                args->seed = (uint32_t)atoi(optarg);
                break;
            case 'n':
                args->no_sort = 1;
                break;
            case 'M':
                args->no_metrics = 1;
                break;
            case 'h':
                sample_usage();
                exit(0);
            default:
                sample_usage();
                return -1;
        }
    }
    
    return 0;
}

/*
 * Main entry point.
 */
int main(int argc, char *argv[]) {
    if (argc < 2) {
        print_usage();
        return 1;
    }
    
    const char *command = argv[1];
    
    /* Check for help flag */
    if (strcmp(command, "-h") == 0 || strcmp(command, "--help") == 0) {
        print_usage();
        return 0;
    }
    
    /* Check for version flag */
    if (strcmp(command, "-v") == 0 || strcmp(command, "--version") == 0) {
        printf("samsampleX v%s\n", SAMSAMPLEX_VERSION);
        return 0;
    }
    
    /* Dispatch to subcommand */
    if (strcmp(command, "map") == 0) {
        map_args_t args;
        if (parse_map_args(argc - 1, argv + 1, &args) != 0) {
            return 1;
        }
        return map_run(&args);
    }
    else if (strcmp(command, "sample") == 0) {
        sample_args_t args;
        if (parse_sample_args(argc - 1, argv + 1, &args) != 0) {
            return 1;
        }
        return sample_run(&args);
    }
    else {
        fprintf(stderr, "Error: Unknown command '%s'\n\n", command);
        print_usage();
        return 1;
    }
}

