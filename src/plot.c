/*
 * plot.c - Plot subcommand for visualizing depth of coverage
 *
 * Uses pbPlots (header-only C plotting library) to generate PNG line plots
 * comparing depth of coverage between source, template, and output BAM files.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "plot.h"
#include "depth.h"
#include "bed.h"
#include "samsampleX.h"

/* pbPlots library for PNG plot generation */
#include "pbPlots.h"
#include "supportLib.h"

/*
 * Print plot subcommand usage.
 */
void plot_usage(void) {
    fprintf(stderr, "Usage: %s plot [options]\n\n", SAMSAMPLEX_NAME);
    fprintf(stderr, "Compare depth of coverage and output as PNG plot or TSV data.\n\n");
    fprintf(stderr, "Required options:\n");
    fprintf(stderr, "  --source-bam FILE     Source BAM file\n");
    fprintf(stderr, "  --out-bam FILE        Output BAM file (from sampling)\n");
    fprintf(stderr, "  --region REGION       Target region (samtools-style, e.g. chr1:1000-2000)\n\n");
    fprintf(stderr, "Template options (one required):\n");
    fprintf(stderr, "  --template-bam FILE   Template BAM file\n");
    fprintf(stderr, "  --template-bed FILE   Template BED file\n\n");
    fprintf(stderr, "Output options (one required, mutually exclusive):\n");
    fprintf(stderr, "  --out-png FILE        Output PNG plot file\n");
    fprintf(stderr, "  --out-tsv FILE        Output TSV data file\n\n");
    fprintf(stderr, "Examples:\n");
    fprintf(stderr, "  # Generate PNG plot\n");
    fprintf(stderr, "  %s plot --source-bam source.bam --template-bam template.bam \\\n", SAMSAMPLEX_NAME);
    fprintf(stderr, "              --out-bam sampled.bam --region chr1:1000-2000 --out-png plot.png\n\n");
    fprintf(stderr, "  # Export TSV for custom plotting\n");
    fprintf(stderr, "  %s plot --source-bam source.bam --template-bed template.bed \\\n", SAMSAMPLEX_NAME);
    fprintf(stderr, "              --out-bam sampled.bam --region chr1:1000-2000 --out-tsv depths.tsv\n");
}

/*
 * Create an array of x-coordinates (genomic positions) for plotting.
 * Uses 1-based coordinates for display (matching user input).
 */
static double *create_x_coords(int64_t start_0based, size_t length) {
    double *xs = malloc(length * sizeof(double));
    if (!xs) {
        fprintf(stderr, "Error: Memory allocation failed for x coordinates\n");
        return NULL;
    }
    /* Use relative positions (0 to length-1) for better numeric stability */
    for (size_t i = 0; i < length; i++) {
        xs[i] = (double)i;
    }
    return xs;
}

/*
 * Convert depth array (int32_t) to double array for plotting.
 */
static double *depths_to_double(int32_t *depths, size_t length) {
    double *ys = malloc(length * sizeof(double));
    if (!ys) {
        fprintf(stderr, "Error: Memory allocation failed for y coordinates\n");
        return NULL;
    }
    for (size_t i = 0; i < length; i++) {
        ys[i] = (double)depths[i];
    }
    return ys;
}

/*
 * Write depth data to TSV file.
 * Format: position<TAB>source<TAB>template<TAB>output
 */
static int write_tsv(const char *filename, region_t *region,
                     depth_array_t *source, depth_array_t *template_d, 
                     depth_array_t *out) {
    FILE *fp = fopen(filename, "w");
    if (!fp) {
        fprintf(stderr, "Error: Cannot open TSV file for writing: %s\n", filename);
        return 1;
    }
    
    /* Write header */
    fprintf(fp, "position\tsource_depth\ttemplate_depth\toutput_depth\n");
    
    /* Write data rows with 1-based genomic positions */
    for (size_t i = 0; i < source->length; i++) {
        int64_t pos = region->start + 1 + (int64_t)i;  /* 1-based */
        fprintf(fp, "%ld\t%d\t%d\t%d\n", 
                pos, source->depths[i], template_d->depths[i], out->depths[i]);
    }
    
    fclose(fp);
    return 0;
}

/*
 * Run the plot subcommand.
 */
int plot_run(plot_args_t *args) {
    int ret = 0;
    
    /* Validate arguments */
    if (!args->source_bam) {
        fprintf(stderr, "Error: --source-bam is required\n");
        plot_usage();
        return 1;
    }
    if (!args->out_bam) {
        fprintf(stderr, "Error: --out-bam is required\n");
        plot_usage();
        return 1;
    }
    if (!args->region) {
        fprintf(stderr, "Error: --region is required\n");
        plot_usage();
        return 1;
    }
    if (!args->template_bam && !args->template_bed) {
        fprintf(stderr, "Error: Either --template-bam or --template-bed is required\n");
        plot_usage();
        return 1;
    }
    if (args->template_bam && args->template_bed) {
        fprintf(stderr, "Error: --template-bam and --template-bed are mutually exclusive\n");
        plot_usage();
        return 1;
    }
    if (!args->out_png && !args->out_tsv) {
        fprintf(stderr, "Error: Either --out-png or --out-tsv is required\n");
        plot_usage();
        return 1;
    }
    if (args->out_png && args->out_tsv) {
        fprintf(stderr, "Error: --out-png and --out-tsv are mutually exclusive\n");
        plot_usage();
        return 1;
    }
    
    /* Parse region */
    region_t *region = region_parse(args->region);
    if (!region) {
        fprintf(stderr, "Error: Failed to parse region: %s\n", args->region);
        return 1;
    }
    
    fprintf(stderr, "[plot] Region: %s:%ld-%ld\n", region->contig, region->start + 1, region->end);
    
    /* Load depth arrays */
    fprintf(stderr, "[plot] Loading source depths from: %s\n", args->source_bam);
    depth_array_t *source_depth = depth_from_bam(args->source_bam, region->contig, 
                                                  region->start, region->end);
    if (!source_depth) {
        fprintf(stderr, "Error: Failed to load source BAM depths\n");
        ret = 1;
        goto cleanup_region;
    }
    
    /* Update region end if it was set from BAM */
    if (region->end < 0) {
        region->end = source_depth->end;
    }
    
    /* Load template depths (from BAM or BED) */
    depth_array_t *template_depth = NULL;
    if (args->template_bam) {
        fprintf(stderr, "[plot] Loading template depths from BAM: %s\n", args->template_bam);
        template_depth = depth_from_bam(args->template_bam, region->contig,
                                        region->start, region->end);
    } else {
        fprintf(stderr, "[plot] Loading template depths from BED: %s\n", args->template_bed);
        template_depth = bed_read_depths(args->template_bed, region->contig,
                                         region->start, region->end);
    }
    
    if (!template_depth) {
        fprintf(stderr, "Error: Failed to load template depths\n");
        ret = 1;
        goto cleanup_source;
    }
    
    fprintf(stderr, "[plot] Loading output depths from: %s\n", args->out_bam);
    depth_array_t *out_depth = depth_from_bam(args->out_bam, region->contig,
                                              region->start, region->end);
    if (!out_depth) {
        fprintf(stderr, "Error: Failed to load output BAM depths\n");
        ret = 1;
        goto cleanup_template;
    }
    
    /* Verify all arrays have same length */
    if (source_depth->length != template_depth->length || 
        source_depth->length != out_depth->length) {
        fprintf(stderr, "Error: Depth array length mismatch\n");
        ret = 1;
        goto cleanup_out;
    }
    
    size_t n_points = source_depth->length;
    fprintf(stderr, "[plot] Processing %zu data points\n", n_points);
    
    /* Branch based on output format */
    if (args->out_tsv) {
        /* TSV output mode */
        fprintf(stderr, "[plot] Writing TSV to: %s\n", args->out_tsv);
        ret = write_tsv(args->out_tsv, region, source_depth, template_depth, out_depth);
        if (ret == 0) {
            fprintf(stderr, "[plot] Done.\n");
            fprintf(stderr, "[plot] Columns: position, source_depth, template_depth, output_depth\n");
        }
        goto cleanup_out;
    }
    
    /* PNG output mode */
    /* Initialize pbPlots memory arena */
    StartArenaAllocator();
    
    /* Create x-coordinates (genomic positions) */
    double *xs = create_x_coords(region->start, n_points);
    if (!xs) {
        ret = 1;
        goto cleanup_arena;
    }
    
    /* Convert depths to double arrays */
    double *source_ys = depths_to_double(source_depth->depths, n_points);
    double *template_ys = depths_to_double(template_depth->depths, n_points);
    double *out_ys = depths_to_double(out_depth->depths, n_points);
    
    if (!source_ys || !template_ys || !out_ys) {
        ret = 1;
        goto cleanup_arrays;
    }
    
    /* Set up scatter plot settings */
    ScatterPlotSettings *settings = GetDefaultScatterPlotSettings();
    settings->width = DEFAULT_PLOT_WIDTH;
    settings->height = DEFAULT_PLOT_HEIGHT;
    settings->autoBoundaries = true;
    settings->autoPadding = true;
    settings->title = L"Depth of Coverage Comparison";
    settings->titleLength = wcslen(settings->title);
    settings->xLabel = L"Relative Position (bp)";
    settings->xLabelLength = wcslen(settings->xLabel);
    settings->yLabel = L"Depth";
    settings->yLabelLength = wcslen(settings->yLabel);
    settings->showGrid = true;
    
    /* Create three series for source, template, and output */
    ScatterPlotSeries *series[3];
    settings->scatterPlotSeries = series;
    settings->scatterPlotSeriesLength = 3;
    
    /* Source series (blue) */
    series[0] = GetDefaultScatterPlotSeriesSettings();
    series[0]->xs = xs;
    series[0]->xsLength = n_points;
    series[0]->ys = source_ys;
    series[0]->ysLength = n_points;
    series[0]->linearInterpolation = true;
    series[0]->lineType = L"solid";
    series[0]->lineTypeLength = wcslen(series[0]->lineType);
    series[0]->lineThickness = 1.0;
    series[0]->pointType = L"";
    series[0]->pointTypeLength = 0;
    series[0]->color = CreateRGBColor(0.2, 0.4, 0.8);  /* Blue */
    
    /* Template series (green) */
    series[1] = GetDefaultScatterPlotSeriesSettings();
    series[1]->xs = xs;
    series[1]->xsLength = n_points;
    series[1]->ys = template_ys;
    series[1]->ysLength = n_points;
    series[1]->linearInterpolation = true;
    series[1]->lineType = L"solid";
    series[1]->lineTypeLength = wcslen(series[1]->lineType);
    series[1]->lineThickness = 1.0;
    series[1]->pointType = L"";
    series[1]->pointTypeLength = 0;
    series[1]->color = CreateRGBColor(0.2, 0.7, 0.3);  /* Green */
    
    /* Output series (red) */
    series[2] = GetDefaultScatterPlotSeriesSettings();
    series[2]->xs = xs;
    series[2]->xsLength = n_points;
    series[2]->ys = out_ys;
    series[2]->ysLength = n_points;
    series[2]->linearInterpolation = true;
    series[2]->lineType = L"solid";
    series[2]->lineTypeLength = wcslen(series[2]->lineType);
    series[2]->lineThickness = 1.0;
    series[2]->pointType = L"";
    series[2]->pointTypeLength = 0;
    series[2]->color = CreateRGBColor(0.8, 0.2, 0.2);  /* Red */
    
    /* Create canvas and draw plot */
    RGBABitmapImageReference *canvasRef = CreateRGBABitmapImageReference();
    StringReference *errorMessage = CreateStringReference(L"", 0);
    
    _Bool success = DrawScatterPlotFromSettings(canvasRef, settings, errorMessage);
    
    if (!success) {
        fprintf(stderr, "Error: Failed to draw plot: %ls\n", errorMessage->string);
        ret = 1;
        goto cleanup_arrays;
    }
    
    /* Convert to PNG and write to file */
    ByteArray *pngData = ConvertToPNG(canvasRef->image);
    
    fprintf(stderr, "[plot] Writing PNG to: %s\n", args->out_png);
    WriteToFile(pngData, (char *)args->out_png);
    
    fprintf(stderr, "[plot] Done.\n");
    fprintf(stderr, "[plot] Legend: Blue=Source, Green=Template, Red=Output\n");
    
cleanup_arrays:
    free(xs);
    free(source_ys);
    free(template_ys);
    free(out_ys);
cleanup_arena:
    FreeAllocations();
cleanup_out:
    depth_array_free(out_depth);
cleanup_template:
    depth_array_free(template_depth);
cleanup_source:
    depth_array_free(source_depth);
cleanup_region:
    region_free(region);
    
    return ret;
}

