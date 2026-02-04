/*
 * metrics.c - Metrics for comparing depth distributions
 *
 * Implements:
 * - Wasserstein-1 distance (Earth Mover's Distance) between depth distributions
 * - Total Variation (TV) distance for per-base depth comparison
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "metrics.h"
#include "samsampleX.h"

/*
 * Calculate Wasserstein-1 distance between two 1D discrete distributions.
 *
 * For 1D distributions, W1 = integral |CDF_A(x) - CDF_B(x)| dx
 * For discrete distributions: W1 = sum |cumsum(A) - cumsum(B)| / n
 *
 * We normalize by total depth to make it scale-invariant.
 */
static double calculate_wasserstein(const int32_t *depths_a, const int32_t *depths_b, 
                                    size_t n) {
    if (n == 0) return 0.0;
    
    /* Calculate cumulative sums */
    double *cumsum_a, *cumsum_b;
    SAFE_MALLOC(cumsum_a, (n + 1) * sizeof(double));
    SAFE_MALLOC(cumsum_b, (n + 1) * sizeof(double));
    
    cumsum_a[0] = 0.0;
    cumsum_b[0] = 0.0;
    
    double total_a = 0.0;
    double total_b = 0.0;
    
    for (size_t i = 0; i < n; i++) {
        cumsum_a[i + 1] = cumsum_a[i] + (double)depths_a[i];
        cumsum_b[i + 1] = cumsum_b[i] + (double)depths_b[i];
        total_a += (double)depths_a[i];
        total_b += (double)depths_b[i];
    }
    
    /* Normalize cumulative sums to [0, 1] (CDFs) */
    if (total_a > 0) {
        for (size_t i = 0; i <= n; i++) {
            cumsum_a[i] /= total_a;
        }
    }
    if (total_b > 0) {
        for (size_t i = 0; i <= n; i++) {
            cumsum_b[i] /= total_b;
        }
    }
    
    /* Calculate W1 = sum of |CDF_A - CDF_B| */
    double w1 = 0.0;
    for (size_t i = 0; i <= n; i++) {
        double diff = cumsum_a[i] - cumsum_b[i];
        if (diff < 0) diff = -diff;
        w1 += diff;
    }
    
    /* Normalize by number of positions */
    w1 /= (double)(n + 1);
    
    free(cumsum_a);
    free(cumsum_b);
    
    return w1;
}

/*
 * Calculate Total Variation distance between two depth arrays.
 * TV = (1/2) * sum|a[i] - b[i]| / n
 * 
 * This is half of MAE and represents the proportion of depth that differs
 * between the two arrays at each position on average.
 */
static double calculate_tv(const int32_t *depths_a, const int32_t *depths_b,
                           size_t n) {
    if (n == 0) return 0.0;
    
    double sum = 0.0;
    for (size_t i = 0; i < n; i++) {
        int32_t diff = depths_a[i] - depths_b[i];
        if (diff < 0) diff = -diff;
        sum += (double)diff;
    }
    
    return sum / (2.0 * (double)n);
}

/*
 * Calculate mean depth of an array.
 */
static double calculate_mean(const int32_t *depths, size_t n) {
    if (n == 0) return 0.0;
    
    double sum = 0.0;
    for (size_t i = 0; i < n; i++) {
        sum += (double)depths[i];
    }
    
    return sum / (double)n;
}

/*
 * Calculate metrics comparing template and output depth arrays.
 */
int metrics_calculate(const depth_array_t *template_depth,
                      const depth_array_t *output_depth,
                      metrics_result_t *result) {
    if (!template_depth || !output_depth || !result) {
        return -1;
    }
    
    /* Arrays must have same length */
    if (template_depth->length != output_depth->length) {
        fprintf(stderr, "Error: Depth arrays have different lengths (%zu vs %zu)\n",
                template_depth->length, output_depth->length);
        return -1;
    }
    
    size_t n = template_depth->length;
    
    /* Calculate metrics */
    result->wasserstein = calculate_wasserstein(template_depth->depths, 
                                                 output_depth->depths, n);
    result->tv = calculate_tv(template_depth->depths, output_depth->depths, n);
    result->mean_template = calculate_mean(template_depth->depths, n);
    result->mean_output = calculate_mean(output_depth->depths, n);
    
    return 0;
}

/*
 * Print metrics to stderr.
 */
void metrics_print(const metrics_result_t *result) {
    fprintf(stderr, "Mean Template Depth:     %.2f\n", result->mean_template);
    fprintf(stderr, "Mean Output Depth:       %.2f\n", result->mean_output);
    fprintf(stderr, "Total Variation:         %.4f\n", result->tv);
    fprintf(stderr, "Norm. Wasserstein Dist.:    %.6f\n", result->wasserstein);
    fprintf(stderr, "=======================================\n");
}

