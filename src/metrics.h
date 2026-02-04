/*
 * metrics.h - Metrics for comparing depth distributions
 */

#ifndef METRICS_H
#define METRICS_H

#include "samsampleX.h"

/*
 * Results from metric calculations
 */
typedef struct {
    double wasserstein;     /* Wasserstein-1 distance */
    double tv;              /* Total Variation distance (per-position) */
    double mean_template;   /* Mean depth in template */
    double mean_output;     /* Mean depth in output */
} metrics_result_t;

/*
 * Calculate metrics comparing template and output depth arrays.
 *
 * @template: Template depth array
 * @output: Output depth array
 * @result: Output structure for results
 * @return: 0 on success, non-zero on error
 */
int metrics_calculate(const depth_array_t *template_depth,
                      const depth_array_t *output_depth,
                      metrics_result_t *result);

/*
 * Print metrics to stderr.
 */
void metrics_print(const metrics_result_t *result);

#endif /* METRICS_H */

