/*
 * test_metrics.c - Unit tests for metrics calculations (Wasserstein, MAE)
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "test_framework.h"
#include "../include/samsampleX.h"
#include "../src/bed.h"
#include "../src/metrics.h"

/* ============================================================
 * MAE Tests
 * ============================================================ */

TEST_BEGIN(mae_identical_arrays)
    depth_array_t *a = depth_array_alloc("chr1", 0, 5);
    depth_array_t *b = depth_array_alloc("chr1", 0, 5);
    
    for (int i = 0; i < 5; i++) {
        a->depths[i] = 100;
        b->depths[i] = 100;
    }
    
    metrics_result_t result;
    int ret = metrics_calculate(a, b, &result);
    
    ASSERT_EQ(0, ret);
    ASSERT_EQ_DBL(0.0, result.mae, 0.0001);
    
    depth_array_free(a);
    depth_array_free(b);
TEST_END()

TEST_BEGIN(mae_uniform_difference)
    depth_array_t *a = depth_array_alloc("chr1", 0, 4);
    depth_array_t *b = depth_array_alloc("chr1", 0, 4);
    
    /* All positions differ by 10 */
    a->depths[0] = 100; b->depths[0] = 110;
    a->depths[1] = 50;  b->depths[1] = 60;
    a->depths[2] = 75;  b->depths[2] = 85;
    a->depths[3] = 200; b->depths[3] = 210;
    
    metrics_result_t result;
    int ret = metrics_calculate(a, b, &result);
    
    ASSERT_EQ(0, ret);
    ASSERT_EQ_DBL(10.0, result.mae, 0.0001);
    
    depth_array_free(a);
    depth_array_free(b);
TEST_END()

TEST_BEGIN(mae_mixed_difference)
    depth_array_t *a = depth_array_alloc("chr1", 0, 4);
    depth_array_t *b = depth_array_alloc("chr1", 0, 4);
    
    /* Differences: 10, 20, 0, 30 => MAE = 60/4 = 15 */
    a->depths[0] = 100; b->depths[0] = 110;
    a->depths[1] = 50;  b->depths[1] = 70;
    a->depths[2] = 75;  b->depths[2] = 75;
    a->depths[3] = 200; b->depths[3] = 230;
    
    metrics_result_t result;
    int ret = metrics_calculate(a, b, &result);
    
    ASSERT_EQ(0, ret);
    ASSERT_EQ_DBL(15.0, result.mae, 0.0001);
    
    depth_array_free(a);
    depth_array_free(b);
TEST_END()

TEST_BEGIN(mae_negative_differences)
    depth_array_t *a = depth_array_alloc("chr1", 0, 4);
    depth_array_t *b = depth_array_alloc("chr1", 0, 4);
    
    /* Some b values less than a */
    a->depths[0] = 100; b->depths[0] = 90;   /* |100-90| = 10 */
    a->depths[1] = 50;  b->depths[1] = 70;   /* |50-70| = 20 */
    a->depths[2] = 75;  b->depths[2] = 55;   /* |75-55| = 20 */
    a->depths[3] = 200; b->depths[3] = 210;  /* |200-210| = 10 */
    /* MAE = 60/4 = 15 */
    
    metrics_result_t result;
    int ret = metrics_calculate(a, b, &result);
    
    ASSERT_EQ(0, ret);
    ASSERT_EQ_DBL(15.0, result.mae, 0.0001);
    
    depth_array_free(a);
    depth_array_free(b);
TEST_END()

/* ============================================================
 * Wasserstein Tests
 * ============================================================ */

TEST_BEGIN(wasserstein_identical)
    depth_array_t *a = depth_array_alloc("chr1", 0, 5);
    depth_array_t *b = depth_array_alloc("chr1", 0, 5);
    
    for (int i = 0; i < 5; i++) {
        a->depths[i] = 100;
        b->depths[i] = 100;
    }
    
    metrics_result_t result;
    int ret = metrics_calculate(a, b, &result);
    
    ASSERT_EQ(0, ret);
    ASSERT_EQ_DBL(0.0, result.wasserstein, 0.0001);
    
    depth_array_free(a);
    depth_array_free(b);
TEST_END()

TEST_BEGIN(wasserstein_scaled)
    /* If one array is exactly half the other, normalized CDFs should match */
    depth_array_t *a = depth_array_alloc("chr1", 0, 4);
    depth_array_t *b = depth_array_alloc("chr1", 0, 4);
    
    a->depths[0] = 100; b->depths[0] = 50;
    a->depths[1] = 200; b->depths[1] = 100;
    a->depths[2] = 100; b->depths[2] = 50;
    a->depths[3] = 200; b->depths[3] = 100;
    
    metrics_result_t result;
    int ret = metrics_calculate(a, b, &result);
    
    ASSERT_EQ(0, ret);
    /* Proportionally identical distributions should have W1 = 0 */
    ASSERT_EQ_DBL(0.0, result.wasserstein, 0.0001);
    
    depth_array_free(a);
    depth_array_free(b);
TEST_END()

TEST_BEGIN(wasserstein_different)
    /* Different distributions should have non-zero W1 */
    depth_array_t *a = depth_array_alloc("chr1", 0, 4);
    depth_array_t *b = depth_array_alloc("chr1", 0, 4);
    
    /* Front-loaded vs back-loaded */
    a->depths[0] = 100; b->depths[0] = 0;
    a->depths[1] = 100; b->depths[1] = 0;
    a->depths[2] = 0;   b->depths[2] = 100;
    a->depths[3] = 0;   b->depths[3] = 100;
    
    metrics_result_t result;
    int ret = metrics_calculate(a, b, &result);
    
    ASSERT_EQ(0, ret);
    /* Should be non-zero for different distributions */
    ASSERT_TRUE(result.wasserstein > 0.0);
    
    depth_array_free(a);
    depth_array_free(b);
TEST_END()

/* ============================================================
 * Mean Depth Tests
 * ============================================================ */

TEST_BEGIN(mean_depth_uniform)
    depth_array_t *a = depth_array_alloc("chr1", 0, 5);
    
    for (int i = 0; i < 5; i++) {
        a->depths[i] = 100;
    }
    
    depth_array_t *b = depth_array_alloc("chr1", 0, 5);
    for (int i = 0; i < 5; i++) {
        b->depths[i] = 50;
    }
    
    metrics_result_t result;
    int ret = metrics_calculate(a, b, &result);
    
    ASSERT_EQ(0, ret);
    ASSERT_EQ_DBL(100.0, result.mean_template, 0.0001);
    ASSERT_EQ_DBL(50.0, result.mean_output, 0.0001);
    
    depth_array_free(a);
    depth_array_free(b);
TEST_END()

TEST_BEGIN(mean_depth_varied)
    depth_array_t *a = depth_array_alloc("chr1", 0, 4);
    depth_array_t *b = depth_array_alloc("chr1", 0, 4);
    
    a->depths[0] = 10;  b->depths[0] = 5;
    a->depths[1] = 20;  b->depths[1] = 15;
    a->depths[2] = 30;  b->depths[2] = 25;
    a->depths[3] = 40;  b->depths[3] = 35;
    /* Mean a = 25, mean b = 20 */
    
    metrics_result_t result;
    int ret = metrics_calculate(a, b, &result);
    
    ASSERT_EQ(0, ret);
    ASSERT_EQ_DBL(25.0, result.mean_template, 0.0001);
    ASSERT_EQ_DBL(20.0, result.mean_output, 0.0001);
    
    depth_array_free(a);
    depth_array_free(b);
TEST_END()

/* ============================================================
 * Error Cases
 * ============================================================ */

TEST_BEGIN(metrics_null_inputs)
    metrics_result_t result;
    
    /* NULL template */
    ASSERT_EQ(-1, metrics_calculate(NULL, NULL, &result));
    
    depth_array_t *a = depth_array_alloc("chr1", 0, 5);
    ASSERT_EQ(-1, metrics_calculate(a, NULL, &result));
    ASSERT_EQ(-1, metrics_calculate(NULL, a, &result));
    ASSERT_EQ(-1, metrics_calculate(a, a, NULL));
    
    depth_array_free(a);
TEST_END()

TEST_BEGIN(metrics_length_mismatch)
    depth_array_t *a = depth_array_alloc("chr1", 0, 5);
    depth_array_t *b = depth_array_alloc("chr1", 0, 10);
    
    metrics_result_t result;
    /* Should return error for mismatched lengths */
    ASSERT_EQ(-1, metrics_calculate(a, b, &result));
    
    depth_array_free(a);
    depth_array_free(b);
TEST_END()

/* ============================================================
 * Main
 * ============================================================ */

int main(void) {
    TEST_SUITE_BEGIN("Mean Absolute Error");
    RUN_TEST(mae_identical_arrays);
    RUN_TEST(mae_uniform_difference);
    RUN_TEST(mae_mixed_difference);
    RUN_TEST(mae_negative_differences);
    
    TEST_SUITE_BEGIN("Wasserstein Distance");
    RUN_TEST(wasserstein_identical);
    RUN_TEST(wasserstein_scaled);
    RUN_TEST(wasserstein_different);
    
    TEST_SUITE_BEGIN("Mean Depth Calculation");
    RUN_TEST(mean_depth_uniform);
    RUN_TEST(mean_depth_varied);
    
    TEST_SUITE_BEGIN("Error Cases");
    RUN_TEST(metrics_null_inputs);
    RUN_TEST(metrics_length_mismatch);
    
    TEST_SUMMARY();
}

