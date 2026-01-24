/*
 * test_depth.c - Unit tests for depth array operations and BED file utilities
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "test_framework.h"
#include "../include/samsampleX.h"
#include "../src/bed.h"

/* ============================================================
 * Depth Array Allocation Tests
 * ============================================================ */

TEST_BEGIN(depth_array_alloc_basic)
    depth_array_t *arr = depth_array_alloc("chr1", 0, 100);
    ASSERT_NOT_NULL(arr);
    ASSERT_NOT_NULL(arr->depths);
    ASSERT_STR_EQ("chr1", arr->contig);
    ASSERT_EQ(0, arr->start);
    ASSERT_EQ(100, arr->end);
    ASSERT_EQ(100, arr->length);
    
    /* Check that depths are initialized to 0 (calloc) */
    for (size_t i = 0; i < arr->length; i++) {
        ASSERT_EQ(0, arr->depths[i]);
    }
    
    depth_array_free(arr);
TEST_END()

TEST_BEGIN(depth_array_alloc_offset_region)
    depth_array_t *arr = depth_array_alloc("chr21", 1000, 2000);
    ASSERT_NOT_NULL(arr);
    ASSERT_STR_EQ("chr21", arr->contig);
    ASSERT_EQ(1000, arr->start);
    ASSERT_EQ(2000, arr->end);
    ASSERT_EQ(1000, arr->length);
    depth_array_free(arr);
TEST_END()

TEST_BEGIN(depth_array_free_null)
    /* Should not crash on NULL */
    depth_array_free(NULL);
    ASSERT_TRUE(1);  /* If we got here, test passed */
TEST_END()

/* ============================================================
 * Depth Array Combine Tests
 * ============================================================ */

TEST_BEGIN(combine_depths_single_array)
    depth_array_t *arr = depth_array_alloc("chr1", 0, 5);
    arr->depths[0] = 10;
    arr->depths[1] = 20;
    arr->depths[2] = 30;
    arr->depths[3] = 40;
    arr->depths[4] = 50;
    
    depth_array_t *arrays[1] = { arr };
    depth_array_t *result = bed_combine_depths(arrays, 1, MODE_MIN, 42);
    
    ASSERT_NOT_NULL(result);
    ASSERT_EQ(5, result->length);
    ASSERT_EQ(10, result->depths[0]);
    ASSERT_EQ(20, result->depths[1]);
    ASSERT_EQ(30, result->depths[2]);
    ASSERT_EQ(40, result->depths[3]);
    ASSERT_EQ(50, result->depths[4]);
    
    depth_array_free(arr);
    depth_array_free(result);
TEST_END()

TEST_BEGIN(combine_depths_mode_min)
    depth_array_t *a = depth_array_alloc("chr1", 0, 4);
    depth_array_t *b = depth_array_alloc("chr1", 0, 4);
    
    a->depths[0] = 10; b->depths[0] = 20;
    a->depths[1] = 30; b->depths[1] = 15;
    a->depths[2] = 25; b->depths[2] = 25;
    a->depths[3] = 40; b->depths[3] = 5;
    
    depth_array_t *arrays[2] = { a, b };
    depth_array_t *result = bed_combine_depths(arrays, 2, MODE_MIN, 42);
    
    ASSERT_NOT_NULL(result);
    ASSERT_EQ(10, result->depths[0]);  /* min(10, 20) */
    ASSERT_EQ(15, result->depths[1]);  /* min(30, 15) */
    ASSERT_EQ(25, result->depths[2]);  /* min(25, 25) */
    ASSERT_EQ(5, result->depths[3]);   /* min(40, 5) */
    
    depth_array_free(a);
    depth_array_free(b);
    depth_array_free(result);
TEST_END()

TEST_BEGIN(combine_depths_mode_max)
    depth_array_t *a = depth_array_alloc("chr1", 0, 4);
    depth_array_t *b = depth_array_alloc("chr1", 0, 4);
    
    a->depths[0] = 10; b->depths[0] = 20;
    a->depths[1] = 30; b->depths[1] = 15;
    a->depths[2] = 25; b->depths[2] = 25;
    a->depths[3] = 40; b->depths[3] = 5;
    
    depth_array_t *arrays[2] = { a, b };
    depth_array_t *result = bed_combine_depths(arrays, 2, MODE_MAX, 42);
    
    ASSERT_NOT_NULL(result);
    ASSERT_EQ(20, result->depths[0]);  /* max(10, 20) */
    ASSERT_EQ(30, result->depths[1]);  /* max(30, 15) */
    ASSERT_EQ(25, result->depths[2]);  /* max(25, 25) */
    ASSERT_EQ(40, result->depths[3]);  /* max(40, 5) */
    
    depth_array_free(a);
    depth_array_free(b);
    depth_array_free(result);
TEST_END()

TEST_BEGIN(combine_depths_mode_mean)
    depth_array_t *a = depth_array_alloc("chr1", 0, 4);
    depth_array_t *b = depth_array_alloc("chr1", 0, 4);
    
    a->depths[0] = 10; b->depths[0] = 20;  /* mean = 15 */
    a->depths[1] = 30; b->depths[1] = 10;  /* mean = 20 */
    a->depths[2] = 25; b->depths[2] = 25;  /* mean = 25 */
    a->depths[3] = 41; b->depths[3] = 9;   /* mean = 25 (integer division) */
    
    depth_array_t *arrays[2] = { a, b };
    depth_array_t *result = bed_combine_depths(arrays, 2, MODE_MEAN, 42);
    
    ASSERT_NOT_NULL(result);
    ASSERT_EQ(15, result->depths[0]);
    ASSERT_EQ(20, result->depths[1]);
    ASSERT_EQ(25, result->depths[2]);
    ASSERT_EQ(25, result->depths[3]);
    
    depth_array_free(a);
    depth_array_free(b);
    depth_array_free(result);
TEST_END()

TEST_BEGIN(combine_depths_three_arrays)
    depth_array_t *a = depth_array_alloc("chr1", 0, 3);
    depth_array_t *b = depth_array_alloc("chr1", 0, 3);
    depth_array_t *c = depth_array_alloc("chr1", 0, 3);
    
    a->depths[0] = 10; b->depths[0] = 20; c->depths[0] = 30;
    a->depths[1] = 50; b->depths[1] = 25; c->depths[1] = 10;
    a->depths[2] = 15; b->depths[2] = 15; c->depths[2] = 15;
    
    depth_array_t *arrays[3] = { a, b, c };
    
    /* Test min */
    depth_array_t *result_min = bed_combine_depths(arrays, 3, MODE_MIN, 42);
    ASSERT_NOT_NULL(result_min);
    ASSERT_EQ(10, result_min->depths[0]);
    ASSERT_EQ(10, result_min->depths[1]);
    ASSERT_EQ(15, result_min->depths[2]);
    
    /* Test max */
    depth_array_t *result_max = bed_combine_depths(arrays, 3, MODE_MAX, 42);
    ASSERT_NOT_NULL(result_max);
    ASSERT_EQ(30, result_max->depths[0]);
    ASSERT_EQ(50, result_max->depths[1]);
    ASSERT_EQ(15, result_max->depths[2]);
    
    /* Test mean: (10+20+30)/3=20, (50+25+10)/3=28, (15+15+15)/3=15 */
    depth_array_t *result_mean = bed_combine_depths(arrays, 3, MODE_MEAN, 42);
    ASSERT_NOT_NULL(result_mean);
    ASSERT_EQ(20, result_mean->depths[0]);
    ASSERT_EQ(28, result_mean->depths[1]);
    ASSERT_EQ(15, result_mean->depths[2]);
    
    depth_array_free(a);
    depth_array_free(b);
    depth_array_free(c);
    depth_array_free(result_min);
    depth_array_free(result_max);
    depth_array_free(result_mean);
TEST_END()

TEST_BEGIN(combine_depths_null_input)
    depth_array_t *result = bed_combine_depths(NULL, 0, MODE_MIN, 42);
    ASSERT_NULL(result);
TEST_END()

/* ============================================================
 * Main
 * ============================================================ */

int main(void) {
    TEST_SUITE_BEGIN("Depth Array Allocation");
    RUN_TEST(depth_array_alloc_basic);
    RUN_TEST(depth_array_alloc_offset_region);
    RUN_TEST(depth_array_free_null);
    
    TEST_SUITE_BEGIN("Depth Array Combine");
    RUN_TEST(combine_depths_single_array);
    RUN_TEST(combine_depths_mode_min);
    RUN_TEST(combine_depths_mode_max);
    RUN_TEST(combine_depths_mode_mean);
    RUN_TEST(combine_depths_three_arrays);
    RUN_TEST(combine_depths_null_input);
    
    TEST_SUMMARY();
}

