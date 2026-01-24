/*
 * test_parsing.c - Unit tests for region parsing and combine mode parsing
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "test_framework.h"
#include "../include/samsampleX.h"
#include "../src/depth.h"

/* ============================================================
 * Region Parsing Tests
 * ============================================================ */

TEST_BEGIN(region_parse_contig_only)
    region_t *r = region_parse("chr1");
    ASSERT_NOT_NULL(r);
    ASSERT_STR_EQ("chr1", r->contig);
    ASSERT_EQ(-1, r->start);
    ASSERT_EQ(-1, r->end);
    region_free(r);
TEST_END()

TEST_BEGIN(region_parse_contig_with_start)
    region_t *r = region_parse("chr1:1000");
    ASSERT_NOT_NULL(r);
    ASSERT_STR_EQ("chr1", r->contig);
    ASSERT_EQ(999, r->start);  /* 0-based, so 1000 -> 999 */
    ASSERT_EQ(-1, r->end);
    region_free(r);
TEST_END()

TEST_BEGIN(region_parse_contig_with_range)
    region_t *r = region_parse("chr21:1000-2000");
    ASSERT_NOT_NULL(r);
    ASSERT_STR_EQ("chr21", r->contig);
    ASSERT_EQ(999, r->start);   /* 0-based: 1000 -> 999 */
    ASSERT_EQ(2000, r->end);    /* End is exclusive */
    region_free(r);
TEST_END()

TEST_BEGIN(region_parse_numeric_contig)
    region_t *r = region_parse("22:5000000-10000000");
    ASSERT_NOT_NULL(r);
    ASSERT_STR_EQ("22", r->contig);
    ASSERT_EQ(4999999, r->start);
    ASSERT_EQ(10000000, r->end);
    region_free(r);
TEST_END()

TEST_BEGIN(region_parse_start_at_one)
    region_t *r = region_parse("chrX:1-100");
    ASSERT_NOT_NULL(r);
    ASSERT_STR_EQ("chrX", r->contig);
    ASSERT_EQ(0, r->start);     /* 1-based 1 -> 0-based 0 */
    ASSERT_EQ(100, r->end);
    region_free(r);
TEST_END()

/* ============================================================
 * Combine Mode Parsing Tests
 * ============================================================ */

TEST_BEGIN(parse_combine_mode_min)
    ASSERT_EQ(MODE_MIN, parse_combine_mode("min"));
TEST_END()

TEST_BEGIN(parse_combine_mode_max)
    ASSERT_EQ(MODE_MAX, parse_combine_mode("max"));
TEST_END()

TEST_BEGIN(parse_combine_mode_mean)
    ASSERT_EQ(MODE_MEAN, parse_combine_mode("mean"));
TEST_END()

TEST_BEGIN(parse_combine_mode_random)
    ASSERT_EQ(MODE_RANDOM, parse_combine_mode("random"));
TEST_END()

TEST_BEGIN(parse_combine_mode_invalid)
    /* Invalid mode should default to MODE_MIN with a warning */
    ASSERT_EQ(MODE_MIN, parse_combine_mode("invalid"));
TEST_END()

/* ============================================================
 * Mode to String Tests
 * ============================================================ */

TEST_BEGIN(combine_mode_to_string_min)
    ASSERT_STR_EQ("min", combine_mode_to_string(MODE_MIN));
TEST_END()

TEST_BEGIN(combine_mode_to_string_max)
    ASSERT_STR_EQ("max", combine_mode_to_string(MODE_MAX));
TEST_END()

TEST_BEGIN(combine_mode_to_string_mean)
    ASSERT_STR_EQ("mean", combine_mode_to_string(MODE_MEAN));
TEST_END()

TEST_BEGIN(combine_mode_to_string_random)
    ASSERT_STR_EQ("random", combine_mode_to_string(MODE_RANDOM));
TEST_END()

/* ============================================================
 * Main
 * ============================================================ */

int main(void) {
    TEST_SUITE_BEGIN("Region Parsing");
    RUN_TEST(region_parse_contig_only);
    RUN_TEST(region_parse_contig_with_start);
    RUN_TEST(region_parse_contig_with_range);
    RUN_TEST(region_parse_numeric_contig);
    RUN_TEST(region_parse_start_at_one);
    
    TEST_SUITE_BEGIN("Combine Mode Parsing");
    RUN_TEST(parse_combine_mode_min);
    RUN_TEST(parse_combine_mode_max);
    RUN_TEST(parse_combine_mode_mean);
    RUN_TEST(parse_combine_mode_random);
    RUN_TEST(parse_combine_mode_invalid);
    
    TEST_SUITE_BEGIN("Mode to String");
    RUN_TEST(combine_mode_to_string_min);
    RUN_TEST(combine_mode_to_string_max);
    RUN_TEST(combine_mode_to_string_mean);
    RUN_TEST(combine_mode_to_string_random);
    
    TEST_SUMMARY();
}

