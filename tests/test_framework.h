/*
 * test_framework.h - Minimal unit testing framework for samsampleX
 *
 * A header-only test framework that requires no external dependencies.
 * Provides basic assertions and test organization.
 */

#ifndef TEST_FRAMEWORK_H
#define TEST_FRAMEWORK_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* Colors for terminal output */
#define COLOR_RED     "\033[1;31m"
#define COLOR_GREEN   "\033[1;32m"
#define COLOR_YELLOW  "\033[1;33m"
#define COLOR_RESET   "\033[0m"

/* Global test counters */
static int _tests_run = 0;
static int _tests_passed = 0;
static int _tests_failed = 0;
static int _current_test_failed = 0;

/* Test function type */
typedef void (*test_fn)(void);

/* Begin a test case */
#define TEST_BEGIN(name) \
    static void test_##name(void) { \
        _tests_run++; \
        _current_test_failed = 0; \
        printf("  Running: %s ... ", #name);

/* End a test case */
#define TEST_END() \
        if (_current_test_failed) { \
            _tests_failed++; \
            printf(COLOR_RED "FAILED" COLOR_RESET "\n"); \
        } else { \
            _tests_passed++; \
            printf(COLOR_GREEN "OK" COLOR_RESET "\n"); \
        } \
    }

/* Run a test */
#define RUN_TEST(name) test_##name()

/* Assertion macros */
#define ASSERT_TRUE(cond) do { \
    if (!(cond)) { \
        fprintf(stderr, "\n    Assertion failed: %s (line %d)\n", #cond, __LINE__); \
        _current_test_failed = 1; \
        return; \
    } \
} while(0)

#define ASSERT_FALSE(cond) ASSERT_TRUE(!(cond))

#define ASSERT_EQ(expected, actual) do { \
    if ((expected) != (actual)) { \
        fprintf(stderr, "\n    Assertion failed: expected %ld, got %ld (line %d)\n", \
                (long)(expected), (long)(actual), __LINE__); \
        _current_test_failed = 1; \
        return; \
    } \
} while(0)

#define ASSERT_EQ_DBL(expected, actual, epsilon) do { \
    double _e = (expected), _a = (actual), _eps = (epsilon); \
    if (fabs(_e - _a) > _eps) { \
        fprintf(stderr, "\n    Assertion failed: expected %.6f, got %.6f (line %d)\n", \
                _e, _a, __LINE__); \
        _current_test_failed = 1; \
        return; \
    } \
} while(0)

#define ASSERT_STR_EQ(expected, actual) do { \
    const char *_e = (expected), *_a = (actual); \
    if (_e == NULL && _a == NULL) break; \
    if (_e == NULL || _a == NULL || strcmp(_e, _a) != 0) { \
        fprintf(stderr, "\n    Assertion failed: expected \"%s\", got \"%s\" (line %d)\n", \
                _e ? _e : "NULL", _a ? _a : "NULL", __LINE__); \
        _current_test_failed = 1; \
        return; \
    } \
} while(0)

#define ASSERT_NOT_NULL(ptr) do { \
    if ((ptr) == NULL) { \
        fprintf(stderr, "\n    Assertion failed: %s is NULL (line %d)\n", #ptr, __LINE__); \
        _current_test_failed = 1; \
        return; \
    } \
} while(0)

#define ASSERT_NULL(ptr) do { \
    if ((ptr) != NULL) { \
        fprintf(stderr, "\n    Assertion failed: %s is not NULL (line %d)\n", #ptr, __LINE__); \
        _current_test_failed = 1; \
        return; \
    } \
} while(0)

/* Test suite organization */
#define TEST_SUITE_BEGIN(name) \
    printf("\n" COLOR_YELLOW "Test Suite: %s" COLOR_RESET "\n", name);

/* Print summary and exit */
#define TEST_SUMMARY() do { \
    printf("\n========================================\n"); \
    printf("Tests run:    %d\n", _tests_run); \
    printf("Tests passed: " COLOR_GREEN "%d" COLOR_RESET "\n", _tests_passed); \
    printf("Tests failed: %s%d" COLOR_RESET "\n", \
           _tests_failed > 0 ? COLOR_RED : "", _tests_failed); \
    printf("========================================\n"); \
    return _tests_failed > 0 ? EXIT_FAILURE : EXIT_SUCCESS; \
} while(0)

#endif /* TEST_FRAMEWORK_H */

