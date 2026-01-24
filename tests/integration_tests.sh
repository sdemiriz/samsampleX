#!/bin/bash
#
# integration_tests.sh - Integration tests for samsampleX CLI
#
# These tests verify that the CLI interface works correctly.
# They use temporary files and don't require real BAM data for basic checks.

set -e

# Colors
RED='\033[1;31m'
GREEN='\033[1;32m'
YELLOW='\033[1;33m'
NC='\033[0m'

# Test counters
TESTS_RUN=0
TESTS_PASSED=0
TESTS_FAILED=0

# Get script directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"
SAMSAMPLEX="$PROJECT_DIR/samsampleX"

# Temp directory
TMPDIR=$(mktemp -d)
trap "rm -rf $TMPDIR" EXIT

# Test helper functions
pass() {
    echo -e "  ${GREEN}PASS${NC}: $1"
    TESTS_PASSED=$((TESTS_PASSED + 1))
    TESTS_RUN=$((TESTS_RUN + 1))
}

fail() {
    echo -e "  ${RED}FAIL${NC}: $1"
    if [ -n "$2" ]; then
        echo -e "        $2"
    fi
    TESTS_FAILED=$((TESTS_FAILED + 1))
    TESTS_RUN=$((TESTS_RUN + 1))
}

test_suite() {
    echo -e "\n${YELLOW}Test Suite: $1${NC}"
}

# ============================================================
# Basic CLI Tests
# ============================================================

test_suite "Basic CLI"

# Test: Version flag
if "$SAMSAMPLEX" --version 2>&1 | grep -q "samsampleX"; then
    pass "--version shows program name"
else
    fail "--version should show program name"
fi

# Test: Help flag
if "$SAMSAMPLEX" --help 2>&1 | grep -q "Usage"; then
    pass "--help shows usage"
else
    fail "--help should show usage"
fi

# Test: Help flag short form
if "$SAMSAMPLEX" -h 2>&1 | grep -q "Usage"; then
    pass "-h shows usage"
else
    fail "-h should show usage"
fi

# Test: No arguments shows usage
if "$SAMSAMPLEX" 2>&1 | grep -q "Usage"; then
    pass "No arguments shows usage"
else
    fail "No arguments should show usage"
fi

# ============================================================
# Map Subcommand Tests
# ============================================================

test_suite "Map Subcommand"

# Test: Map help
if "$SAMSAMPLEX" map --help 2>&1 | grep -q -i "map"; then
    pass "map --help shows map usage"
else
    fail "map --help should show map usage"
fi

# Test: Map without required args
if ! "$SAMSAMPLEX" map 2>/dev/null; then
    pass "map without args returns error"
else
    fail "map without args should return error"
fi

# Test: Map with missing input file
if ! "$SAMSAMPLEX" map --input nonexistent.bam --region chr1 --output /dev/null 2>/dev/null; then
    pass "map with missing input returns error"
else
    fail "map with missing input should return error"
fi

# ============================================================
# Sample Subcommand Tests
# ============================================================

test_suite "Sample Subcommand"

# Test: Sample help
if "$SAMSAMPLEX" sample --help 2>&1 | grep -q -i "sample"; then
    pass "sample --help shows sample usage"
else
    fail "sample --help should show sample usage"
fi

# Test: Sample without required args
if ! "$SAMSAMPLEX" sample 2>/dev/null; then
    pass "sample without args returns error"
else
    fail "sample without args should return error"
fi

# Test: Sample mode parsing
for mode in min max mean random; do
    if "$SAMSAMPLEX" sample --help 2>&1 | grep -q "$mode"; then
        pass "sample help mentions $mode mode"
    else
        fail "sample help should mention $mode mode"
    fi
done

# ============================================================
# Invalid Subcommand Tests
# ============================================================

test_suite "Invalid Commands"

# Test: Invalid subcommand
if ! "$SAMSAMPLEX" invalid_command 2>/dev/null; then
    pass "invalid subcommand returns error"
else
    fail "invalid subcommand should return error"
fi

# ============================================================
# BED File Tests (using synthetic data)
# ============================================================

test_suite "BED File Handling"

# Create a test BED file
cat > "$TMPDIR/test.bed" << 'EOF'
chr1	0	100	10
chr1	100	200	20
chr1	200	300	30
chr2	0	100	15
EOF

# Test: BED file is readable (basic syntax check)
if [ -f "$TMPDIR/test.bed" ]; then
    lines=$(wc -l < "$TMPDIR/test.bed")
    if [ "$lines" -eq 4 ]; then
        pass "Test BED file created successfully"
    else
        fail "Test BED file should have 4 lines"
    fi
else
    fail "Could not create test BED file"
fi

# ============================================================
# Summary
# ============================================================

echo ""
echo "========================================"
echo "Integration Tests Complete"
echo "========================================"
echo "Tests run:    $TESTS_RUN"
echo -e "Tests passed: ${GREEN}$TESTS_PASSED${NC}"
if [ $TESTS_FAILED -gt 0 ]; then
    echo -e "Tests failed: ${RED}$TESTS_FAILED${NC}"
    exit 1
else
    echo -e "Tests failed: $TESTS_FAILED"
    exit 0
fi

