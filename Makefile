# samsampleX-c Makefile
# Requires: htslib (install via apt, brew, or from source)

CC ?= gcc
CFLAGS = -Wall -Wextra -O2 -g
PREFIX ?= /usr/local

# Program name
PROG = samsampleX

# Directories
SRC_DIR = src
INC_DIR = include
BUILD_DIR = build

# Source files
SRCS = $(SRC_DIR)/main.c \
       $(SRC_DIR)/map.c \
       $(SRC_DIR)/sample.c \
       $(SRC_DIR)/bed.c \
       $(SRC_DIR)/depth.c \
       $(SRC_DIR)/metrics.c

# Object files
OBJS = $(SRCS:$(SRC_DIR)/%.c=$(BUILD_DIR)/%.o)

# Include paths
INCLUDES = -I$(INC_DIR) -I$(SRC_DIR)

# Check for htslib using pkg-config (with --static for full dependencies)
HTSLIB_CFLAGS := $(shell pkg-config --cflags htslib 2>/dev/null)
HTSLIB_LIBS := $(shell pkg-config --libs --static htslib 2>/dev/null)

# Fallback if pkg-config doesn't find htslib
ifeq ($(HTSLIB_LIBS),)
    # Try common locations
    ifneq ($(wildcard /usr/include/htslib/hts.h),)
        HTSLIB_CFLAGS =
        HTSLIB_LIBS = -lhts
    else ifneq ($(wildcard /usr/local/include/htslib/hts.h),)
        HTSLIB_CFLAGS = -I/usr/local/include
        HTSLIB_LIBS = -L/usr/local/lib -lhts
    else ifneq ($(wildcard $(HOME)/local/include/htslib/hts.h),)
        HTSLIB_CFLAGS = -I$(HOME)/local/include
        HTSLIB_LIBS = -L$(HOME)/local/lib -lhts
    else
        $(error htslib not found. Install with: apt install libhts-dev, brew install htslib, or build from source)
    endif
endif

# Additional libraries that htslib may depend on (for static linking)
# These are commonly needed on HPC clusters with static htslib builds
EXTRA_LIBS = -lm -lpthread -lz -lbz2 -llzma -lcurl -lssl -lcrypto -ldeflate

CFLAGS += $(HTSLIB_CFLAGS) $(INCLUDES)
LDFLAGS += $(HTSLIB_LIBS) $(EXTRA_LIBS)

# Default target
all: check-htslib $(BUILD_DIR) $(PROG)

# Check htslib is available
check-htslib:
	@echo "Checking for htslib..."
	@echo "  HTSLIB_CFLAGS: $(HTSLIB_CFLAGS)"
	@echo "  HTSLIB_LIBS: $(HTSLIB_LIBS)"

# Create build directory
$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

# Link program
$(PROG): $(OBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)
	@echo "Built $(PROG) successfully"

# Compile source files
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.c | $(BUILD_DIR)
	$(CC) $(CFLAGS) -c -o $@ $<

# Install
install: $(PROG)
	install -d $(DESTDIR)$(PREFIX)/bin
	install -m 755 $(PROG) $(DESTDIR)$(PREFIX)/bin/

# Uninstall
uninstall:
	rm -f $(DESTDIR)$(PREFIX)/bin/$(PROG)

# Clean
clean:
	rm -rf $(BUILD_DIR) $(PROG) tests/test_parsing tests/test_depth tests/test_metrics

# ============================================================
# Testing
# ============================================================

TEST_DIR = tests

# Test executables
TEST_PARSING = $(TEST_DIR)/test_parsing
TEST_DEPTH = $(TEST_DIR)/test_depth
TEST_METRICS = $(TEST_DIR)/test_metrics

# Test object files (reuse some from main build)
TEST_OBJS_COMMON = $(BUILD_DIR)/depth.o $(BUILD_DIR)/bed.o

# Build test executables
$(TEST_PARSING): $(TEST_DIR)/test_parsing.c $(BUILD_DIR)/depth.o $(BUILD_DIR)/bed.o | $(BUILD_DIR)
	$(CC) $(CFLAGS) -o $@ $< $(BUILD_DIR)/depth.o $(BUILD_DIR)/bed.o $(LDFLAGS)

$(TEST_DEPTH): $(TEST_DIR)/test_depth.c $(BUILD_DIR)/bed.o | $(BUILD_DIR)
	$(CC) $(CFLAGS) -o $@ $< $(BUILD_DIR)/bed.o $(LDFLAGS)

$(TEST_METRICS): $(TEST_DIR)/test_metrics.c $(BUILD_DIR)/metrics.o $(BUILD_DIR)/bed.o | $(BUILD_DIR)
	$(CC) $(CFLAGS) -o $@ $< $(BUILD_DIR)/metrics.o $(BUILD_DIR)/bed.o $(LDFLAGS)

# Run unit tests
test-unit: $(TEST_PARSING) $(TEST_DEPTH) $(TEST_METRICS)
	@echo ""
	@echo "Running unit tests..."
	@echo "========================================"
	@$(TEST_PARSING) && $(TEST_DEPTH) && $(TEST_METRICS)

# Run integration tests (requires main program to be built)
test-integration: $(PROG)
	@echo ""
	@echo "Running integration tests..."
	@echo "========================================"
	@$(TEST_DIR)/integration_tests.sh

# Run all tests
test: test-unit test-integration
	@echo ""
	@echo "All tests completed successfully!"

# Phony targets
.PHONY: all check-htslib install uninstall clean test test-unit test-integration

# Dependencies (auto-generated would be better, but keep it simple)
$(BUILD_DIR)/main.o: $(SRC_DIR)/main.c $(INC_DIR)/samsampleX.h $(SRC_DIR)/map.h $(SRC_DIR)/sample.h
$(BUILD_DIR)/map.o: $(SRC_DIR)/map.c $(SRC_DIR)/map.h $(INC_DIR)/samsampleX.h $(SRC_DIR)/bed.h $(SRC_DIR)/depth.h
$(BUILD_DIR)/sample.o: $(SRC_DIR)/sample.c $(SRC_DIR)/sample.h $(INC_DIR)/samsampleX.h $(SRC_DIR)/bed.h $(SRC_DIR)/depth.h $(SRC_DIR)/metrics.h
$(BUILD_DIR)/bed.o: $(SRC_DIR)/bed.c $(SRC_DIR)/bed.h $(INC_DIR)/samsampleX.h
$(BUILD_DIR)/depth.o: $(SRC_DIR)/depth.c $(SRC_DIR)/depth.h $(INC_DIR)/samsampleX.h
$(BUILD_DIR)/metrics.o: $(SRC_DIR)/metrics.c $(SRC_DIR)/metrics.h $(INC_DIR)/samsampleX.h

