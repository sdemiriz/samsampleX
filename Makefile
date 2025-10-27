.PHONY: run-example run-example-prg example example-map example-prg example-map-prg figure-1 single-interval multi-interval activate install test test-filehandler test-loader test-mapper test-intervals test-comparator test-plotter

PYTHON := python
PIP := $(PYTHON) -m pip
SR := subsample-reads
VENV_ACTIVATE := source venv/bin/activate
EXAMPLE_DIR := examples
BENCHMARKS_DIR := benchmarks
PUBLICATION_DIR := publication

install: requirements.txt
	@echo "[INFO] Creating virtual environment and installing dependencies..."
	python3 -m venv venv
	$(VENV_ACTIVATE) && $(PIP) install -r requirements.txt

example: install
	@echo "[INFO] Generating example BAM file to downsample..."
	$(VENV_ACTIVATE) && $(PYTHON) $(EXAMPLE_DIR)/make_example.py

example-map: install
	@echo "[INFO] Generating example BAM file used for mapping..."
	$(VENV_ACTIVATE) && $(PYTHON) $(EXAMPLE_DIR)/make_example_map.py

example-prg: install
	@echo "[INFO] Generating example BAM file with HLA*LA PRG to downsample..."
	$(VENV_ACTIVATE) && $(PYTHON) $(EXAMPLE_DIR)/make_example_prg.py

example-map-prg: install
	@echo "[INFO] Generating example BAM file used for mapping with HLA*LA PRG..."
	$(VENV_ACTIVATE) && $(PYTHON) $(EXAMPLE_DIR)/make_example_map_prg.py


run-example: IN_BAM = $(EXAMPLE_DIR)/example.bam
run-example: MAP_BAM = $(EXAMPLE_DIR)/example-map.bam
run-example: OUT_BAM = $(EXAMPLE_DIR)/example-out.bam
run-example: BED = $(EXAMPLE_DIR)/example-map.bed
run-example: PLT = $(EXAMPLE_DIR)/example-out.png
run-example: example example-map
	
	@echo "[STEP 1] Mapping example..."
	$(SR) map --in-bam $(MAP_BAM) --contig chr1 --start 100 --end 1100 --interval-count 10 --bed $(BED)

	@echo "[STEP 2] Sampling example..."
	$(SR) sample --in-bam $(IN_BAM) --bed $(BED) --out-bam $(OUT_BAM)

	@echo "[STEP 3] Plotting example..."
	$(SR) plot --in-bam $(IN_BAM) --map-bam $(MAP_BAM) --out-bam $(OUT_BAM) --bed $(BED) --out-plt $(PLT)

run-example-prg: IN_BAM = $(EXAMPLE_DIR)/example-prg.bam
run-example-prg: MAP_BAM = $(EXAMPLE_DIR)/example-prg-map.bam
run-example-prg: OUT_BAM = $(EXAMPLE_DIR)/example-prg-out.bam
run-example-prg: BED = $(EXAMPLE_DIR)/example-prg-map.bed
run-example-prg: PLT = $(EXAMPLE_DIR)/example-prg-out.png
run-example-prg: example-prg example-map-prg

	@echo "[STEP 1] Mapping PRG example..."
	$(SR) map --in-bam $(MAP_BAM) --contig chr6 --start 29941260 --end 29942260 --interval-count 10 --bed $(BED)

	@echo "[STEP 2] Sampling PRG example..."
	$(SR) sample --prg GRCh38 --in-bam $(IN_BAM) --bed $(BED) --out-bam $(OUT_BAM)

	@echo "[STEP 3] Plotting PRG example..."
	$(SR) plot --in-bam $(IN_BAM) --map-bam $(MAP_BAM) --out-bam $(OUT_BAM) --bed $(BED) --out-plt $(PLT)

algo-demo: IN_BAM = $(EXAMPLE_DIR)/algo-demo.bam
algo-demo: BED = $(EXAMPLE_DIR)/algo-demo.bed
algo-demo: OUT_BAM = $(EXAMPLE_DIR)/algo-demo-out.bam
algo-demo:
	@echo "[INFO] Running algorithm demo..."
	$(VENV_ACTIVATE) && $(PYTHON) $(EXAMPLE_DIR)/algo_demo.py
	$(SR) sample --in-bam $(IN_BAM) --bed $(BED) --out-bam $(OUT_BAM)
	
	@echo "[INFO] The following step will fail without samtools being installed"
	@echo "[INFO] Reads shown in Figure 1 (r2, r4, r5, r7, r10, r14, r15) should be present: "
	samtools view $(OUT_BAM) | cut -f 1 | tr '\n' ' '

download-files:
	@echo "[INFO] Downloading files for publication..."
	bash $(PUBLICATION_DIR)/download-files.sh

single-interval:
	@echo "[INFO] Running single interval benchmark..."
	bash $(BENCHMARKS_DIR)/single-interval-benchmark.sh

multi-interval:
	@echo "[INFO] Running multi interval benchmark..."
	bash $(BENCHMARKS_DIR)/multi-interval-benchmark.sh

figure-1:
	@echo "[INFO] Running figure 1 example..."
	bash $(PUBLICATION_DIR)/figure-1.sh

clean-examples:
	@echo "[INFO] Cleaning up example files..."
	rm $(EXAMPLE_DIR)/*.bam
	rm $(EXAMPLE_DIR)/*.bam.bai
	rm $(EXAMPLE_DIR)/*.bed
	rm $(EXAMPLE_DIR)/*.png

clean-benchmarks:
	@echo "[INFO] Cleaning up single-interval benchmark files..."
	rm $(BENCHMARKS_DIR)/single-interval-*.log
	rm $(BENCHMARKS_DIR)/gatk-inputs/single-interval-subset.bam*
	rm $(BENCHMARKS_DIR)/outputs/*.bam*

	@echo "[INFO] Cleaning up multi-interval benchmark files..."
	rm $(BENCHMARKS_DIR)/multi-interval-*.log
	rm -f $(BENCHMARKS_DIR)/gatk-inputs/*.bam*
	rm $(BENCHMARKS_DIR)/outputs/*.bam*
	
test-filehandler:
	@echo "[INFO] Running FileHandler tests..."
	$(VENV_ACTIVATE) && $(PYTHON) -m pytest test/test_filehandler.py -v

test-loader:
	@echo "[INFO] Running Loader tests..."
	$(VENV_ACTIVATE) && $(PYTHON) -m pytest test/test_loader.py -v

test-mapper:
	@echo "[INFO] Running Mapper tests..."
	$(VENV_ACTIVATE) && $(PYTHON) -m pytest test/test_mapper.py -v

test-intervals:
	@echo "[INFO] Running Intervals tests..."
	$(VENV_ACTIVATE) && $(PYTHON) -m pytest test/test_intervals.py -v

test-comparator:
	@echo "[INFO] Running Comparator tests..."
	$(VENV_ACTIVATE) && $(PYTHON) -m pytest test/test_comparator.py -v

test-plotter:
	@echo "[INFO] Running Plotter tests..."
	$(VENV_ACTIVATE) && $(PYTHON) -m pytest test/test_plotter.py -v

test:
	@echo "[INFO] Running all tests..."
	$(VENV_ACTIVATE) && $(PYTHON) -m pytest test/ -v

