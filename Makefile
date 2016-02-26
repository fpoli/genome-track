
linter:
	@pep8 .

test:
	@echo "=== Python tests..."
	@nosetests-3.4

	@echo "=== Compare output (1/2)..."
	@./gtrack.py sequence_coverage \
		-i sample/alignments.sam \
		-l 123000 \
		| tail -n +2 | md5sum
	@./gtrack.py sequence_coverage_lowmem \
		-i sample/alignments.sam \
		-l 123000 \
		| tail -n +2 | md5sum

	@echo "=== Compare output (2/2)..."
	@./gtrack.py physical_coverage \
		-i sample/alignments.sam \
		-l 123000 \
		| tail -n +2 | md5sum
	@./gtrack.py physical_coverage_lowmem \
		-i sample/alignments.sam \
		-l 123000 \
		| tail -n +2 | md5sum

	@echo "=== Done. Checksums should be equal."

clean:
	@find . -name "*.pyc" -delete
	@find . -name "__pycache__" -type d -exec rm -r {} +

.PHONY: linter test
