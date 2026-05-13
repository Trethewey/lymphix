# Lymphix — build and test helpers.
#
# Override REGISTRY and TAG when pushing:
#     make build push REGISTRY=ghcr.io/myorg/lymphix TAG=0.1.0

REGISTRY ?= ghcr.io/trethewey/lymphix
TAG      ?= 0.1.0
IMAGES   = fastp fgbio trust4 igblast clonality

.PHONY: help build push test test-unit test-smoke clean

help:
	@echo "Lymphix targets:"
	@echo "  make build              Build all $(words $(IMAGES)) container images locally"
	@echo "  make build-IMAGE        Build one image (fastp|fgbio|trust4|igblast|clonality)"
	@echo "  make push               Push all images to \$$REGISTRY"
	@echo "  make test               Run unit + smoke tests"
	@echo "  make test-unit          Python unit tests for clonality metrics"
	@echo "  make test-smoke         End-to-end mock pipeline (no Docker needed)"
	@echo "  make clean              Remove work/ and results/ artefacts"

build: $(addprefix build-,$(IMAGES))

build-%:
	docker build -f containers/$*/Dockerfile -t $(REGISTRY)/$*:$(TAG) .

push: $(addprefix push-,$(IMAGES))

push-%:
	docker push $(REGISTRY)/$*:$(TAG)

test: test-unit test-smoke

test-unit:
	@command -v pytest >/dev/null 2>&1 || { echo "pytest not found; pip install pytest"; exit 1; }
	pytest tests/test_clonality_metrics.py -v

test-smoke:
	bash tests/test_smoke.sh

clean:
	rm -rf work/ results/ results_*/ .nextflow* __pycache__/
