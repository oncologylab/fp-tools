.PHONY: test build release-check paper-pdf paper-smoke clean-paper

PYTHON ?= .venv/bin/python
LATEXMK ?= latexmk

test:
	$(PYTHON) -m unittest discover -s tests -v

build:
	cd /tmp && $(abspath $(PYTHON)) -m build $(CURDIR) --outdir $(CURDIR)/dist

release-check: build
	$(PYTHON) -m twine check dist/*

paper-pdf:
	cd manuscript && $(LATEXMK) -pdf -shell-escape -interaction=nonstopmode main.tex

paper-smoke:
	$(PYTHON) manuscript/scripts/plot_denovo_motif_validation.py
	cd manuscript && $(LATEXMK) -pdf -shell-escape -interaction=nonstopmode main.tex

clean-paper:
	cd manuscript && $(LATEXMK) -c
