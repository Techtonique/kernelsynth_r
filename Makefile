.PHONY: clean docs build check install test help
.DEFAULT_GOAL := help

PKG := kernelsynth
VERSION := $(shell grep '^Version:' DESCRIPTION | sed 's/Version: //')

help: ## show this help
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | awk 'BEGIN {FS = ":.*?## "}; {printf "%-10s %s\n", $$1, $$2}'

docs: ## generate documentation (roxygen2)
	Rscript -e 'roxygen2::roxygenise()'

build: docs ## build source package
	cd .. && R CMD build $(PKG)

check: build ## R CMD check
	cd .. && R CMD check --as-cran --no-manual $(PKG)_$(VERSION).tar.gz

install: docs ## install package
	R CMD INSTALL .

test: ## run tests
	Rscript -e 'devtools::test()'

clean: ## remove build artifacts
	rm -rf ../$(PKG)_*.tar.gz ../$(PKG).Rcheck
