.PHONY: conda-info conda-list setup-build setup-tests conda-rerender \
	conda-build conda-upload conda-dump-env \
	run-tests run-coveralls dev-install dev-environment help

.DEFAULT_GOAL: help

SHELL := /bin/bash
os = $(shell uname)
pkg_name = genutil
repo_name = genutil
build_script = conda-recipes/build_tools/conda_build.py

test_pkgs = testsrunner
last_stable ?= 8.2

conda_env ?= genutil-build
workdir ?= $(PWD)/workspace
branch ?= $(shell git rev-parse --abbrev-ref HEAD)
extra_channels ?= cdat/label/nightly conda-forge
conda ?= $(or $(CONDA_EXE),$(shell find /opt/*conda*/bin $(HOME)/*conda*/bin -type f -iname conda))
artifact_dir ?= $(PWD)/artifacts
conda_env_filename ?= spec-file
build_version ?= 3.7

ifneq ($(coverage),)
coverage = -c tests/coverage.json --coverage-from-egg
endif

conda_recipes_branch ?= master

conda_base = $(patsubst %/bin/conda,%,$(conda))
conda_activate = $(conda_base)/bin/activate

conda_build_extra = --copy_conda_package $(artifact_dir)/

ifndef $(local_repo)
local_repo = $(dir $(realpath $(firstword $(MAKEFILE_LIST))))
endif

help: ## Prints help
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | \
		sort | awk 'BEGIN {FS = ":.*?## "}; {printf "\033[36m%-30s\033[0m %s\n", $$1, $$2}'

dev-docker:
	docker run -d --name genutil-dev -v $(PWD):/src -w /src continuumio/miniconda3 /bin/sleep infinity || exit 0
	docker start genutil-dev
	docker exec -it genutil-dev /bin/bash -c "apt update; apt install -y make"
	docker exec -it genutil-dev /bin/bash -c "make dev-environment"
	docker exec -it genutil-dev /bin/bash -c "conda init bash; echo 'conda activate genutil-dev' >> ~/.bashrc"
	docker exec -it genutil-dev /bin/bash

dev-environment: conda_channels := -c conda-forge
dev-environment: gcc := $(or $(if $(findstring os,Darwin),clang_osx-64), gcc_linux-64)
dev-environment: conda_pkgs := $(gcc) "numpy>=1.18" udunits expat pytest ipython cdms2
dev-environment: export conda_env := genutil-dev
dev-environment: ## Creates dev environment and installs genutil in "Development Mode". If you modify c code you will need to run "make dev-install".
ifeq ($(os),Darwin)
	$(error dev-environment on OSX is not support)
endif

	source $(conda_activate) base; conda create -y -n $(conda_env) \
		$(conda_channels) $(conda_pkgs)

	$(MAKE) dev-install

dev-install: export conda_env := genutil-dev
dev-install: ## Installs genutil in "Development Mode", will recompile c code each invocation."
	source $(conda_activate) $(conda_env); \
		python setup.py build -gf; \
		python setup.py install

conda-info: ## Prints conda info for environment
	source $(conda_activate) $(conda_env); conda info

conda-list: ## Prints packages installed in environment
	source $(conda_activate) $(conda_env); conda list

setup-build: ## Setup build
ifeq ($(wildcard $(workdir)/conda-recipes),)
	git clone -b $(conda_recipes_branch) https://github.com/CDAT/conda-recipes $(workdir)/conda-recipes
else
	cd $(workdir)/conda-recipes; git pull
endif

setup-tests: ## Setup test environment
	source $(conda_activate) base; conda create -y -n $(conda_env) --use-local \
		$(foreach x,$(extra_channels),-c $(x)) $(pkg_name) $(foreach x,$(test_pkgs),"$(x)") \
		$(foreach x,$(extra_pkgs),"$(x)")

conda-rerender: setup-build ## Rerender conda recipe using conda-smithy
	python $(workdir)/$(build_script) -w $(workdir) -l $(last_stable) -B 0 -p $(pkg_name) -r $(repo_name) \
		-b $(branch) --do_rerender --conda_env $(conda_env) --ignore_conda_missmatch \
		--conda_activate $(conda_activate)

conda-build: ## Builds conda recipe
	mkdir -p $(artifact_dir)

	python $(workdir)/$(build_script) -w $(workdir) -p $(pkg_name) --build_version $(build_version) \
		--do_build --conda_env $(conda_env) --extra_channels $(extra_channels) \
		--conda_activate $(conda_activate) $(conda_build_extra)

conda-upload: ## Upload conda packages in artifcat directory
	source $(conda_activate) $(conda_env); \
		anaconda -t $(conda_upload_token) upload -u $(user) -l $(label) --force $(artifact_dir)/*.tar.bz2

conda-dump-env: ## Dumps conda environment
	mkdir -p $(artifact_dir)

	source $(conda_activate) $(conda_env); conda list --explicit > $(artifact_dir)/$(conda_env_filename).txt

run-tests: ## Runs the tests using environment
	source $(conda_activate) $(conda_env); python run_tests.py -H -v2 $(coverage)

run-coveralls: ## Runs coveralls using environment
	source $(conda_activate) $(conda_env); coveralls;
