.PHONY: conda-info conda-list setup-build setup-tests conda-rerender \
	conda-build conda-upload conda-dump-env \
	run-tests run-coveralls dev-install dev-environment help

.DEFAULT_GOAL: help

SHELL := /bin/bash
os = $(shell uname)
pkg_name = genutil
repo_name = genutil

user ?= cdat
label ?= nightly

build_script = conda-recipes/build_tools/conda_build.py

test_pkgs = testsrunner
last_stable ?= 8.2

conda_test_env = test-$(pkg_name)
conda_build_env = build-$(pkg_name)

branch ?= $(shell git rev-parse --abbrev-ref HEAD)
extra_channels ?= cdat/label/nightly conda-forge
conda ?= $(or $(CONDA_EXE),$(shell find /opt/*conda*/bin $(HOME)/*conda*/bin -type f -iname conda))

conda_env_filename ?= spec-file
build_version ?= 3.7

# Only populate if workdir is not defined
ifeq ($(workdir),)
# Create .tempdir if it doesn't exist
ifeq ($(wildcard $(PWD)/.tempdir),)
workdir := $(shell mktemp -d -t "build_$(pkg_name).XXXXXXXX")
$(shell echo $(workdir) > $(PWD)/.tempdir)
endif

# Read tempdir
workdir := $(shell cat $(PWD)/.tempdir)
endif

artif_dir = $(workdir)/$(artifact_dir)

ifneq ($(coverage),)
coverage = -c tests/coverage.json --coverage-from-egg
endif

conda_recipes_branch ?= master

conda_base = $(patsubst %/bin/conda,%,$(conda))
conda_activate = $(conda_base)/bin/activate

conda_build_extra = --copy_conda_package $(artif_dir)/

# Is this needed?
# ifndef $(local_repo)
# local_repo = $(dir $(realpath $(firstword $(MAKEFILE_LIST))))
# endif

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
dev-environment: conda_pkgs := $(gcc) "numpy>=1.18" udunits expat pytest ipython cdms2 $(test_pkgs)
dev-environment: export conda_env := dev-$(pkg_name)
dev-environment: ## Creates dev environment and installs genutil. Will need to run dev-install after any code change.
ifeq ($(os),Darwin)
	$(error dev-environment on OSX is not support)
endif

	source $(conda_activate) base; conda create -y -n $(conda_env) \
		$(conda_channels) $(conda_pkgs)

	$(MAKE) dev-install

dev-install: export conda_env := genutil-dev
dev-install: ## Installs genutil in conda environment "genutil-dev"
	source $(conda_activate) $(conda_env); \
		python setup.py build -gf; \
		python setup.py install

dev-conda-build: conda_build_extra := --local_repo $(PWD)
dev-conda-build: setup-build conda-rerender conda-build ## Build conda package

conda-info: ## Prints conda info for environment
	source $(conda_activate) base; conda info

conda-list: ## Prints packages installed in environment
	source $(conda_activate) $(conda_test_env); conda list

setup-build: ## Setup build
ifeq ($(wildcard $(workdir)/conda-recipes),)
	git clone -b $(conda_recipes_branch) https://github.com/CDAT/conda-recipes $(workdir)/conda-recipes
else
	cd $(workdir)/conda-recipes; git pull
endif

setup-tests: ## Setup test environment
	source $(conda_activate) base; conda create -y -n $(conda_test_env) --use-local \
		$(foreach x,$(extra_channels),-c $(x)) $(pkg_name) $(foreach x,$(test_pkgs),"$(x)") \
		$(foreach x,$(extra_pkgs),"$(x)")

conda-rerender: setup-build ## Rerender conda recipe using conda-smithy
	python $(workdir)/$(build_script) -w $(workdir) -l $(last_stable) -B 0 -p $(pkg_name) -r $(repo_name) \
		-b $(branch) --do_rerender --conda_env $(conda_build_env) --ignore_conda_missmatch \
		--conda_activate $(conda_activate)

conda-build: ## Builds conda recipe
	mkdir -p $(artif_dir)

	python $(workdir)/$(build_script) -w $(workdir) -p $(pkg_name) --build_version $(build_version) \
		--do_build --conda_env $(conda_build_env) --extra_channels $(extra_channels) \
		--conda_activate $(conda_activate) $(conda_build_extra)

conda-upload: ## Upload conda packages in artifcat directory
	source $(conda_activate) $(conda_build_env); \
		anaconda -t $(conda_upload_token) upload -u $(user) -l $(label) --force $(artif_dir)/*.tar.bz2

conda-dump-env: ## Dumps conda environment
	mkdir -p $(artifact_dir)

	source $(conda_activate) $(conda_test_env); conda list --explicit > $(artifact_dir)/$(conda_env_filename).txt

run-tests: ## Runs the tests using environment
	source $(conda_activate) $(conda_test_env); python run_tests.py -H -v2 $(coverage)

run-coveralls: ## Runs coveralls using environment
	source $(conda_activate) $(conda_test_env); coveralls;
