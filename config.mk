# config.mk -- resolve Python and GSL from environment
# Source activate_conda + activate_spack before using.

# GSL: resolve from spack num_gcc environment
ifdef SPACK_ENV
  GSL_ROOT := $(SPACK_ENV)/.spack-env/view
  GSL_INCLUDE := -I$(GSL_ROOT)/include
  GSL_LIBS    := -L$(GSL_ROOT)/lib -lgsl -lgslcblas
else
  GSL_INCLUDE := $(shell pkg-config --cflags gsl 2>/dev/null)
  GSL_LIBS    := $(shell pkg-config --libs gsl gslcblas 2>/dev/null)
endif

# Python for bindings: prefer $CONDA_PREFIX from conda activate, fall back to PATH
ifdef CONDA_PREFIX
  PYTHON = $(CONDA_PREFIX)/bin/python3
else
  PYTHON := $(shell which python3 2>/dev/null)
endif

# nanobind cmake dir
NANOBIND_CMAKE_DIR := $(shell $(PYTHON) -c "import nanobind; print(nanobind.cmake_dir())" 2>/dev/null)
