cmake_minimum_required(VERSION 3.18)
project(bns_core LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(Python COMPONENTS Interpreter Development REQUIRED)
find_package(nanobind CONFIG REQUIRED)

nanobind_add_module(bns_core bns_bindings.cpp)

target_include_directories(bns_core PRIVATE
  ${CMAKE_SOURCE_DIR}/../include
)

# GSL: resolved via -DGSL_ROOT from Makefile (set by activate_spack)
if(DEFINED GSL_ROOT)
  find_library(GSL_LIB gsl
    HINTS ${GSL_ROOT}/lib
    NO_DEFAULT_PATH)
  find_library(GSLCBLAS_LIB gslcblas
    HINTS ${GSL_ROOT}/lib
    NO_DEFAULT_PATH)
else()
  find_library(GSL_LIB gsl)
  find_library(GSLCBLAS_LIB gslcblas)
endif()
target_link_libraries(bns_core PRIVATE ${GSL_LIB} ${GSLCBLAS_LIB})

# Set RPATH so the built .so finds libgsl at runtime without LD_PRELOAD
if(GSL_LIB)
  get_filename_component(GSL_LIB_DIR ${GSL_LIB} DIRECTORY)
  set_target_properties(bns_core PROPERTIES
    INSTALL_RPATH "${GSL_LIB_DIR}"
    BUILD_RPATH "${GSL_LIB_DIR}"
  )
endif()

# Disabled: BS_ASSERT checks pair/brem kernel validity but abort() 
# corrupts nanobind's Python runtime. GR-Athena++ also disables this.
# target_compile_definitions(bns_core PRIVATE BS_DEBUG)
target_compile_options(bns_core PRIVATE
  -O2 -Wall -Wextra -pedantic
  -Wno-unused-parameter -Wno-unused-variable
)
