#!/bin/bash

# clean compile
rm -rf build
mkdir -p build
cd build

cmake -DENABLE_KOKKOS=OFF ..
cmake --build . --target MWE_profile

# run
./MWE_profile 100

cd ..
