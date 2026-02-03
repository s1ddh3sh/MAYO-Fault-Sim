#!/usr/bin/env bash
set -euo pipefail

rm -rf build
mkdir build

cd build || { echo "build directory not found"; exit 1; }
cmake -DMAYO_VARIANT=MAYO_1 ..
make
./apps/test_mayo_mayo_1
