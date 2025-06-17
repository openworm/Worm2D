#!/bin/bash
set -ex

quick_test=0

if [[ ($# -eq 1) && ($1 == '-q') ]]; then
    quick_test=1
fi

cd src/CE_locomotion/

make clean
make tests

rm -rf test_output/*.dat
./tests

make

echo "Finished running Worm2D tests"
