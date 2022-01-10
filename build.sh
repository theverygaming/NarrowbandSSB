#!/bin/bash
rm build/decimate
make main
cd build
./decimate
cd ..
