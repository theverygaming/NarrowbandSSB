#!/bin/bash
make clean
make
cd build
./modulate
./demodulate
cd ..
