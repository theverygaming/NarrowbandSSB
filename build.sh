#!/bin/bash
rm build/modulate
rm build/demodulate
make main
cd build
./modulate
./demodulate
cd ..
