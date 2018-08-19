#! /usr/bin/env sh
g++ -c pmepot.cpp
gfortran -c pmepotFDriver.f03
gfortran -o fdriver pmepotFDriver.o pmepot.o -lfftw3 -lstdc++
./fdriver
