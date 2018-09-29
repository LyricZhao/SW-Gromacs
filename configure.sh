#!/usr/bin/env bash
cd ~/zhaocg/
rm -rf build
mkdir build
cd build
LD=mpiCC CC=mpicc CXX=mpiCC cmake ../code/gromacs-5.1.5 -DGMX_OPENMP=off -DGMX_FFT_LIBRARY=fftpack -DGMX_MPI=ON -DGMX_DOUBLE=ON -DGMX_BUILD_MDRUN_ONLY=ON -DBUILD_SHARED_LIBS=off -LH
