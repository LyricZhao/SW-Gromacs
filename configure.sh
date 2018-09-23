#!/usr/bin/env bash
cd ~/zhaocg/gromacs-5.1.5/
mkdir build
cd build
LD=mpiicc CC=mpiicc CXX=mpiicpc cmake .. -DGMX_FFT_LIBRARY=fftpack -DGMX_MPI=on -DGMX_DOUBLE=ON -DGMX_BUILD_MDRUN_ONLY=ON -DBUILD_SHARED_LIBS=off -LH
