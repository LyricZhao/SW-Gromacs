#!/usr/bin/env bash
cd sw_slave_kernel
sw5cc -slave -c -I../gromacs-5.1.5/src -I../gromacs-5.1.5/src/gromacs/mdlib/nbnxn_kernels -msimd sw_slave_kernel.c -o sw_slave_kernel.o
sw5ar cr libmdrun-sw.a sw_slave_kernel.o

# CC=mpicc CXX=mpiCC LDFLAGS="-Wl,--whole-archive,`readlink -f cpelib`/libmdrun-sw.a,--no-whole-archive" cmake .. -DGMX_FFT_LIBRARY=fftpack -DGMX_MPI=on -DGMX_DOUBLE=ON -DGMX_BUILD_MDRUN_ONLY=ON -DBUILD_SHARED_LIBS=off -LH
