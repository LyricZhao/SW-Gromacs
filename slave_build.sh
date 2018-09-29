#!/usr/bin/env bash
cd sw_slave_kernel
sw5cc -slave -c -I../gromacs-5.1.5/src -I../gromacs-5.1.5/src/gromacs/mdlib/nbnxn_kernels -msimd sw_slave_kernel.c -o ../../sw_build/sw_slave_kernel.o
sw5ar cr ../../sw_build/libmdrun-sw.a ../../sw_build/sw_slave_kernel.o
