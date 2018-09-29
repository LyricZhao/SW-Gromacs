#!/usr/bin/env bash
cd sw_slave_kernel
sw5cc -slave -c -I../gromacs-5.1.5/src -I. -I../gromacs-5.1.5/src/gromacs/mdlib/nbnxn_kernels -I../gromacs-5.1.5/src/external/thread_mpi/include -msimd -DGMX_DOUBLE sw_slave_kernel.c -o ../../sw_build/sw_slave_kernel.o
sw5cc -host  -c -msimd sw_cpe_extern.c -o ../../sw_build/sw_cpe_kernel.o
sw5ar cr ../../sw_build/libmdrun-sw.a ../../sw_build/sw_slave_kernel.o ../../sw_build/sw_cpe_kernel.o
