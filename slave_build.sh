#!/usr/bin/env bash
my_flags='-DBOOST_NO_TYPEID -DGMX_DOUBLE -DHAVE_CONFIG_H -Wundef -Wall -Wno-unused -Wunused-value -Wunused-parameter -Wno-unknown-pragmas -O3 -DNDEBUG -funroll-all-loops   -I/home/export/base/CPC/cpc051/zhaocg/build/src/external/tng_io/include -I/home/export/base/CPC/cpc051/zhaocg/code/gromacs-5.1.5/src/external/tng_io/include -isystem /home/export/base/CPC/cpc051/zhaocg/code/gromacs-5.1.5/src/external/boost -I/home/export/base/CPC/cpc051/zhaocg/build/src -I/home/export/base/CPC/cpc051/zhaocg/code/gromacs-5.1.5/src/external/thread_mpi/include -I/home/export/base/CPC/cpc051/zhaocg/code/gromacs-5.1.5/src -I/home/export/base/CPC/cpc051/zhaocg/code/gromacs-5.1.5/src/gromacs/mdlib/nbnxn_kernels/ -I./'
cd sw_slave_kernel
sw5cc -slave -c $my_flags sw_slave_kernel.c -o ../../sw_build/sw_slave_kernel.o
sw5cc -host  -c $my_flags sw_cpe_extern.c -o ../../sw_build/sw_cpe_kernel.o
sw5ar cr ../../sw_build/libmdrun-sw.a ../../sw_build/sw_slave_kernel.o ../../sw_build/sw_cpe_kernel.o
