#!/usr/bin/env bash
./slave_build.sh
cd ~/zhaocg/
rm -rf build
mkdir build
cd build
LD=mpiCC CC=mpicc CXX=mpiCC LDFLAGS="-Wl,--whole-archive,/home/export/online1/cpc051/zhaocg/sw_build/libmdrun-sw.a,--no-whole-archive -Wl,--whole-archive,-wrap,athread_init,-wrap,__expt_handler,-wrap,__real_athread_spawn /home/export/online1/swmore/release/lib/libspc.a -Wl,--no-whole-archive -lm_slave" cmake ../code/gromacs-5.1.5 -DGMX_FFT_LIBRARY=fftpack -DGMX_MPI=ON -DGMX_DOUBLE=ON -DGMX_BUILD_MDRUN_ONLY=ON -DBUILD_SHARED_LIBS=off -LH
# LD=mpiCC CC=mpicc CXX=mpiCC cmake ../code/gromacs-5.1.5 -DGMX_FFT_LIBRARY=fftpack -DGMX_MPI=ON -DGMX_DOUBLE=ON -DGMX_BUILD_MDRUN_ONLY=ON -DBUILD_SHARED_LIBS=off -LH
