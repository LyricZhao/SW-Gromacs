#!/usr/bin/env bash
./slave_build.sh
cd ~/zhaocg/build
rm -rf ./bin/mdrun_mpi_d
make -j48
