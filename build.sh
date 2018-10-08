#!/usr/bin/env bash
cd ~/cpu_version/build
rm -rf ./bin/mdrun_mpi_d
make -j48
