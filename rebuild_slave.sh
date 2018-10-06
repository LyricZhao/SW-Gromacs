#!/usr/bin/env bash
echo 'rebuilding slave library ... ok!'
rm ../build/bin/mdrun_mpi_d
./slave_build.sh
./build.sh
