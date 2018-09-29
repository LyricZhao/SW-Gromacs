#!/usr/bin/env bash
bsub -I -b -q q_sw_cpc -cgsp 64 -n 64 -np 4 -share_size 6500 -host_stack 500 -J GROMACS_LYRIC /home/export/online1/cpc051/zhaocg/build/bin/mdrun_mpi_d -s /home/export/online1/cpc051/cases/lignocellulose-rf.BGQ-st.tpr -v
