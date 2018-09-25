#!/usr/bin/env bash
bsub -I -b -q q_sw_cpc -cgsp 64 -n 512 -np 4 -share_size 6500 -host_stack 500 -J GROMACS_LYRIC ~/zhaocg/gromacs-5.1.5/build/bin/mdrun_mpi_d -s ~/zhaocg/origin/ion_channel-st.tpr -v
