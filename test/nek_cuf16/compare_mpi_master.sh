#!/usr/bin/env bash

 make -f makefile.mpi.jlse clean

 make -f makefile.mpi.jlse > compile.mpi.output 2> compile.mpi.error || cat compile.mpi.error

 mpiexec -n 1 ./nekbone_mpi data > run.mpi.output 2> run.mpi.error || cat run.mpi.error

 vimdiff run.master_branch.output run.mpi.output 

