#!/usr/bin/env bash

 make -f makefile.mpi.jlse clean
 make -f makefile.acc.jlse clean

 make -f makefile.mpi.jlse > compile.mpi.output 2> compile.mpi.error || cat compile.mpi.error
 make -f makefile.acc.jlse > compile.acc.output 2> compile.acc.error || cat compile.acc.error

 mpiexec -n 1 ./nekbone_mpi data > run.mpi.output 2> run.mpi.error || cat run.mpi.error
 mpiexec -n 1 ./nekbone_acc data > run.acc.output 2> run.acc.error || cat run.acc.error

 vimdiff run.mpi.output run.acc.output

