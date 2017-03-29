#!/usr/bin/env bash

 make -f makefile.mpi.jlse  clean
 make -f makefile.acc.jlse clean

 make -f makefile.mpi.jlse  > compile.mpi.stdout  2> compile.mpi.stderr  || cat compile.mpi.stderr
 make -f makefile.acc.jlse > compile.acc.stdout 2> compile.acc.stderr || cat compile.acc.stderr

 mpiexec -n 1 ./nekbone_mpi data  > run.mpi.stdout  2> run.mpi.stderr  || cat run.mpi.stderr
 mpiexec -n 1 ./nekbone_acc data > run.acc.stdout 2> run.acc.stderr || cat run.acc.stderr

 vimdiff run.mpi.stdout run.acc.stdout

