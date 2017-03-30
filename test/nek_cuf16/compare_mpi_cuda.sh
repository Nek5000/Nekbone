#!/usr/bin/env bash

 make -f makefile.mpi.jlse  clean
 make -f makefile.cuda.jlse clean

 make -f makefile.mpi.jlse  > compile.mpi.output  2> compile.mpi.error  || cat compile.mpi.error
 make -f makefile.cuda.jlse > compile.cuda.output 2> compile.cuda.error || cat compile.cuda.error

 mpiexec -n 1 ./nekbone_mpi data  > run.mpi.output  2> run.mpi.error  || cat run.mpi.error
 mpiexec -n 1 ./nekbone_cuda data > run.cuda.output 2> run.cuda.error || cat run.cuda.error

 vimdiff run.mpi.output run.cuda.output

