#!/usr/bin/env bash

 make -f makefile.mpi.jlse  clean
 make -f makefile.cuda.jlse clean

 make -f makefile.mpi.jlse  > compile.mpi.stdout  2> compile.mpi.stderr  || cat compile.mpi.stderr
 make -f makefile.cuda.jlse > compile.cuda.stdout 2> compile.cuda.stderr || cat compile.cuda.stderr

 mpiexec -n 1 ./nekbone_mpi data  > run.mpi.stdout  2> run.mpi.stderr  || cat run.mpi.stderr
 mpiexec -n 1 ./nekbone_cuda data > run.cuda.stdout 2> run.cuda.stderr || cat run.cuda.stderr

 vimdiff run.mpi.stdout run.cuda.stdout

