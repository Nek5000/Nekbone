#!/usr/bin/env bash

 make -f makefile.cuda.jlse clean

 make -f makefile.cuda.jlse > compile.cuda.output 2> compile.cuda.error || cat compile.cuda.error

 mpiexec -n 1 ./nekbone_cuda data > run.cuda.output 2> run.cuda.error || cat run.cuda.error

 vimdiff run.master_branch.output run.cuda.output

