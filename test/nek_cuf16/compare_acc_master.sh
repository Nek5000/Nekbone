#!/usr/bin/env bash

 make -f makefile.acc.jlse clean

 make -f makefile.acc.jlse > compile.acc.output 2> compile.acc.error || cat compile.acc.error

 mpiexec -n 1 ./nekbone_acc data > run.acc.output 2> run.acc.error || cat run.acc.error

 vimdiff run.master_branch.output run.acc.output

