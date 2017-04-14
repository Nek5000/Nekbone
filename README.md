Nekbone
=======

Nekbone solves a standard Poisson equation using a conjugate gradient iteration
with a simple or spectral element multigrid preconditioner on a block or linear
geometry. It exposes the principal computational kernel to reveal the essential
elements of the algorithmic- architectural coupling that is pertinent to
Nek5000.

CUDA/OpenACC Branch
-------------------

The CUDA/OpenACC branch contains GPU implementations of the conjugate gradient
solver. This includes a pure OpenACC implementation as well as a hybrid
OpenACC/CUDA implementation with a CUDA kernel for matrix-vector multiplication,
based on previous works [1][2]. This implementation can also be compiled with 
without OpenACC or CUDA.  These implementations were tested with the nek\_gpu1
example in `tests/nek_gpu1`.  They were also tested with the PGI compilers.  

### Compiling the nek\_gpu1 example

The nek\_gpu1 example contains three makenek scripts to compile Nekbone with
various configurations:

* `makenek.mpi`: MPI only, without OpenACC or CUDA
* `makenek.acc`: MPI and OpenACC
* `makenek.cuda`: MPI and hybrid OpenACC/CUDA

The different build configurations are determined by the presence/absence of
the -acc and -Mcuda flags, as demonstrated in the scripts.  These scripts will
work without out-of-the-box provided that:

* `mpicc`, `mpif77`, and `mpif90` set to suitable compilers in your environment
  (i.e., MPI compilers that use PGI or Cray)
* You are compiling for a compute-capability 5.0 or 6.0 GPU

If these defaults are not suitable, you may manually edit the relevant settings
in script.  You may also compile without MPI by setting the `IFMPI` option in
makenek.

To use the compilation scripts (for example, with OpenACC only):

```
$ cd test/nek_gpu1
$ ./makenek.acc clean
$ ./makenek.acc
```

The clean step is necessary if you are switching between MPI, OpenACC, and
OpenACC/CUDA builds.  

### Running nek\_gpu1

Successful compilation will produce an executable called `nekbone`.  You may
run the executable through an MPI process manager, for example:

```
$ mpiexec -n 1 ./nekbone data
```

The `data` argument specifies that you will be using the runtime settings from
`data.rea`.  See `USERGUIDE.pdf` for a detailed description of the `.rea` file
format.

### Testing nek\_gpu1

The `tests/nek_gpu1/scripts/` folder contains simple scripts for comparing
results from the different implementations:
* `compare_mpi_master.sh`: Runs and compares MPI solution to reference solution
  from the Nekbone master branch
* `compare_acc_master.sh`: Runs and compares OpenACC solution to reference
  solution 
* `compare_cuda_master.sh`: Runs and compares OpenACC/CUDA solution to
  reference solution
* `compare_mpi_acc.sh`: Runs and compares MPI and OpenACC solutions
* `compare_mpi_cuda.sh`: Runs and compares MPI and CUDA solutions

To use them:

```
$ cd test/nek_gpu1
$ scripts/compare_acc_master.sh
```

### Modifying nek\_gpu1

The `nek_gpu1` problem is setup for 256 elements and a polynomial order of 16.
You may alter this setup in the SIZE` and `data.rea` files as described in
`USERGUIDE.pdf`.  

### Compile/Runs on Other Systems 

* Compile/Run for MPI only on Titan (16 MPI ranks, 1 node):

```
$ cd test/nek_gpu1
$ cp ../../bin/makenek.mpi.titan .
$ cp ../../bin/nekpmpi.titan .
$ ./makenek.mpi.titan
$ ./nekpmpi.titan data 16 1
```

* Compile/Run for MPI and OpenACC on Titan (1 GPU, 1 node): 

```
$ cd test/nek_gpu1
$ cp ../../bin/makenek.acc.titan .
$ cp ../../bin/nekpgpu.titan .
$ ./makenek.acc.titan
$ ./nekpgpu.titan data 1 1
```

* Compile/Run for MPI and hybrid OpenACC/CUDA on Titan (1 GPU, 1 node): 

```
$ cd test/nek_gpu1
$ cp ../../bin/makenek.cuda.titan .
$ cp ../../bin/nekpgpu.titan .
$ ./makenek.cuda.titan
$ ./nekpgpu.titan data 1 1
```
---

### References

[1] Jing Gong, Stefano Markidis, Erwin Laure, Matthew Otten, Paul Fischer, Misun Min,   
    Nekbone performance on GPUs with OpenACC and CUDA Fortran implementations, 
    The jounral of Supercomputing, Vol. 72, pp. 4160-4180, 2016. 

[2] Matthew Otten, Jing Gong, Azamat Mametjanov, Aaron Vose, John Levesque, Paul  
    Fischer, and Misun Min, An MPI/OpenACC implementation of a high order 
    electromagnetics solver with GPUDirect communication, The International Journal 
    of High Performance Computing Application, Vol. 30, No. 3, pp. 320-334, 2016. 

