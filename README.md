Nekbone
=======

Nekbone captures the basic structure and user interface of the ex-
tensive Nek5000 software. Nek5000 is a high order, incompressible
Navier-Stokes solver based on the spectral element method. It has
a wide range of applications and intricate customizations available
to users. Nekbone, on the other hand, solves a Helmholtz equation
in a box, using the spectral element method. It is pared down to
include only the necessary features to compile, run, and solve the
applications found in the test/ directory. Since almost all prac-
tical applications are in the three dimensional space, the solver is
set to work with three dimensional geometries as default. Nekbone
solves a standard Poisson equation using a conjugate gradient it-
eration with a simple or spectral element multigrid preconditioner
on a block or linear geometry (parameters are set within the test
directory of the simulation). Nekbone exposes the principal compu-
tational kernel to reveal the essential elements of the algorithmic-
architectural coupling that is pertinent to Nek5000.
