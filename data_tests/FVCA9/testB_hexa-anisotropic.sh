#!/bin/bash
#

# Edge and cell degrees
kmax=3
lmax=3

# Choice basis
#choice_basis="Mon"
choice_basis="ON"

# Boundary conditions (0=Dirichlet, 1=Neumann)
bc=0

# Test case
tcsol=2
tcdiff=1

# Solver
solver_type="bicgstab"

mesh[1]="anisotropic/hexa20x20"
mesh[2]="anisotropic/hexa40x80"
mesh[3]="anisotropic/hexa80x320"
mesh[4]="anisotropic/hexa120x720"

