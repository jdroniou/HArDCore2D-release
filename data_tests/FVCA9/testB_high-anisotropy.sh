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
#solver_type="ma41"

# Meshes (without .typ2)
mesh[1]="anisotropic/cart50_a10"
mesh[2]="anisotropic/cart200_a20"
mesh[3]="anisotropic/cart400_a50"


