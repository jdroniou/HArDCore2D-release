#!/bin/bash
#

# In this test case, the number of elements is doubled at each iteration, and the anisotropy of the 
# cells in the lower half of the domain is also doubled.

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


mesh[1]="anisotropic/cart25_a1"
mesh[2]="anisotropic/cart50_a2"
mesh[3]="anisotropic/cart100_a4"
mesh[4]="anisotropic/cart200_a8"


