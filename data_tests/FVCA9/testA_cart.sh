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
tcdiff=3

# Solver
solver_type="bicgstab"
#solver_type="ma41"

# Meshes (without .typ2)
mesh[1]="cart10x10"
mesh[2]="cart20x20"
mesh[3]="cart40x40"
mesh[4]="cart80x80"


