#!/bin/bash
#

# Edge and cell degrees
kmax=3
lmax=3

# Use threads
use_threads="true"

# Export matrix
export_matrix="true"

# Boundary conditions (0=Dirichlet, 1=Neumann)
bc="D"

# Test case
tcsol=2
tcdiff=5

# Solver
solver_type="bicgstab"

# Meshes (without .typ2)
mesh[1]="anisotropic/hexa_straight10x10"
mesh[2]="anisotropic/hexa_straight20x40"
mesh[3]="anisotropic/hexa_straight40x160"
mesh[4]="anisotropic/hexa_straight80x640"

