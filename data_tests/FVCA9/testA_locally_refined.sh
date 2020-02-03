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
tcdiff=3

# Solver
solver_type="bicgstab"

# Meshes (without .typ2)
mesh[1]="mesh3_2"
mesh[2]="mesh3_3"
mesh[3]="mesh3_4"
mesh[4]="mesh3_5"


