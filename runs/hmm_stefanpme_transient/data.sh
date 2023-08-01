#!/bin/bash
#

# Use threads
use_threads="true"

# Boundary conditions (D, N, Mx)
bc="D"  # N may not work...

# Test case
tcsol=1
tcdiff=1

# Type of nonlinearity
tcNL=2
mPME=2

# Weight of mass-lumping on edges
weight=0

# Source: exact (1) or zero (0)
source=1

# Final time and initial time step
FinalTps=1
dt_initial=0.1

# Solver
solver_type="bicgstab"
#solver_type="ma41"

# Meshes (without .typ2)
#mesh[1]="cart10x10"
#mesh[2]="cart20x20"
#mesh[3]="cart40x40"
#mesh[4]="cart80x80"
#mesh[5]="cart160x160"

#mesh[1]="mesh1_3"
#mesh[2]="mesh1_4"
#mesh[3]="mesh1_5"
#mesh[4]="mesh1_6"
#mesh[5]="mesh1_7"

#mesh[1]="mesh4_1_1"
#mesh[2]="mesh4_1_2"
#mesh[3]="mesh4_1_3"
#mesh[4]="mesh4_1_4"
#mesh[5]="mesh4_1_5"
#mesh[6]="mesh4_1_6"

#mesh[1]="mesh3_1"
#mesh[2]="mesh3_2"
#mesh[3]="mesh3_3"
#mesh[4]="mesh3_4"
#mesh[5]="mesh3_5"

mesh[1]="hexa1_1"
mesh[2]="hexa1_2"
mesh[3]="hexa1_3"
#mesh[4]="hexa1_4"
#mesh[5]="hexa1_5"


