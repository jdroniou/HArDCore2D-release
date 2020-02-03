#!/bin/bash
#

# Edge and cell degrees
k=2
l=1

# Use threads
use_threads="true"

# Export matrix
export_matrix="true"

# Boundary conditions (D, N, Mx)
bc="M0"

# Test case
tcsol=2
tcdiff=2

# Solver
solver_type="bicgstab"
#solver_type="ma41"

# Meshes (without .typ2)
#mesh[1]="cart10x10"
#mesh[2]="cart20x20"
#mesh[3]="cart40x40"
#mesh[4]="cart80x80"

#mesh[1]="mesh1_2"
#mesh[2]="mesh1_3"
#mesh[3]="mesh1_4"
#mesh[4]="mesh1_5"

mesh[1]="hexa1_1"
mesh[2]="hexa1_2"
mesh[3]="hexa1_3"
#mesh[4]="hexa1_4"
#mesh[5]="hexa1_5"

#mesh[1]="mesh3_2"
#mesh[2]="mesh3_3"
#mesh[3]="mesh3_4"
#mesh[4]="mesh3_5"

#mesh[1]="mesh4_1_1"
#mesh[2]="mesh4_1_2"
#mesh[3]="mesh4_1_3"
#mesh[4]="mesh4_1_4"
#mesh[5]="mesh4_1_5"
#mesh[6]="mesh4_1_6"


###

#mesh[1]="anisotropic/cart25_a1"
#mesh[2]="anisotropic/cart50_a2"
#mesh[3]="anisotropic/cart100_a4"
#mesh[4]="anisotropic/cart200_a8"

#mesh[1]="anisotropic/cart50_a10"
#mesh[2]="anisotropic/cart200_a20"
#mesh[3]="anisotropic/cart400_a50"

#mesh[1]="anisotropic/tri20x20"
#mesh[2]="anisotropic/tri40x80"
#mesh[3]="anisotropic/tri80x320"
#mesh[4]="anisotropic/tri100x500"

#mesh[1]="anisotropic/hexa20x20"
#mesh[2]="anisotropic/hexa40x80"
#mesh[3]="anisotropic/hexa80x320"
#mesh[4]="anisotropic/hexa120x720"

#mesh[1]="anisotropic/cart160x160"
#mesh[2]="anisotropic/cart80x320"
#mesh[3]="anisotropic/cart40x640"
#mesh[4]="anisotropic/cart20x1280"

#mesh[1]="anisotropic/two_quads120x120"
#mesh[2]="anisotropic/two_quads60x240"
#mesh[3]="anisotropic/two_quads30x480"
#mesh[4]="anisotropic/two_quads15x960"

#mesh[1]="anisotropic/tri_aniso60x60"
#mesh[2]="anisotropic/tri_aniso30x120"
#mesh[3]="anisotropic/tri_aniso15x240"

#mesh[1]="anisotropic/quads120x120"
#mesh[2]="anisotropic/quads60x240"
#mesh[3]="anisotropic/quads30x480"
#mesh[4]="anisotropic/quads15x960"

#mesh[1]="anisotropic/hexa_straight10x10"
#mesh[2]="anisotropic/hexa_straight20x40"
#mesh[3]="anisotropic/hexa_straight40x160"
#mesh[4]="anisotropic/hexa_straight80x640"

#mesh[1]="anisotropic/hexa_straight10x10"
#mesh[2]="anisotropic/hexa_straight20x20"
#mesh[3]="anisotropic/hexa_straight40x40"
#mesh[4]="anisotropic/hexa_straight80x80"
#mesh[5]="anisotropic/hexa_straight160x160"

