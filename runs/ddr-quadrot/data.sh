#!/bin/bash
#

# Degree
degrees[1]=0
degrees[2]=1
degrees[3]=2

# Solution
solution=4

# Stabilization parameter
stab_par=1

# Use threads
use_threads="true"

# Meshes
mesh_family=cart

case ${mesh_family} in
    cart)
        mesh[1]="cart10x10"
        mesh[2]="cart20x20"
        mesh[3]="cart40x40"
        #  mesh[4]="cart80x80"
        ;;    
    hexa)
        mesh[1]="hexa1_1"
        mesh[2]="hexa1_2"
        mesh[3]="hexa1_3"
        mesh[4]="hexa1_4"
        #  mesh[5]="hexa1_5"
        #  mesh[6]="hexa1_6"
        #  mesh[7]="hexa1_7"
        ;;
    locref)
        mesh[1]="mesh3_2"
        mesh[2]="mesh3_3"
        mesh[3]="mesh3_4"
        mesh[4]="mesh3_5"
        ;;
    tri)
        mesh[1]="mesh1_1"
        mesh[2]="mesh1_2"
        mesh[3]="mesh1_3"
        mesh[4]="mesh1_4"
        mesh[5]="mesh1_5"
        # mesh[6]="mesh1_6"
        # mesh[7]="mesh1_7"
        ;;    
    tri2)
        mesh[1]="tri1_1"
        mesh[2]="tri1_2"
        mesh[3]="tri1_3"
        mesh[4]="tri1_4"
        mesh[5]="tri1_5"
        ;;    
    cart_refined_boundary)
        mesh[1]="cart_refined_boundary1"
        mesh[2]="cart_refined_boundary2"
        mesh[3]="cart_refined_boundary3"
        #    mesh[4]="cart_refined_boundary4"
        ;;
    tri2_refined_boundary)
        mesh[1]="tri2_refined_boundary1"
        mesh[2]="tri2_refined_boundary2"
        mesh[3]="tri2_refined_boundary3"
        #    mesh[4]="tri2_refined_boundary4"
        ;;    
    tri_agglo)
        mesh[1]="agglomerated/mesh1_1.coarse.0"
        mesh[2]="agglomerated/mesh1_2.coarse.2"
        mesh[3]="agglomerated/mesh1_3.coarse.3"
esac
