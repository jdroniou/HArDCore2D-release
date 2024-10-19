# Library: ddrcore

FILE(GLOB source_list ddrcore.* platescore.* xgrad.* xcurl.* excurl.* operatorfillers.* xdivdiv.* xrotrot.* xrot.* compute_rotrot_norm.*)
ADD_LIBRARY(ddrcore ${source_list})
TARGET_LINK_LIBRARIES(ddrcore basis mesh quadrature Eigen3::Eigen)
