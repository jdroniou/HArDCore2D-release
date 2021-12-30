# Library: ddrcore

FILE(GLOB source_list ddrcore.* platescore.* xgrad.* xcurl.* excurl.* xdivdiv.*)
ADD_LIBRARY(ddrcore ${source_list})
TARGET_LINK_LIBRARIES(ddrcore basis mesh quadrature)
