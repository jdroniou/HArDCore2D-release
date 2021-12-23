# Library: ddrcore

FILE(GLOB source_list ddrcore.* xgrad.* xcurl.* excurl.*)
ADD_LIBRARY(ddrcore ${source_list})
TARGET_LINK_LIBRARIES(ddrcore basis mesh quadrature)
