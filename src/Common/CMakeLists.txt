# Library: ddrcore

FILE(GLOB source_list basis.* dofspace.*)
ADD_LIBRARY(basis ${source_list})
TARGET_LINK_LIBRARIES(basis mesh quadrature)
