# Library: hybridcore

FILE(GLOB source_list hybridcore.* elementquad.*)
ADD_LIBRARY(hybridcore ${source_list})
TARGET_LINK_LIBRARIES(hybridcore mesh quadrature)


