# Library: hhocore

FILE(GLOB source_list hhospace.*)
ADD_LIBRARY(hhocore ${source_list})
TARGET_LINK_LIBRARIES(hhocore basis mesh quadrature Eigen3::Eigen)
