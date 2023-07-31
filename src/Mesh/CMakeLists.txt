# Library: mesh

FILE(GLOB source_list mesh.*)
ADD_LIBRARY(mesh ${source_list})
TARGET_LINK_LIBRARIES(mesh Eigen3::Eigen)

