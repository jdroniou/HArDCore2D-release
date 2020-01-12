# Library: mesh

FILE(GLOB source_list mesh.* mesh_builder.* import_mesh.* cell.* edge.* vertex.*)
ADD_LIBRARY(mesh ${source_list})

