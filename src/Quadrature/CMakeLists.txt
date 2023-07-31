# Library: quadrature

FILE(GLOB source_list quadraturerule.* triangle_dunavant_rule.* quad1d.* quad2d.* GMpoly_cell.* GMpoly_edge.*)
ADD_LIBRARY(quadrature ${source_list})
TARGET_LINK_LIBRARIES(quadrature Eigen3::Eigen)

