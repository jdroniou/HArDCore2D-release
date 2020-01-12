# Library: quadrature

FILE(GLOB source_list quadraturerule.* triangle_dunavant_rule.* quad1d.* quad2d.*)
ADD_LIBRARY(quadrature ${source_list})

