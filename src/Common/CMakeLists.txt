# Library: basis

FILE(GLOB source_list basis.* localdofspace.* globaldofspace.*)
ADD_LIBRARY(basis ${source_list})
TARGET_LINK_LIBRARIES(basis mesh quadrature)
