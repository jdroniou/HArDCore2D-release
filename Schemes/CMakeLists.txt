# TestCase
ADD_LIBRARY(TestCase OBJECT TestCase/TestCase.cpp)

# HHO Diffusion
ADD_EXECUTABLE(hho-diffusion $<TARGET_OBJECTS:TestCase> HHO-diffusion/HHO_Diffusion.cpp)
TARGET_LINK_LIBRARIES(hho-diffusion mesh quadrature hybridcore plot ${Boost_LIBRARIES})
TARGET_INCLUDE_DIRECTORIES(hho-diffusion PUBLIC ${EIGEN3_INCLUDE_DIR})

# HHO Locally variable diffusion
ADD_EXECUTABLE(hho-locvardiff $<TARGET_OBJECTS:TestCase> HHO-locvardiff/HHO_LocVarDiff.cpp)
TARGET_LINK_LIBRARIES(hho-locvardiff mesh quadrature hybridcore plot ${Boost_LIBRARIES})
TARGET_INCLUDE_DIRECTORIES(hho-locvardiff PUBLIC ${EIGEN3_INCLUDE_DIR})

