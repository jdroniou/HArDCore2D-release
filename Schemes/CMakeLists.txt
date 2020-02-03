# TestCase
ADD_LIBRARY(TestCase OBJECT TestCase/TestCase.cpp TestCase/BoundaryConditions.cpp)
ADD_LIBRARY(TestCaseStefanPME OBJECT TestCase/TestCaseNonLinearity.cpp TestCase/TestCaseStefanPME.cpp)

# HHO Diffusion
ADD_EXECUTABLE(hho-diffusion $<TARGET_OBJECTS:TestCase> HHO-diffusion/HHO_Diffusion.cpp)
TARGET_LINK_LIBRARIES(hho-diffusion mesh quadrature hybridcore plot ${Boost_LIBRARIES})
TARGET_INCLUDE_DIRECTORIES(hho-diffusion PUBLIC ${EIGEN3_INCLUDE_DIR})

# HHO Locally variable diffusion
ADD_EXECUTABLE(hho-locvardiff $<TARGET_OBJECTS:TestCase> HHO-locvardiff/HHO_LocVarDiff.cpp)
TARGET_LINK_LIBRARIES(hho-locvardiff mesh quadrature hybridcore plot ${Boost_LIBRARIES})
TARGET_INCLUDE_DIRECTORIES(hho-locvardiff PUBLIC ${EIGEN3_INCLUDE_DIR})

# Bubble Polytopal Non-conforming for Stefan/PME
ADD_EXECUTABLE(bpnc-stefanpme $<TARGET_OBJECTS:TestCase> $<TARGET_OBJECTS:TestCaseStefanPME> BPNC/BPNC_StefanPME.cpp)
TARGET_LINK_LIBRARIES(bpnc-stefanpme mesh quadrature hybridcore plot ${Boost_LIBRARIES})
TARGET_INCLUDE_DIRECTORIES(bpnc-stefanpme PUBLIC ${EIGEN3_INCLUDE_DIR})

