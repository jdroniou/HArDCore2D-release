# Libraries
ADD_LIBRARY(BoundaryConditions BoundaryConditions/BoundaryConditions.cpp)
ADD_LIBRARY(TestCase OBJECT TestCase/TestCase.cpp)
ADD_LIBRARY(TestCaseStefanPME OBJECT TestCase/TestCaseNonLinearity.cpp TestCase/TestCaseStefanPME.cpp)
ADD_LIBRARY(TestCaseTransient OBJECT TestCase/TestCaseTransient.cpp)
ADD_LIBRARY(hho-general OBJECT HHO-general/HHO2D.cpp)
ADD_LIBRARY(DirectedGraph OBJECT MeshCoarsen/DirectedGraph.cpp)

# Mesh Coarsen
ADD_EXECUTABLE(mesh-coarsen $<TARGET_OBJECTS:DirectedGraph> MeshCoarsen/MeshCoarsen.cpp)
TARGET_LINK_LIBRARIES(mesh-coarsen ${Boost_LIBRARIES})
TARGET_INCLUDE_DIRECTORIES(mesh-coarsen PUBLIC ${EIGEN3_INCLUDE_DIR})

# HHO Diffusion
ADD_EXECUTABLE(hho-diffusion $<TARGET_OBJECTS:TestCase> $<TARGET_OBJECTS:BoundaryConditions> HHO-diffusion/HHO_Diffusion.cpp)
TARGET_LINK_LIBRARIES(hho-diffusion mesh quadrature hybridcore plot ${Boost_LIBRARIES})
TARGET_INCLUDE_DIRECTORIES(hho-diffusion PUBLIC ${EIGEN3_INCLUDE_DIR})

# HHO Diffusion Advection Reaction 
ADD_EXECUTABLE(hho-diff-advec-reac $<TARGET_OBJECTS:hho-general> $<TARGET_OBJECTS:TestCase>  $<TARGET_OBJECTS:BoundaryConditions> HHO-diff-advec-reac/HHO_DiffAdvecReac.cpp)
TARGET_LINK_LIBRARIES(hho-diff-advec-reac mesh quadrature hybridcore plot ${Boost_LIBRARIES})
TARGET_INCLUDE_DIRECTORIES(hho-diff-advec-reac PUBLIC ${EIGEN3_INCLUDE_DIR})

# HHO Locally variable diffusion
ADD_EXECUTABLE(hho-locvardiff $<TARGET_OBJECTS:TestCase> $<TARGET_OBJECTS:BoundaryConditions> HHO-locvardiff/HHO_LocVarDiff.cpp)
TARGET_LINK_LIBRARIES(hho-locvardiff mesh quadrature hybridcore plot ${Boost_LIBRARIES})
TARGET_INCLUDE_DIRECTORIES(hho-locvardiff PUBLIC ${EIGEN3_INCLUDE_DIR})

# Locally Enriched Polytopal Non-Conforming for diffusion
ADD_EXECUTABLE(lepnc-diffusion $<TARGET_OBJECTS:TestCase> $<TARGET_OBJECTS:BoundaryConditions> LEPNC/LEPNC_diffusion.cpp)
TARGET_LINK_LIBRARIES(lepnc-diffusion mesh quadrature hybridcore plot ${Boost_LIBRARIES})
TARGET_INCLUDE_DIRECTORIES(lepnc-diffusion PUBLIC ${EIGEN3_INCLUDE_DIR})

# Locally Enriched Polytopal Non-Conforming for Stefan/PME
ADD_EXECUTABLE(lepnc-stefanpme $<TARGET_OBJECTS:TestCase> $<TARGET_OBJECTS:TestCaseStefanPME> $<TARGET_OBJECTS:BoundaryConditions> LEPNC/LEPNC_StefanPME.cpp)
TARGET_LINK_LIBRARIES(lepnc-stefanpme mesh quadrature hybridcore plot ${Boost_LIBRARIES})
TARGET_INCLUDE_DIRECTORIES(lepnc-stefanpme PUBLIC ${EIGEN3_INCLUDE_DIR})

# Locally Enriched Polytopal Non-Conforming for transient Stefan/PME
ADD_EXECUTABLE(lepnc-stefanpme-transient $<TARGET_OBJECTS:TestCase> $<TARGET_OBJECTS:TestCaseStefanPME> $<TARGET_OBJECTS:BoundaryConditions> $<TARGET_OBJECTS:TestCaseTransient> LEPNC/LEPNC_StefanPME_transient.cpp)
TARGET_LINK_LIBRARIES(lepnc-stefanpme-transient mesh quadrature hybridcore plot ${Boost_LIBRARIES})
TARGET_INCLUDE_DIRECTORIES(lepnc-stefanpme-transient PUBLIC ${EIGEN3_INCLUDE_DIR})

# HMM for transient Stefan/PME
ADD_EXECUTABLE(hmm-stefanpme-transient $<TARGET_OBJECTS:TestCase> $<TARGET_OBJECTS:TestCaseStefanPME> $<TARGET_OBJECTS:BoundaryConditions> $<TARGET_OBJECTS:TestCaseTransient> HMM/HMM_StefanPME_transient.cpp)
TARGET_LINK_LIBRARIES(hmm-stefanpme-transient mesh quadrature hybridcore plot ${Boost_LIBRARIES})
TARGET_INCLUDE_DIRECTORIES(hmm-stefanpme-transient PUBLIC ${EIGEN3_INCLUDE_DIR})

# DDR Reissner-Mindlin
ADD_EXECUTABLE(ddr-rmplate $<TARGET_OBJECTS:BoundaryConditions> DDR-rmplate/ddr-rmplate.cpp)
TARGET_LINK_LIBRARIES(ddr-rmplate mesh quadrature ddrcore plot ${EXTERNAL_LIBRARIES})
TARGET_INCLUDE_DIRECTORIES(ddr-rmplate PUBLIC ${EXTERNAL_INCLUDES})


