# Libraries
ADD_LIBRARY(TestCase OBJECT TestCase/TestCase.cpp)
TARGET_LINK_LIBRARIES(TestCase Eigen3::Eigen)

ADD_LIBRARY(BoundaryConditions BoundaryConditions/BoundaryConditions.cpp)
TARGET_LINK_LIBRARIES(BoundaryConditions Eigen3::Eigen)

ADD_LIBRARY(TestCaseStefanPME OBJECT TestCase/TestCaseNonLinearity.cpp TestCase/TestCaseStefanPME.cpp)
TARGET_LINK_LIBRARIES(TestCaseStefanPME Eigen3::Eigen)

ADD_LIBRARY(TestCaseTransient OBJECT TestCase/TestCaseTransient.cpp)
TARGET_LINK_LIBRARIES(TestCaseTransient Eigen3::Eigen)

ADD_LIBRARY(hho-general OBJECT HHO-general/HHO2D.cpp)
TARGET_LINK_LIBRARIES(hho-general Eigen3::Eigen)

ADD_LIBRARY(DirectedGraph OBJECT MeshCoarsen/DirectedGraph.cpp)
TARGET_LINK_LIBRARIES(DirectedGraph Eigen3::Eigen)

# Mesh Coarsen
ADD_EXECUTABLE(mesh-coarsen $<TARGET_OBJECTS:DirectedGraph> MeshCoarsen/MeshCoarsen.cpp)
TARGET_LINK_LIBRARIES(mesh-coarsen ${Boost_LIBRARIES} Eigen3::Eigen)
TARGET_INCLUDE_DIRECTORIES(mesh-coarsen PUBLIC ${EIGEN3_INCLUDE_DIR})

# HHO Diffusion
ADD_EXECUTABLE(hho-diffusion $<TARGET_OBJECTS:TestCase> $<TARGET_OBJECTS:BoundaryConditions> HHO-diffusion/HHO_Diffusion.cpp)
TARGET_LINK_LIBRARIES(hho-diffusion hybridcore plot ${Boost_LIBRARIES} Eigen3::Eigen)

# HHO Diffusion Advection Reaction 
ADD_EXECUTABLE(hho-diff-advec-reac $<TARGET_OBJECTS:hho-general> $<TARGET_OBJECTS:TestCase>  $<TARGET_OBJECTS:BoundaryConditions> HHO-diff-advec-reac/HHO_DiffAdvecReac.cpp)
TARGET_LINK_LIBRARIES(hho-diff-advec-reac hybridcore plot ${Boost_LIBRARIES} Eigen3::Eigen)

# HHO Full Gradient Diffusion
ADD_EXECUTABLE(hho-fullgradientdiff HHO-full-gradient-diff/hho-fullgradientdiff.cpp)
TARGET_LINK_LIBRARIES(hho-fullgradientdiff hhocore plot BoundaryConditions ${EXTERNAL_LIBRARIES})
TARGET_INCLUDE_DIRECTORIES(hho-fullgradientdiff PUBLIC ${EXTERNAL_INCLUDES})

# HHO Locally variable diffusion
ADD_EXECUTABLE(hho-locvardiff $<TARGET_OBJECTS:TestCase> $<TARGET_OBJECTS:BoundaryConditions> HHO-locvardiff/HHO_LocVarDiff.cpp)
TARGET_LINK_LIBRARIES(hho-locvardiff hybridcore plot ${Boost_LIBRARIES} Eigen3::Eigen)

# Locally Enriched Polytopal Non-Conforming for diffusion
ADD_EXECUTABLE(lepnc-diffusion $<TARGET_OBJECTS:TestCase> $<TARGET_OBJECTS:BoundaryConditions> LEPNC/LEPNC_diffusion.cpp)
TARGET_LINK_LIBRARIES(lepnc-diffusion hybridcore plot ${Boost_LIBRARIES} Eigen3::Eigen)

# Locally Enriched Polytopal Non-Conforming for Stefan/PME
ADD_EXECUTABLE(lepnc-stefanpme $<TARGET_OBJECTS:TestCase> $<TARGET_OBJECTS:TestCaseStefanPME> $<TARGET_OBJECTS:BoundaryConditions> LEPNC/LEPNC_StefanPME.cpp)
TARGET_LINK_LIBRARIES(lepnc-stefanpme hybridcore plot ${Boost_LIBRARIES} Eigen3::Eigen)

# Locally Enriched Polytopal Non-Conforming for transient Stefan/PME
ADD_EXECUTABLE(lepnc-stefanpme-transient $<TARGET_OBJECTS:TestCase> $<TARGET_OBJECTS:TestCaseStefanPME> $<TARGET_OBJECTS:BoundaryConditions> $<TARGET_OBJECTS:TestCaseTransient> LEPNC/LEPNC_StefanPME_transient.cpp)
TARGET_LINK_LIBRARIES(lepnc-stefanpme-transient hybridcore plot ${Boost_LIBRARIES} Eigen3::Eigen)

# HMM for transient Stefan/PME
ADD_EXECUTABLE(hmm-stefanpme-transient $<TARGET_OBJECTS:TestCase> $<TARGET_OBJECTS:TestCaseStefanPME> $<TARGET_OBJECTS:BoundaryConditions> $<TARGET_OBJECTS:TestCaseTransient> HMM/HMM_StefanPME_transient.cpp)
TARGET_LINK_LIBRARIES(hmm-stefanpme-transient hybridcore plot ${Boost_LIBRARIES} Eigen3::Eigen)

# DDR Reissner-Mindlin
ADD_EXECUTABLE(ddr-rmplate $<TARGET_OBJECTS:BoundaryConditions> DDR-rmplate/ddr-rmplate.cpp)
TARGET_LINK_LIBRARIES(ddr-rmplate ddrcore plot ${EXTERNAL_LIBRARIES} Eigen3::Eigen)
TARGET_INCLUDE_DIRECTORIES(ddr-rmplate PUBLIC ${EXTERNAL_INCLUDES})

# DDR Kirchhoff-Love
ADD_EXECUTABLE(ddr-klplate $<TARGET_OBJECTS:BoundaryConditions> DDR-klplate/ddr-klplate.cpp)
TARGET_LINK_LIBRARIES(ddr-klplate ddrcore plot ${EXTERNAL_LIBRARIES})
TARGET_INCLUDE_DIRECTORIES(ddr-klplate PUBLIC ${EXTERNAL_INCLUDES})

# DDR quad-rot
ADD_EXECUTABLE(ddr-quadrot $<TARGET_OBJECTS:BoundaryConditions> DDR-quadrot/ddr-quadrot.cpp)
TARGET_LINK_LIBRARIES(ddr-quadrot ddrcore plot ${EXTERNAL_LIBRARIES})
TARGET_INCLUDE_DIRECTORIES(ddr-quadrot PUBLIC ${EXTERNAL_INCLUDES})
