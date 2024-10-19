//
// Solver for a quad-rot problem
// Author: Daniele Di Pietro (daniele.di-pietro@umontpellier.fr)
//
#include <fstream>
#include <iomanip>
#include <thread>

#include "ddr-quadrot.hpp"
#include "compute_rotrot_norm.hpp"

#include <parallel_for.hpp>
#include "vtu_writer.hpp"

#include <boost/program_options.hpp>
#include <boost/timer/timer.hpp>

#ifdef WITH_UMFPACK
#include <Eigen/UmfPackSupport>
#endif

#ifdef WITH_MKL
#include <Eigen/PardisoSupport>
#include <mkl.h>
#endif

#define FORMAT(W)                                                       \
  std::setiosflags(std::ios_base::left) << std::setw(W) << std::setfill(' ')

#include <BoundaryConditions/BChandlers.hpp>

using namespace HArDCore2D;

//------------------------------------------------------------------------------
// Mesh filenames
//------------------------------------------------------------------------------

const std::string mesh_dir = "../../typ2_meshes/";
std::string default_mesh = mesh_dir + "hexa1_1.typ2";

// Max number of cells to plot a graph
constexpr size_t max_nb_cells_for_plot = 15000;

//------------------------------------------------------------------------------

int main(int argc, const char* argv[])
{
  // Program options
  boost::program_options::options_description desc("Allowed options");
  desc.add_options()
    ("help,h", "Display this help message")
    ("mesh,m", boost::program_options::value<std::string>(), "Set the mesh")
    ("degree,k", boost::program_options::value<size_t>()->default_value(0), "The polynomial degree of the sequence")
    ("sequential-execution,q", "Use thread-based parallelism")
    ("solution,s", boost::program_options::value<int>()->default_value(1), "Select the solution")
    ("plot", boost::program_options::value<std::string>()->default_value("displacement"), "Save plot of the solution to the given filename")
    ("export-matrix,e", "Export matrix to Matrix Market format")
    ("iterative-solver,i", "Use iterative linear solver")
    ("stabilization-parameter,x", boost::program_options::value<double>()->default_value(1.), "Set the stabilization parameter")
    ("test-reconstructions,r", "Test the convergence of reconstructions");

  boost::program_options::variables_map vm;
  boost::program_options::store(boost::program_options::parse_command_line(argc, argv, desc), vm);
  boost::program_options::notify(vm);

  // Display the help options
  if (vm.count("help")) {
    std::cout << desc << std::endl;
    return 0;
  }

  // Select the mesh
  std::string mesh_file = (vm.count("mesh") ? vm["mesh"].as<std::string>() : default_mesh);
  std::cout << FORMAT(50) << "[main] Mesh file" << mesh_file << std::endl;
  
  // Select the degree 
       size_t K = vm["degree"].as<size_t>();
  std::cout << FORMAT(50) << "[main] Degree" << K << std::endl;

  // Initialize the exact solution
  int solution = (vm.count("solution") ? vm["solution"].as<int>() : 0);
  QuadRot::PotentialType u;
  QuadRot::RotorType rot_u;
  QuadRot::RotRotType rotrot_u;
  QuadRot::LagrangeMultiplierType p;
  QuadRot::ForcingTermType f;

  switch (solution) {
  case 0:
    std::cout << "[main] Trigonometric solution" << std::endl;    
    f = trigonometric_f;
    u = trigonometric_u;
    rot_u = trigonometric_rot_u;
    rotrot_u = trigonometric_rotrot_u;
    p = trigonometric_p;
    break;
    
  case 1:
    std::cout << "[main] Linear solution" << std::endl;    
    f = linear_f;
    u = linear_u;
    rot_u = linear_rot_u;
    rotrot_u = linear_rotrot_u;
    p = linear_p;
    break;

  case 2:
    std::cout << "[main] Quadratic solution" << std::endl;
    f = quadratic_f;
    u = quadratic_u;
    rot_u = quadratic_rot_u;
    rotrot_u = quadratic_rotrot_u;
    p = quadratic_p;
    break;

  case 3:
    std::cout << "[main] Cubic solution" << std::endl;
    f = cubic_f;
    u = cubic_u;
    rot_u = cubic_rot_u;
    rotrot_u = cubic_rotrot_u;
    p = cubic_p;
    break;
    
  case 4:
    std::cout << "[main] Quartic solution" << std::endl;
    f = quartic_f;
    u = quartic_u;
    rot_u = quartic_rot_u;
    rotrot_u = quartic_rotrot_u;
    p = quartic_p;
    break;

  default:
    std::cerr << "[main] ERROR: Unknown exact solution" << std::endl;
    exit(1);  
  }
  
  // Build the mesh and re-order vertices and edges to put the boundary ones first
  MeshBuilder builder = MeshBuilder(mesh_file);
  std::unique_ptr<Mesh> mesh_ptr = builder.build_the_mesh();
  std::cout << FORMAT(50) << "[main] Mesh size" << mesh_ptr->h_max() << std::endl;

  // Handle the boundary conditions
  BoundaryConditions bc_u("D", *mesh_ptr.get());
  bc_u.reorder_vertices("start");
  bc_u.reorder_edges("start");
  BoundaryConditions bc_p("D", *mesh_ptr.get());
  bc_p.reorder_vertices("start");
  bc_p.reorder_edges("start");
  
  boost::timer::cpu_timer timer;
  
  // Create DDR core
  timer.start();
  bool use_threads = vm.count("sequential-execution") ? false : true;
  std::cout << "[main] " << (use_threads ? "Parallel execution" : "Sequential execution") << std:: endl;
  DDRCore ddr_core(*mesh_ptr, K, use_threads);
  timer.stop();
  double t_wall_ddrcore = double(timer.elapsed().wall) * pow(10, -9);
  double t_proc_ddrcore = double(timer.elapsed().user + timer.elapsed().system) * pow(10, -9);
  std::cout << "[main] Time DDRCore (wall/proc) " << t_wall_ddrcore << "/" << t_proc_ddrcore << std::endl;

  // Assemble the problem
  timer.start();
  
  QuadRot quadrot(ddr_core, bc_u, bc_p, use_threads);
  if(vm.count("stabilization-parameter")) {
    quadrot.stabilizationParameter() = vm["stabilization-parameter"].as<double>();    
  }
  
  std::cout << FORMAT(50) << "[main] Stabilization parameter" << quadrot.stabilizationParameter() << std::endl;
  std::cout << FORMAT(50) << "[main] Dimension of the spaces " << quadrot.dimensionSpace() << std::endl;
  std::cout << FORMAT(50) << "[main] Dimension of the linear system" << quadrot.dimensionLinearSystem() << std::endl;

  Eigen::VectorXd uI = quadrot.assembleLinearSystem(f, u, rot_u, p);
  timer.stop();
  double t_wall_model = double(timer.elapsed().wall) * pow(10, -9);
  double t_proc_model = double(timer.elapsed().user + timer.elapsed().system) * pow(10, -9);
  std::cout << FORMAT(50) << "[main] Wall time model " << t_wall_model << std::endl;
  std::cout << FORMAT(50) << "[main] Proc time model " << t_proc_model << std::endl;

  // Interpolate the exact solution
  std::cout << FORMAT(50) << "[main] Interpolating the exact solution" << std::endl;
  Eigen::VectorXd rot_uI = quadrot.xRot().interpolate(rot_u);

  // Export matrix if requested  
  if (vm.count("export-matrix")) {
    std::cout << "[main] Exporting matrix to Matrix Market format" << std::endl;
    saveMarket(quadrot.systemMatrix(), "A_quadrot.mtx");
    saveMarket(quadrot.systemVector(), "b_quadrot.mtx");
  }
  
  // Solve the problem
  timer.start();
  Eigen::VectorXd uh_int;
  if (vm.count("iterative-solver")) {
    std::cout << "[main] Solving the linear system using BiCGSTAB" << std::endl;
    
    Eigen::BiCGSTAB<QuadRot::SystemMatrixType, Eigen::IncompleteLUT<double> > solver;
    // solver.preconditioner().setFillfactor(2);
    solver.compute(quadrot.systemMatrix());
    if (solver.info() != Eigen::Success) {
      std::cerr << "[main] ERROR: Could not factorize matrix" << std::endl;
      exit(1);
    }
    uh_int = solver.solve(quadrot.systemVector());
    if (solver.info() != Eigen::Success) {
      std::cerr << "[main] ERROR: Could not solve direct system" << std::endl;
      exit(1);
    }
  } else { 
#ifdef WITH_MKL
    std::cout << "[main] Solving the linear system using Pardiso" << std::endl;    
    unsigned nb_threads_hint = std::thread::hardware_concurrency();
    mkl_set_dynamic(0);
    mkl_set_num_threads(nb_threads_hint);
    Eigen::PardisoLU<QuadRot::SystemMatrixType> solver;
#elif WITH_UMFPACK
    std::cout << "[main] Solving the linear system using Umfpack" << std::endl;    
    Eigen::UmfPackLU<QuadRot::SystemMatrixType> solver;
#else
    std::cout << "[main] Solving the linear system using direct solver" << std::endl;    
    Eigen::SparseLU<QuadRot::SystemMatrixType> solver;
#endif
    solver.compute(quadrot.systemMatrix());
    if (solver.info() != Eigen::Success) {
      std::cerr << "[main] ERROR: Could not factorize matrix" << std::endl;
    }
    uh_int = solver.solve(quadrot.systemVector());
    if (solver.info() != Eigen::Success) {
      std::cerr << "[main] ERROR: Could not solve linear system" << std::endl;
    }
  }
  // Norm of residual
  std::cout << FORMAT(50) << "[main] Norm of residual after solving " << (quadrot.systemMatrix()*uh_int-quadrot.systemVector()).norm() << std::endl;

  // Reconstruct the exact solution from the boundary datum lifting uI
  // and uh_int
  Eigen::VectorXd uh = Eigen::VectorXd::Zero(quadrot.dimensionSpace());
  for (size_t i = 0; i < quadrot.dimensionSpace(); i++) {
    int unknown_i = quadrot.dofToUnknown(i);
    if (unknown_i >= 0) {
      uh(i) = uh_int(unknown_i);
    } else {
      uh(i) = uI(i);
    }
  } // for i

  // Compute the errors
  Eigen::VectorXd errh = uh - uI;
  const auto & err_u  = errh.head(quadrot.xRotRot().dimension());
  const auto & err_p  = errh.tail(quadrot.xGrad().dimension());

  double err_u_l2     = quadrot.xRotRot().computeL2Norm(err_u);
  double err_u_rotrot = compute_rotrot_norm(err_u, quadrot.xRot(), quadrot.xRotRot());
  double err_p_l2     = quadrot.xGrad().computeL2Norm(err_p);
  double err_p_grad   = quadrot.xRotRot().computeGradientL2Norm(err_p, &quadrot.xGrad());

  //------------------------------------------------------------------------------

  std::cout << FORMAT(50) << "[main] degree" << quadrot.xRotRot().degree() << std::endl;

  std::cout << FORMAT(50) << "[main] meshsize" << mesh_ptr->h_max() << std::endl;
  std::cout << FORMAT(50) << "[main] nb_cells" << mesh_ptr->n_cells() << std::endl;
  std::cout << FORMAT(50) << "[main] nb_edges" << mesh_ptr->n_edges() << std::endl;
  std::cout << FORMAT(50) << "[main] nb_vertices" << mesh_ptr->n_vertices() << std::endl;
  std::cout << FORMAT(50) << "[main] dim_spaces " << quadrot.dimensionSpace() << std::endl;
  std::cout << FORMAT(50) << "[main] dim_linsys" << quadrot.dimensionLinearSystem() << std::endl;
  
 
  std::cout << FORMAT(50) << "[main] err_u_l2" << std::scientific << err_u_l2 << std::endl;
  std::cout << FORMAT(50) << "[main] err_u_rotrot" << std::scientific << err_u_rotrot << std::endl;
  std::cout << FORMAT(50) << "[main] err_p_l2" << std::scientific << err_p_l2 << std::endl;
  std::cout << FORMAT(50) << "[main] err_p_grad" << std::scientific << err_p_grad << std::endl;

  
  // Export data and results
  std::ofstream out("results.txt");
  out<< "Degree: " << quadrot.xRotRot().degree() << std::endl;

  out<< "MeshSize: " << mesh_ptr->h_max() << std::endl;
  out<< "NbCells: " << mesh_ptr->n_cells() << std::endl;
  out<< "NbEdges: " << mesh_ptr->n_edges() << std::endl;
  out<< "NbVertices: " << mesh_ptr->n_vertices() << std::endl;

  out << "DimSpaces: " << quadrot.dimensionSpace() << std::endl;
  out << "DimLinSys: " << quadrot.dimensionLinearSystem() << std::endl;  
 
  out<< "ErrUL2: " << std::scientific << err_u_l2 << std::endl;
  out<< "ErrURotRot: " << std::scientific << err_u_rotrot << std::endl;
  out<< "ErrPL2: " << std::scientific << err_p_l2 << std::endl;
  out<< "ErrPGrad: " << std::scientific << err_p_grad << std::endl;
  out.close();

  //------------------------------------------------------------------------------
  // Test reconstructions
  //------------------------------------------------------------------------------

  if (vm.count("test-reconstructions")) {
    double err_pot_uI = 0.;        // Potential reconstruction in XRotRot
    double err_rot_uI = 0.;        // Scalar rotor of the interpolate of u
    double err_commut_rot_uI = 0.; // To check the commutativity RT IT = lprojkT
    double err_rotrot_uI = 0.;     // Vector rotor of the scalar rotor of the interpolate of u
    double err_rot_rot_uI = 0.;    // Vector rotor of the interpolate of rot(u)
    double err_pot_rot_uI = 0.;    // Potential of the interpolate of rot(u)
    double err_rotor_basis = 0.;   // Check that rotor basis behaves as expected
  
    for (size_t iT = 0; iT < mesh_ptr->n_cells(); iT++) {
      const Cell & T = *mesh_ptr->cell(iT);
   
      Eigen::VectorXd uIT = quadrot.xRotRot().restrictCell(iT, uI.head(quadrot.xRotRot().dimension()));
      Eigen::VectorXd rot_uIT = quadrot.xRot().restrictCell(iT, rot_uI);
       
      QuadratureRule quad_2k_T = generate_quadrature_rule(T, 2 * quadrot.xRotRot().degree());
      QuadratureRule quad_2kp2_T = generate_quadrature_rule(T, 2 * quadrot.xRotRot().degree() + 2);
  
      // Compute the error on the rotor
      size_t dim_Pk = quadrot.xRotRot().cellBases(iT).Polyk->dimension();
      size_t dim_Pk2 = quadrot.xRot().cellBases(iT).Polyk2->dimension();
      
      Eigen::VectorXd PT_uI = quadrot.xRotRot().cellOperators(iT).potential * uIT;
      Eigen::VectorXd RT_uI = quadrot.xRotRot().cellRotor(iT) * uIT;
      Eigen::VectorXd diff_RT_uI_lproj_rot_u = RT_uI - rot_uIT.segment(quadrot.xRot().localOffset(T), dim_Pk);
 
      auto basis_Pk_T_quad = evaluate_quad<Function>::compute(*quadrot.xRotRot().cellBases(iT).Polyk, quad_2k_T);
      auto basis_Pk2_T_quad = evaluate_quad<Function>::compute(*quadrot.xRot().cellBases(iT).Polyk2, quad_2k_T);
      
      for (size_t iqn = 0; iqn < quad_2k_T.size(); iqn++) {
        VectorRd PT_uI_iqn = VectorRd::Zero();
        for (size_t i = 0; i < dim_Pk2; i++) {
          PT_uI_iqn += basis_Pk2_T_quad[i][iqn] * PT_uI(i);
        } // for i
        err_pot_uI += quad_2k_T[iqn].w * (PT_uI_iqn - u(quad_2k_T[iqn].vector())).squaredNorm();
        
        double RT_uI_iqn = 0.;
        double diff_RT_uI_lproj_rot_u_iqn = 0.;
        for (size_t i = 0; i < dim_Pk; i++) {
          RT_uI_iqn += basis_Pk_T_quad[i][iqn] * RT_uI(i);
          diff_RT_uI_lproj_rot_u_iqn += basis_Pk_T_quad[i][iqn] * diff_RT_uI_lproj_rot_u(i);
        } // for i
        double rot_u_iqn = rot_u(quad_2k_T[iqn].vector());

        err_rot_uI += quad_2k_T[iqn].w * std::pow(RT_uI_iqn - rot_u_iqn, 2);

        err_commut_rot_uI += quad_2k_T[iqn].w * std::pow(diff_RT_uI_lproj_rot_u_iqn, 2);
      } // for iqn
    
      // Compute the error on the potential reconstruction of rot(u)
      Eigen::VectorXd PT_rot_uI = quadrot.xRot().cellOperators(iT).potential * rot_uIT;
      // std::cout << PT_rot_uI.transpose() << std::endl;

      size_t dim_Pkpo = quadrot.xRot().cellBases(iT).Polykpo->dimension();
      auto basis_Pkpo_T_quad = evaluate_quad<Function>::compute(*quadrot.xRotRot().cellBases(iT).Polykpo, quad_2kp2_T);

      for (size_t iqn = 0; iqn < quad_2kp2_T.size(); iqn++) {
        double PT_rot_uI_iqn = 0.;
        for (size_t i = 0; i < dim_Pkpo; i++) {
          PT_rot_uI_iqn += basis_Pkpo_T_quad[i][iqn] * PT_rot_uI(i);
        } // for i
        double rot_u_iqn = rot_u(quad_2kp2_T[iqn].vector());

        err_pot_rot_uI += quad_2kp2_T[iqn].w * std::pow(PT_rot_uI_iqn - rot_u_iqn, 2);
      } // for iqn
    
      // Compute the error on rotrot
      Eigen::VectorXd RRT_uI = quadrot.xRot().cellOperators(iT).rotor * quadrot.xRotRot().cellOperators(iT).rotor * uIT;
      Eigen::VectorXd R_rot_uI = quadrot.xRot().cellOperators(iT).rotor * rot_uIT;
    
      for (size_t iqn = 0; iqn < quad_2k_T.size(); iqn++) {
        VectorRd RRT_uI_iqn = VectorRd::Zero();
        VectorRd R_rot_uI_iqn = VectorRd::Zero();
        for (size_t i = 0; i < dim_Pk2; i++) {
          RRT_uI_iqn   += basis_Pk2_T_quad[i][iqn] * RRT_uI(i);
          R_rot_uI_iqn += basis_Pk2_T_quad[i][iqn] * R_rot_uI(i);
        } // for i
      
        VectorRd rotrot_u_iqn = rotrot_u(quad_2k_T[iqn].vector());
      
        err_rotrot_uI  += quad_2k_T[iqn].w * (RRT_uI_iqn - rotrot_u_iqn).squaredNorm();
        err_rot_rot_uI += quad_2k_T[iqn].w * (R_rot_uI_iqn - rotrot_u_iqn).squaredNorm();
      } // for iqn

      // Test rotor basis
      if (quadrot.xRot().degree() >= 1) {
        RotorBasis<DDRCore::Poly2BasisCellType> basis_rot_Pk2(*quadrot.xRot().cellBases(iT).Polyk2);
        auto basis_rot_Pk2_T_quad = evaluate_quad<Function>::compute(basis_rot_Pk2, quad_2k_T);
        Eigen::VectorXd pik_u_T = l2_projection(u, *quadrot.xRot().cellBases(iT).Polyk2, quad_2k_T, basis_Pk2_T_quad);
      
        for (size_t iqn = 0; iqn < quad_2k_T.size(); iqn++) {
          double rot_pik_u_T_iqn = 0.;
          for (size_t i = 0; i < dim_Pk2; i++) {
            rot_pik_u_T_iqn += basis_rot_Pk2_T_quad[i][iqn] * pik_u_T(i);
          } // for i
          err_rotor_basis += quad_2k_T[iqn].w * pow(rot_pik_u_T_iqn - rot_u(quad_2kp2_T[iqn].vector()), 2);
        } // for iqn
      } // if
      
    } // for iT

    err_pot_uI = std::sqrt(err_pot_uI);
    err_rot_uI = std::sqrt(err_rot_uI);
    err_rotrot_uI = std::sqrt(err_rotrot_uI);
    err_rot_rot_uI = std::sqrt(err_rot_rot_uI);
    err_pot_rot_uI = std::sqrt(err_pot_rot_uI);
    err_rotor_basis = std::sqrt(err_rotor_basis);

    std::cout << FORMAT(50) << "[main] err_pot_uI" << std::scientific << err_pot_uI << std::endl;
    std::cout << FORMAT(50) << "[main] err_rot_uI" << std::scientific << err_rot_uI << std::endl;
    std::cout << FORMAT(50) << "[main] err_commut_rot_uI" << std::scientific << err_commut_rot_uI << std::endl;
    std::cout << FORMAT(50) << "[main] err_rotrot_uI" << std::scientific << err_rotrot_uI << std::endl;
    std::cout << FORMAT(50) << "[main] err_rot_rot_uI" << std::scientific << err_rot_rot_uI << std::endl;
    std::cout << FORMAT(50) << "[main] err_pot_rot_uI" << std::scientific << err_pot_rot_uI << std::endl;
    std::cout << FORMAT(50) << "[main] err_rotor_basis" << std::scientific << err_rotor_basis << std::endl;
  } // if

  std::cout << "[main] Done" << std::endl;
  return 0;
}

//------------------------------------------------------------------------------
// QuadRot
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
// Constructor

QuadRot::QuadRot(
                 const DDRCore & ddrcore,
                 const BoundaryConditions & bc_u,
                 const BoundaryConditions & bc_p,
                 bool use_threads,
                 std::ostream & output
                 )
  : m_ddrcore(ddrcore),
    m_bc_u(bc_u),
    m_bc_p(bc_p),
    m_use_threads(use_threads),
    m_output(output),
    m_xgrad(ddrcore, use_threads),
    m_xrotrot(ddrcore, use_threads),
    m_xrot(ddrcore, use_threads),
    m_unknowns_stride(4)
{
  m_output << "[QuadRot] Initializing" << std::endl;

  // Vertex unknowns for u
  m_unknowns_stride[0].first  = m_bc_u.n_dir_vertices() * m_xrotrot.numLocalDofsVertex();
  m_unknowns_stride[0].second =
    ( m_ddrcore.mesh().n_vertices() - m_bc_u.n_dir_vertices() ) * m_xrotrot.numLocalDofsVertex();
  // Edge and cell unknowns for u
  m_unknowns_stride[1].first  = m_ddrcore.mesh().n_vertices() * m_xrotrot.numLocalDofsVertex()
    + m_bc_u.n_dir_edges() * m_xrotrot.numLocalDofsEdge();
  m_unknowns_stride[1].second = m_xrotrot.dimension()
    - m_ddrcore.mesh().n_vertices() * m_xrotrot.numLocalDofsVertex()
    - m_bc_u.n_dir_edges() * m_xrotrot.numLocalDofsEdge();
  // Vertex unknowns for p
  m_unknowns_stride[2].first  = m_xrotrot.dimension()
    + m_bc_p.n_dir_vertices() * m_xgrad.numLocalDofsVertex();
  m_unknowns_stride[2].second =
    ( m_ddrcore.mesh().n_vertices() - m_bc_p.n_dir_vertices() ) * m_xgrad.numLocalDofsVertex();
  // Edge and cell unknowns for p
  m_unknowns_stride[3].first  = m_xrotrot.dimension()
    + m_ddrcore.mesh().n_vertices() * m_xgrad.numLocalDofsVertex()
    + m_bc_p.n_dir_edges() * m_xgrad.numLocalDofsEdge();
  m_unknowns_stride[3].second = m_xgrad.dimension()
    - m_ddrcore.mesh().n_vertices() * m_xgrad.numLocalDofsVertex()
    - m_bc_p.n_dir_edges() * m_xgrad.numLocalDofsEdge();

  // size_t n_unknowns = 0; for (size_t i = 0; i < 4; i++) n_unknowns += m_unknowns_stride[i].second;
  // assert( n_unknowns == dimensionLinearSystem() );
  
  // Create DOFs-to-unknowns map, i.e., a vector with -1 at position i if DOF i
  // is not an unknown, its number as an unknown otherwise
  size_t s = dimensionLinearSystem();
  m_dof_to_unknown = replaceSectionsVector(
                                           (-Eigen::VectorXi::Ones(dimensionSpace())).eval(),
                                           Eigen::VectorXi::LinSpaced(s, 0, s-1).eval(),
                                           m_unknowns_stride
                                           );
}

//------------------------------------------------------------------------------
// Assemble the linear system

Eigen::VectorXd QuadRot::assembleLinearSystem(
                                              const ForcingTermType & f,
                                              const PotentialType & u,
                                              const RotorType & rot_u,
                                              const LagrangeMultiplierType & p
                                              )
{
  auto assemble_all = [this, f, u, p](
                                      size_t start,
                                      size_t end,
                                      std::list<Eigen::Triplet<double> > * triplets_A,
                                      Eigen::VectorXd * rhs_b,
                                      std::list<Eigen::Triplet<double> > * triplets_B
                                      ) -> void
  {
    for (size_t iT = start; iT < end; iT++) {
      this->_assemble_local_contribution(
                                         iT,
                                         this->_compute_local_contribution(iT, f),
                                         *triplets_A,
                                         *rhs_b,
                                         *triplets_B
                                         );
    } // for iT
  };
  
  if (m_use_threads) {
    m_output << "[QuadRot] Parallel assembly" << std::endl;
  }else{
    m_output << "[QuadRot] Sequential assembly" << std::endl;
  }

  std::tie(m_A, m_b, m_B) = parallel_assembly_system(
                                                     m_ddrcore.mesh().n_cells(),
                                                     this->dimensionLinearSystem(),
                                                     std::make_pair(this->dimensionLinearSystem(), this->dimensionSpace()),
                                                     assemble_all,
                                                     m_use_threads
                                                     );

  // Interpolate the lifting of the boundary datum
  Eigen::VectorXd uI = Eigen::VectorXd::Zero(dimensionSpace());
  uI.head(xRotRot().dimension()) = xRotRot().interpolate(u, rot_u);
  uI.tail(xGrad().dimension())   = xGrad().interpolate(p);

  // Add the boundary contribution to the RHS
  m_b -= m_B * uI;
  
  return uI;
}

//------------------------------------------------------------------------------
// Create the vector of DOF indices fort the cell T

std::vector<size_t>
QuadRot::globalDOFIndices(const Cell & T) const
{
  std::vector<size_t> I_xrotrot_T = m_xrotrot.globalDOFIndices(T);
  std::vector<size_t> I_xgrad_T = m_xgrad.globalDOFIndices(T);
  size_t dim_T = m_xrotrot.dimensionCell(T) + m_xgrad.dimensionCell(T);
  size_t dim_xrotrot = m_xrotrot.dimension();
  std::vector<size_t> I_T(dim_T);
  auto it_I_T = std::copy(I_xrotrot_T.begin(), I_xrotrot_T.end(), I_T.begin());
  std::transform(I_xgrad_T.begin(), I_xgrad_T.end(), it_I_T, [&dim_xrotrot](const size_t & index) { return index + dim_xrotrot; });
  return I_T;
}

//------------------------------------------------------------------------------
// Compute local contribution

std::pair<Eigen::MatrixXd, Eigen::VectorXd>
QuadRot::_compute_local_contribution(
                                     size_t iT,
                                     const ForcingTermType & f
                                     )
{
  const Cell & T = *m_ddrcore.mesh().cell(iT);
  
  size_t dim_xrotrot_T = xRotRot().dimensionCell(iT);
  size_t dim_xgrad_T = xGrad().dimensionCell(iT);
  size_t dim_T = dim_xrotrot_T + dim_xgrad_T;

  // Local matrix
  Eigen::MatrixXd MT = Eigen::MatrixXd::Zero(dim_T, dim_T);
  const Eigen::MatrixXd & RT = xRotRot().cellOperators(iT).rotor;
  MT.topLeftCorner(dim_xrotrot_T, dim_xrotrot_T)
    = RT.transpose() * xRot().computeRotorL2Product(iT, m_stab_par, m_stab_par) * RT;
  MT.bottomLeftCorner(dim_xgrad_T, dim_xrotrot_T)
    = xRotRot().computeGradientPotentialL2Product(iT, &m_xgrad, m_stab_par);
  MT.topRightCorner(dim_xrotrot_T, dim_xgrad_T)
    =   MT.bottomLeftCorner(dim_xgrad_T, dim_xrotrot_T).transpose();
  
  // Local vector
  Eigen::VectorXd lT = Eigen::VectorXd::Zero(dim_T);
  QuadratureRule quad_2k_T = generate_quadrature_rule(T, 2 * m_ddrcore.degree());
  lT.head(dim_xrotrot_T)
    = m_xrotrot.cellOperators(iT).potential.transpose()
    * integrate(f, evaluate_quad<Function>::compute(*m_ddrcore.cellBases(iT).Polyk2, quad_2k_T), quad_2k_T);
  
  return std::make_pair(MT, lT);
}

//------------------------------------------------------------------------------
// Assemble the local contribution to the linear system

void QuadRot::_assemble_local_contribution(
                                           size_t iT,
                                           const std::pair<Eigen::MatrixXd, Eigen::VectorXd> & lsT,
                                           std::list<Eigen::Triplet<double> > & triplets_A,
                                           Eigen::VectorXd & rhs_b,
                                           std::list<Eigen::Triplet<double> > & triplets_B
                                           )
{
  const Cell & T = *m_ddrcore.mesh().cell(iT);

  // Create the vector of global DOF indices
  std::vector<size_t> I_T = globalDOFIndices(T);

  // Assemble local contribution
  const Eigen::MatrixXd & MT = lsT.first;
  const Eigen::VectorXd & lT = lsT.second;

  size_t dim_T = dimensionLocalSpace(iT);
  
  for (size_t i = 0; i < dim_T; i++) {
    int unknown_i = m_dof_to_unknown(I_T[i]);
    
    if (unknown_i >= 0) {
      rhs_b(unknown_i) += lT(i);

      for(size_t j = 0; j < dim_T; j++) {
        int unknown_j = m_dof_to_unknown(I_T[j]);
        if (unknown_j >= 0) {
          triplets_A.emplace_back(unknown_i, unknown_j, MT(i, j));
        } else {
          triplets_B.emplace_back(unknown_i, I_T[j], MT(i, j));
        }
      } // for j
    }
  } // for i
}
