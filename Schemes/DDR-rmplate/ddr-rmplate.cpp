// Author: Jerome Droniou (jerome.droniou@monash.edu)
#include <fstream>
#include <iomanip>
#include <thread>

#include "ddr-rmplate.hpp"

#include <mesh_builder.hpp>
#include <parallel_for.hpp>
#include "vtu_writer.hpp"

#include <boost/program_options.hpp>
#include <boost/timer/timer.hpp>

#ifdef WITH_UMFPACK
#include <Eigen/UmfPackSupport>
#endif

#ifdef WITH_MKL
#include <Eigen/PardisoSupport>
#endif

#define FORMAT(W)                                                       \
  std::setiosflags(std::ios_base::left) << std::setw(W) << std::setfill(' ')

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
    ("degree,k", boost::program_options::value<size_t>()->default_value(1), "The polynomial degree of the sequence")
    ("pthread,p", boost::program_options::value<bool>()->default_value(true), "Use thread-based parallelism")
    ("solution,s", boost::program_options::value<int>()->default_value(1), "Select the solution")
    ("thickness,t", boost::program_options::value<double>()->default_value(0.1), "Thickness of the plate")
    ("young-modulus,E", boost::program_options::value<double>()->default_value(1), "Young modulus")
    ("poisson-ratio,n", boost::program_options::value<double>()->default_value(0.2), "Poisson ratio")
    ("boundary-conditions,b", boost::program_options::value<std::string>()->default_value("C"), "Boundary conditions (clamped C, simply supported SS)")
    ("plot", boost::program_options::value<std::string>()->default_value("displacement"), "Save plot of the displacement to the given filename")
    ("export-matrix,e", "Export matrix to Matrix Market format")
    ("iterative-solver,i", "Use iterative linear solver")
    ("stabilization-parameter,x", boost::program_options::value<double>(), "Set the stabilization parameter");

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
  std::cout << FORMAT(25) << "[main] Mesh file" << mesh_file << std::endl;
  
  // Select the degree 
  size_t K = vm["degree"].as<size_t>();
  std::cout << FORMAT(25) << "[main] Degree" << K << std::endl;

  // Select the solution
  int solution = (vm.count("solution") ? vm["solution"].as<int>() : 0);
  ReissnerMindlin::ForcingTermType f;
  ReissnerMindlin::SolutionRotationType theta;
  ReissnerMindlin::GradientRotationType grad_theta;
  ReissnerMindlin::SolutionDisplacementType u;

  switch (solution) {
  case 0:
    std::cout << "[main] Constant solution" << std::endl;    
    f = constant_f;
    theta = constant_theta;
    grad_theta = constant_grad_theta;
    u = constant_u;
    break;

  case 1:
    std::cout << "[main] Polynomial solution" << std::endl;    
    f = polynomial_f;
    theta = polynomial_theta;
    grad_theta = polynomial_grad_theta;
    u = polynomial_u;
    break;

  case 2:
    std::cout << "[main] Analytical solution" << std::endl;    
    f = analytical_f;
    theta = analytical_theta;
    grad_theta = analytical_grad_theta;
    u = analytical_u;
    break;

  case 3:
    std::cout << "[main] Unknown solution" << std::endl;    
    f = ukn_f;
    theta = ukn_theta;
    grad_theta = ukn_grad_theta;
    u = ukn_u;
    break;

  case 4:
    std::cout << "[main] Kirchoff solution (limit t->0)" << std::endl;    
    f = kir_f;
    theta = kir_theta;
    grad_theta = kir_grad_theta;
    u = kir_u;
    break;

  default:
    std::cerr << "[main] ERROR: Unknown exact solution" << std::endl;
    exit(1);
  }

  // Build the mesh and re-order vertices and edges to put the boundary ones first
  MeshBuilder builder = MeshBuilder(mesh_file);
  std::unique_ptr<Mesh> mesh_ptr = builder.build_the_mesh();
  std::cout << FORMAT(25) << "[main] Mesh size" << mesh_ptr->h_max() << std::endl;

  // Select the boundary conditions
  std::string bctype = vm["boundary-conditions"].as<std::string>();
  std::cout << FORMAT(25) << "[main] Boundary cond." << bctype << std::endl;
  // The BC are recorded for the rotation (theta) and the displacement (u), using a pair of BoundaryConditions
  // (first is theta, second is u)
  // Only purely clamped or simply supported BC are valid for the moment; to deal with mixed BC, 
  // the re-ordering of boundary elements will have to be done more specifically
  BoundaryConditions BC_u("D", *mesh_ptr.get());
  BC_u.reorder_edges("start");
  BC_u.reorder_vertices("start");
  std::unique_ptr<BoundaryConditions> BC_theta_pt;
  if (bctype == "C"){
    BC_theta_pt.reset(new BoundaryConditions("D", *mesh_ptr.get()));
  }else if (bctype == "SS"){
    BC_theta_pt.reset(new BoundaryConditions("N", *mesh_ptr.get()));
  }else{
    std::cout << "Boundary conditions " << bctype << " unknown." << std::endl;
    exit(1);
  }
  BoundaryConditions BC_theta = *BC_theta_pt;

  boost::timer::cpu_timer timer;
  // Create DDR core
  timer.start();
  bool use_threads = (vm.count("pthread") ? vm["pthread"].as<bool>() : true);
  std::cout << "[main] " << (use_threads ? "Parallel execution" : "Sequential execution") << std:: endl;
  DDRCore ddr_core(*mesh_ptr, K, use_threads);
  timer.stop();
  double t_wall_ddrcore = double(timer.elapsed().wall) * pow(10, -9);
  double t_proc_ddrcore = double(timer.elapsed().user + timer.elapsed().system) * pow(10, -9);
  std::cout << "[main] Time DDRCore (wall/proc) " << t_wall_ddrcore << "/" << t_proc_ddrcore << std::endl;

  // Assemble the problem
  timer.start();
  double thickness = vm["thickness"].as<double>();
  double young_modulus = vm["young-modulus"].as<double>();
  double poisson_ratio = vm["poisson-ratio"].as<double>();
  std::cout << FORMAT(25) << "[main] Parameters t/E/nu: " << thickness << "/" << young_modulus << "/" << poisson_ratio << std::endl;
  RMParameters para(thickness, young_modulus, poisson_ratio);
  
  ReissnerMindlin rm(ddr_core, para, BC_theta, BC_u, use_threads);
  if(vm.count("stabilization-parameter")) {
    rm.stabilizationParameter() = vm["stabilization-parameter"].as<double>();
  }
  rm.assembleLinearSystem(f, theta, grad_theta, u);
  
  timer.stop();
  double t_wall_model = double(timer.elapsed().wall) * pow(10, -9);
  double t_proc_model = double(timer.elapsed().user + timer.elapsed().system) * pow(10, -9);
  std::cout << "[main] Time model (wall/proc) " << t_wall_model << "/" << t_proc_model << std::endl;

  // Export matrix if requested  
  if (vm.count("export-matrix")) {
    std::cout << "[main] Exporting matrix to Matrix Market format" << std::endl;
    saveMarket(rm.systemMatrix(), "A_rm.mtx");
    saveMarket(rm.systemVector(), "b_rm.mtx");
  }

  // Solve the problem
  timer.start();
  Eigen::VectorXd uh_int;
  if (vm.count("iterative-solver")) {
    std::cout << "[main] Solving the linear system using BiCGSTAB" << std::endl;
    
    Eigen::BiCGSTAB<ReissnerMindlin::SystemMatrixType, Eigen::IncompleteLUT<double> > solver;
    // solver.preconditioner().setFillfactor(2);
    solver.compute(rm.systemMatrix());
    if (solver.info() != Eigen::Success) {
      std::cerr << "[main] ERROR: Could not factorize matrix" << std::endl;
      exit(1);
    }
    uh_int = solver.solve(rm.systemVector());
    if (solver.info() != Eigen::Success) {
      std::cerr << "[main] ERROR: Could not solve direct system" << std::endl;
      exit(1);
    }
  } else { 
#ifdef WITH_MKL
    std::cout << "[main] Solving the linear system using Pardiso" << std::endl;    
    Eigen::PardisoLU<ReissnerMindlin::SystemMatrixType> solver;
#elif WITH_UMFPACK
    std::cout << "[main] Solving the linear system using Umfpack" << std::endl;    
    Eigen::UmfPackLU<ReissnerMindlin::SystemMatrixType> solver;
#else
    std::cout << "[main] Solving the linear system using direct solver" << std::endl;    
    Eigen::SparseLU<ReissnerMindlin::SystemMatrixType> solver;
#endif
    solver.compute(rm.systemMatrix());
    if (solver.info() != Eigen::Success) {
      std::cerr << "[main] ERROR: Could not factorize matrix" << std::endl;
    }
    uh_int = solver.solve(rm.systemVector());
    if (solver.info() != Eigen::Success) {
      std::cerr << "[main] ERROR: Could not solve linear system" << std::endl;
    }
  }
  // Re-create complete vector with BC
  Eigen::VectorXd uh = replaceSectionsVector(rm.bdryValues(), uh_int, rm.locUKN());

  timer.stop();
  double t_wall_solve = double(timer.elapsed().wall) * pow(10, -9);
  double t_proc_solve = double(timer.elapsed().user + timer.elapsed().system) * pow(10, -9);
  std::cout << "[main] Time solve (wall/proc) " << t_wall_solve << "/" << t_proc_solve << std::endl;

  // Compute energy error and error of displacement at center of plate
  Eigen::VectorXd uI = Eigen::VectorXd::Zero(rm.dimensionSpace());  
  uI.head(rm.exCurl().dimension()) = rm.exCurl().interpolate(rm.contractPara<Vector2d>(theta));
  uI.tail(rm.xGrad().dimension()) = rm.xGrad().interpolate(rm.contractPara<double>(u));

  Eigen::VectorXd eh = uh - uI;
  double en_error = rm.computeNorm(eh) / rm.computeNorm(uI);
  
  VectorRd center(0.,0.);
  double area = 0.;
  for (Cell * T : mesh_ptr->get_cells()){
    area += T->measure();
    center += T->measure() * T->center_mass();
  }
  center /= area;
  double exact_dis = u(para, center);
  size_t iT_center = mesh_ptr->find_cell(center);
  Eigen::VectorXd uT = rm.xGrad().restrictCell(iT_center, uh.tail(rm.xGrad().dimension()));
  double approx_dis = rm.xGrad().evaluatePotential(iT_center, uT, center);
  double dis_error = std::abs(approx_dis - exact_dis) / std::abs(1e-15 + exact_dis);
  
  std::cout << "[main] Energy error " << en_error << std::endl;
  std::cout << "[main] Error displacement center " << dis_error << std::endl;
  std::cout << "[main] Mesh diameter " << mesh_ptr->h_max() << std::endl;
 
  // --------------------------------------------------------------------------
  //                     Creates a .vtu file of the displacement
  // --------------------------------------------------------------------------
  // Only if we do not have too many cells
  if (vm.count("plot") && mesh_ptr->n_cells() <= max_nb_cells_for_plot) {
    std::cout << "[Main] Writing solution to file" << std::endl;
    std::string filename = vm["plot"].as<std::string>() + std::string(".vtu");
    VtuWriter plotdata(mesh_ptr.get());

    plotdata.write_to_vtu(std::string("exact-")+filename,uI.segment(rm.exCurl().dimension(), mesh_ptr->n_vertices()));
    plotdata.write_to_vtu(filename,uh.segment(rm.exCurl().dimension(), mesh_ptr->n_vertices()));

  }

  // --------------------------------------------------------------------------
  //                     Creates .txt file with data and results
  // --------------------------------------------------------------------------
  std::ofstream out("results.txt");
  out << "Solution: " << solution << std::endl;
  out << "Boundary conditions: " << bctype << std::endl;
  out << "Mesh: " << mesh_file << std::endl;
  out << "Degree: " << K << std::endl;
  out << "MeshSize: " << mesh_ptr->h_max() << std::endl;
  out << "NbCells: " << mesh_ptr->n_cells() << std::endl;
  out << "NbEdges: " << mesh_ptr->n_edges() << std::endl;
  out << "NbVertices: " << mesh_ptr->n_vertices() << std::endl;
  out << "DimEXCurl: " << rm.exCurl().dimension() << std::endl;
  out << "DimXGrad: " << rm.xGrad().dimension() << std::endl;
  out << "DimSystem: " << rm.sizeSystem() << std::endl;
  out << "EnError: " << en_error << std::endl;
  out << "DisError: " << dis_error << std::endl;
  out << "thickness (t): " << para.t << std::endl;
  out << "Young modulus (E): " << para.E << std::endl;
  out << "Poisson ratio (nu): " << para.nu << std::endl;
  out << "TwallDDRCore: " << t_wall_ddrcore << std::endl;  
  out << "TprocDDRCore: " << t_proc_ddrcore << std::endl;  
  out << "TwallModel: " << t_wall_model << std::endl;  
  out << "TprocModel: " << t_proc_model << std::endl;  
  out << "TwallSolve: " << t_wall_solve << std::endl;  
  out << "TprocSolve: " << t_proc_solve << std::endl;  
  out << std::flush;
  out.close();

  std::cout << "[main] Done" << std::endl;
  return 0;
}

//------------------------------------------------------------------------------
// ReissnerMindlin
//------------------------------------------------------------------------------

ReissnerMindlin::ReissnerMindlin(
                               const DDRCore & ddrcore,
                               const RMParameters & para,
                               const BoundaryConditions & BC_theta,
                               const BoundaryConditions & BC_u,
                               bool use_threads,
                               std::ostream & output
                               )
  : m_ddrcore(ddrcore),
    m_para(para),
    m_use_threads(use_threads),
    m_output(output),
    m_excurl(ddrcore, use_threads),
    m_xgrad(ddrcore, use_threads),
    m_BC_theta(BC_theta),
    m_BC_u(BC_u),
    m_bdryMatrix(sizeSystem(), dimensionSpace()),
    m_bdryValues(Eigen::VectorXd::Zero(dimensionSpace())),
    m_A(sizeSystem(), sizeSystem()),
    m_b(Eigen::VectorXd::Zero(sizeSystem())),
    m_stab_par(1.)
{  
  m_output << "[ReissnerMindlin] Initializing" << std::endl;
  
  // Locations of the unknowns among the DOFs: they correspond to the sections (for all i) of DOFs starting at m_locUKN[i].first and having length m_locUKN[i].second
  m_locUKN.resize(3);
  // Edge and cell unknowns for theta: all edges after the Dirichlet edges, and all cells
  m_locUKN[0].first = m_BC_theta.n_dir_edges() * m_excurl.numLocalDofsEdge();
  m_locUKN[0].second = m_excurl.dimension() - m_BC_theta.n_dir_edges() * m_excurl.numLocalDofsEdge();
  // Vertex unknowns for u: all vertices after the Dirichlet vertices
  m_locUKN[1].first = m_excurl.dimension() + m_BC_u.n_dir_vertices() * m_xgrad.numLocalDofsVertex();
  m_locUKN[1].second = (m_ddrcore.mesh().n_vertices() - m_BC_u.n_dir_vertices()) * m_xgrad.numLocalDofsVertex();
  // Edge and cell unknowns for u: all edges after the Dirichlet ones, and all cells
  m_locUKN[2].first = m_excurl.dimension() + m_ddrcore.mesh().n_vertices() * m_xgrad.numLocalDofsVertex() + m_BC_u.n_dir_edges() * m_xgrad.numLocalDofsEdge();
  m_locUKN[2].second = dimensionSpace() - m_locUKN[2].first;
  
  // Create the map DOF->unknowns: vector with -1 at position i if DOF i is not an unknown, otherwise its unknown number
  const size_t N=dimensionSpace();
  const size_t s=sizeSystem();
  m_DOFtoUKN = replaceSectionsVector((-Eigen::VectorXi::Ones(N)).eval(), Eigen::VectorXi::LinSpaced(s, 0, s-1).eval(), m_locUKN);
  
}

//------------------------------------------------------------------------------
//      Assemble global linear system and BC
//------------------------------------------------------------------------------

void ReissnerMindlin::assembleLinearSystem(
                                           const ForcingTermType & f,
                                           const SolutionRotationType & theta,
                                           const GradientRotationType & grad_theta,
                                           const SolutionDisplacementType & u
                                           )

{
  // Function to assemble all local contributions (without low-order stabilisation)
  auto assemble_all = [this, f, theta, grad_theta, u](
                                       size_t start,
                                       size_t end,
                                       std::list<Eigen::Triplet<double> > * triplets_A,
                                       Eigen::VectorXd * rhs_b,
                                       std::list<Eigen::Triplet<double> > * triplets_bdryMatrix
                                       )->void
                      {
                        for (size_t iT = start; iT < end; iT++) {
                          this->_assemble_local_contribution(
                                                             iT,
                                                             this->_compute_local_contribution(iT, grad_theta, f),
                                                             *triplets_A,
                                                             *rhs_b,
                                                             *triplets_bdryMatrix
                                                             );
                        } // for iT
                      };
                      
                      
  // Assemble low-order stabilisations
  auto assemble_all_jump_stab = [this](
                                     size_t start,
                                     size_t end,
                                     std::list<Eigen::Triplet<double> > * triplets_jump_stab,
                                     Eigen::VectorXd * dummy_rhs,
                                     std::list<Eigen::Triplet<double> > * triplets_bdryMatrix
                                     )->void
                      {
                        for (size_t iE = start; iE < end; iE++) {
                          this->_assemble_local_jump_stab(
                                                 iE,
                                                 this->_compute_local_jump_stab(iE),
                                                 *triplets_jump_stab,
                                                 *triplets_bdryMatrix
                                                 );
                        } // for iE
                      };

  // Assemble the matrix and rhs (ordering of unknowns: theta, u and (as per ddrspace) vertices, edges, cells for each
  if (m_use_threads) {
    m_output << "[ReissnerMindlin] Parallel assembly" << std::endl;
  }else{
    m_output << "[ReissnerMindlin] Sequential assembly" << std::endl;
  }

  std::tie(m_A, m_b, m_bdryMatrix) = parallel_assembly_system(m_ddrcore.mesh().n_cells(), this->sizeSystem(), std::make_pair(this->sizeSystem(), this->dimensionSpace()), assemble_all, m_use_threads);

  // Adding jump terms if k=0
  if (m_ddrcore.degree()==0){
    SystemMatrixType jump_A(sizeSystem(), sizeSystem());
    SystemMatrixType jump_bdry(sizeSystem(), dimensionSpace());
    Eigen::VectorXd dummy = Eigen::VectorXd::Zero(this->sizeSystem());
    
    std::tie(jump_A, dummy, jump_bdry) = parallel_assembly_system(m_ddrcore.mesh().n_edges(), this->sizeSystem(), std::make_pair(this->sizeSystem(), this->dimensionSpace()), assemble_all_jump_stab, m_use_threads);
    
    double scaling = 1e-1;
    m_A += scaling * jump_A;
    m_bdryMatrix += scaling * jump_bdry;
  }

  // Create boundary values 
  for (Edge *E : m_ddrcore.mesh().get_b_edges()) {
    if (m_BC_theta.type(*E)=="dir"){
      // Dirichlet BC on theta (values in m_bdryValues): normal components first, then tangential
      QuadratureRule quad_2k_E = generate_quadrature_rule(*E, 2 * m_ddrcore.degree());    
      size_t iE = E->global_index();

      size_t dim_Pk_E = m_ddrcore.edgeBases(iE).Polyk->dimension();
      auto basis_Pk_E_quad = evaluate_quad<Function>::compute(*m_ddrcore.edgeBases(iE).Polyk, quad_2k_E);
      Eigen::Vector2d nE = E->normal();
      auto theta_dot_nE = [this, &nE, &theta](const Eigen::Vector2d & x)->double {
	          return theta(m_para, x).dot(nE);
            };
      m_bdryValues.segment(m_excurl.globalOffset(*E), dim_Pk_E)
	          = l2_projection(theta_dot_nE, *m_ddrcore.edgeBases(iE).Polyk, quad_2k_E, basis_Pk_E_quad);

      Eigen::Vector2d tE = E->tangent();
      auto theta_dot_tE = [this, &tE, &theta](const Eigen::Vector2d & x)->double {
	          return theta(m_para, x).dot(tE);
            };
      m_bdryValues.segment(m_excurl.globalOffset(*E) + dim_Pk_E, dim_Pk_E)
	          = l2_projection(theta_dot_tE, *m_ddrcore.edgeBases(iE).Polyk, quad_2k_E, basis_Pk_E_quad);    
    }

    if (m_BC_u.type(*E)=="dir" && (m_ddrcore.degree()>0)){
      // Dirichlet BC for u (values in m_bdryValues)
      QuadratureRule quad_2k_E = generate_quadrature_rule(*E, 2 * m_ddrcore.degree());    
      size_t iE = E->global_index();

      auto basis_Pkmo_E_quad = evaluate_quad<Function>::compute(*m_ddrcore.edgeBases(iE).Polykmo, quad_2k_E);
      m_bdryValues.segment(m_excurl.dimension() + m_xgrad.globalOffset(*E), m_xgrad.numLocalDofsEdge())
            += l2_projection(contractPara<double>(u), *m_ddrcore.edgeBases(iE).Polykmo, quad_2k_E, basis_Pkmo_E_quad);
    }   
  } // for iE
  
  for (Vertex *V : m_ddrcore.mesh().get_b_vertices()) {
    if (m_BC_u.type(*V)=="dir"){
      // vertex BC on u
      m_bdryValues(m_excurl.dimension() + m_xgrad.globalOffset(*V)) = u(m_para, V->coords());
    }
  }
  
  // Adjust the system's RHS to account for them
  m_b -= m_bdryMatrix * m_bdryValues;
  
}

//------------------------------------------------------------------------------
//      Compute and assemble local contributions (no jump stabilisation)
//------------------------------------------------------------------------------
 
std::pair<Eigen::MatrixXd, Eigen::VectorXd>
ReissnerMindlin::_compute_local_contribution(size_t iT, const GradientRotationType & grad_theta, const ForcingTermType & f)
{
  const Cell & T = *m_ddrcore.mesh().cell(iT);

  size_t dim_excurl_T = m_excurl.dimensionCell(iT);
  size_t dim_xgrad_T = m_xgrad.dimensionCell(iT);
  size_t dim_T = dim_excurl_T + dim_xgrad_T;
  
  Eigen::MatrixXd AT = Eigen::MatrixXd::Zero(dim_T, dim_T);
  Eigen::VectorXd lT = Eigen::VectorXd::Zero(dim_T);

  //----------------------------
  // Local matrix

  // Mass matrix for Polyk2x2
  QuadratureRule quad_2k_T = generate_quadrature_rule(T, 2 * m_ddrcore.degree());
  auto basis_Polyk2x2_quad = evaluate_quad<Function>::compute(m_excurl.Polyk2x2(iT), quad_2k_T);
  Eigen::MatrixXd mass_Pk2x2_T = compute_gram_matrix(basis_Polyk2x2_quad, quad_2k_T);
  
  // a_h: contribution symmetric gradient and stabilisation
  Eigen::MatrixXd GsT = m_excurl.hhoOperators(iT).symmetric_gradient;
  AT.topLeftCorner(dim_excurl_T, dim_excurl_T) += 
        m_para.beta0 * ( GsT.transpose() * mass_Pk2x2_T * GsT + m_excurl.hhoOperators(iT).stabilisation );
  
  // a_h: contribution divergence
  auto basis_trace_Polyk2x2_quad = transform_values_quad<double>(basis_Polyk2x2_quad,
                                                                 [](const Eigen::Matrix2d &m)->double { return m.trace();});
  Eigen::MatrixXd mass_trace_Polyk2x2 = compute_gram_matrix(basis_trace_Polyk2x2_quad, quad_2k_T);
  AT.topLeftCorner(dim_excurl_T, dim_excurl_T) +=
        m_para.beta1 * ( GsT.transpose() * mass_trace_Polyk2x2 * GsT );

  // b_h: contribution split by developping (tau-G w,eta-Gv)=(tau,eta) - (Gw, eta) - (tau,Gv) + (Gw, Gv)
  const double kappa_tm2 = m_para.kappa * std::pow(m_para.t, -2);
  Eigen::MatrixXd mass_Pk2_T = 
      compute_gram_matrix(evaluate_quad<Function>::compute(*m_excurl.cellBases(iT).Polyk2, quad_2k_T), quad_2k_T);
  AT.topLeftCorner(dim_excurl_T, dim_excurl_T) += kappa_tm2 * m_excurl.computeL2Product(iT, m_stab_par, mass_Pk2_T);
  
  Eigen::MatrixXd L2prod_grad_left = m_excurl.computeL2ProductGradient(iT, m_xgrad, "left", m_stab_par, mass_Pk2_T);
  AT.bottomLeftCorner(dim_xgrad_T, dim_excurl_T) -= kappa_tm2 *  L2prod_grad_left;
  AT.topRightCorner(dim_excurl_T, dim_xgrad_T) -= kappa_tm2 * L2prod_grad_left.transpose();
  
  Eigen::MatrixXd L2prod_grad_both = m_excurl.computeL2ProductGradient(iT, m_xgrad, "both", m_stab_par, mass_Pk2_T);
  AT.bottomRightCorner(dim_xgrad_T, dim_xgrad_T) += kappa_tm2 * L2prod_grad_both;
    
  //------------------------
  // Local RHS: load
  QuadratureRule quad_2kpt_T = generate_quadrature_rule(T, 2 * m_ddrcore.degree() + 2);
  lT.tail(dim_xgrad_T) = m_xgrad.cellOperators(iT).potential.transpose()
    * integrate(contractPara<double>(f), evaluate_quad<Function>::compute(*m_ddrcore.cellBases(iT).Polykpo, quad_2kpt_T), quad_2kpt_T);
    
 
  // Local RHS: simply supported BC
  if (T.is_boundary()){
    for (size_t iE=0; iE < T.n_edges(); iE++){
      Edge & E = *T.edge(iE);
      if (m_BC_theta.type(E)=="neu"){
        VectorRd nTE = T.edge_normal(iE);
        VectorRd nE = E.normal();
        VectorRd tE = E.tangent();
        
        std::function<double(const VectorRd&)> strain_nTE_dot_nE = [this, &nTE, &nE, &grad_theta](const VectorRd &x)->double {
            return (m_para.beta0 * symmetrise_matrix(grad_theta(m_para, x)) * nTE + m_para.beta1 * (grad_theta(m_para, x).trace()) * nTE).dot(nE);
          };
        std::function<double(const VectorRd&)> strain_nTE_dot_tE = [this, &nTE, &tE, &grad_theta](const VectorRd &x)->double {
            return (m_para.beta0 * symmetrise_matrix(grad_theta(m_para, x)) * nTE + m_para.beta1 * (grad_theta(m_para, x).trace()) * nTE).dot(tE);
          };

        auto basis_Pk_E = *m_ddrcore.edgeBases(E.global_index()).Polyk;
        size_t dim_Pk_E = basis_Pk_E.dimension();
        auto qr_2kpt_E = generate_quadrature_rule(E, 2*m_ddrcore.degree() + 2);
        auto basis_Pk_E_quad = evaluate_quad<Function>::compute(basis_Pk_E, qr_2kpt_E);
        lT.segment(m_excurl.localOffset(T, E), dim_Pk_E) = integrate(strain_nTE_dot_nE, basis_Pk_E_quad, qr_2kpt_E);
        lT.segment(m_excurl.localOffset(T, E) + dim_Pk_E, dim_Pk_E) = integrate(strain_nTE_dot_tE, basis_Pk_E_quad, qr_2kpt_E);
      }
    }
  }
  
  return std::make_pair(AT, lT);
}

//------------------------------------------------------------------------------

void ReissnerMindlin::_assemble_local_contribution(
                                                  size_t iT,
                                                  const std::pair<Eigen::MatrixXd, Eigen::VectorXd> & lsT,
                                                  std::list<Eigen::Triplet<double> > & triplets_A,
                                                  Eigen::VectorXd & rhs_b,
                                                  std::list<Eigen::Triplet<double> > & triplets_bdryMat
                                                  )
{
  const Cell & T = *m_ddrcore.mesh().cell(iT);

  // Create the vector of global DOF indices
  std::vector<size_t> I_T = globalDOFIndices(T);

  // Assemble
  const Eigen::MatrixXd & AT = lsT.first;
  const Eigen::VectorXd & bT = lsT.second;
  size_t dim_T = m_excurl.dimensionCell(iT) + m_xgrad.dimensionCell(iT);
  for (size_t i = 0; i < dim_T; i++) {
    int ukn_i = m_DOFtoUKN(I_T[i]);
    if (ukn_i>=0){
      rhs_b(ukn_i) += bT(i);
      for (size_t j = 0; j < dim_T; j++) {
        int ukn_j = m_DOFtoUKN(I_T[j]);
        if (ukn_j>=0){
          triplets_A.emplace_back(ukn_i, ukn_j, AT(i,j));
        }else{
          triplets_bdryMat.emplace_back(ukn_i, I_T[j], AT(i,j));
        }
      }
    } // for j
  } // for i

}

//------------------------------------------------------------------------------
// Compute and assemble local jump stabilisation (for low-order)
//------------------------------------------------------------------------------
 
Eigen::MatrixXd ReissnerMindlin::_compute_local_jump_stab(size_t iE)
{
  const Edge & E = *m_ddrcore.mesh().edge(iE);
  
  // Different definition of jump for interior or boundary edges
  if (!E.is_boundary()){
    // Cells on each side
    Cell & T1 = *E.get_cells()[0];
    Cell & T2 = *E.get_cells()[1];
    size_t dim_T1 = m_excurl.dimensionCell(T1.global_index());
    size_t dim_T2 = m_excurl.dimensionCell(T2.global_index());

    // Matrix corresponding to jump bilinear form
    Eigen::MatrixXd AT = Eigen::MatrixXd::Zero(dim_T1 + dim_T2, dim_T1 + dim_T2);
    
    // Compute jump penalisation by developping the jumps.
    auto qr_E_2kpt = generate_quadrature_rule(E, 2*(m_ddrcore.degree()+1));
    auto basis_Pkpo2_T1_quad = evaluate_quad<Function>::compute(m_excurl.Polykpo2(T1.global_index()), qr_E_2kpt);
    auto basis_Pkpo2_T2_quad = evaluate_quad<Function>::compute(m_excurl.Polykpo2(T2.global_index()), qr_E_2kpt);
    Eigen::MatrixXd MT1T1 = compute_gram_matrix(basis_Pkpo2_T1_quad, qr_E_2kpt);
    Eigen::MatrixXd MT2T2 = compute_gram_matrix(basis_Pkpo2_T2_quad, qr_E_2kpt);
    Eigen::MatrixXd MT1T2 = compute_gram_matrix(basis_Pkpo2_T1_quad, basis_Pkpo2_T2_quad, qr_E_2kpt);

    Eigen::MatrixXd pT1 = m_excurl.hhoOperators(T1.global_index()).potential;
    Eigen::MatrixXd pT2 = m_excurl.hhoOperators(T2.global_index()).potential;

    AT.topLeftCorner(dim_T1, dim_T1) = pT1.transpose() * MT1T1 * pT1;
    AT.topRightCorner(dim_T1, dim_T2) = -pT1.transpose() * MT1T2 * pT2;
    AT.bottomLeftCorner(dim_T2, dim_T1) = -pT2.transpose() * MT1T2.transpose() * pT1;
    AT.bottomRightCorner(dim_T2, dim_T2) = pT2.transpose() * MT2T2 * pT2;

    return (1./E.diam()) * AT;
  }else if (m_BC_theta.type(E) != "neu"){
    Cell & T = *E.get_cells()[0];
    auto qr_E_2kpt = generate_quadrature_rule(E, 2*(m_ddrcore.degree()+1));
    auto basis_Pkpo2_T_quad = evaluate_quad<Function>::compute(m_excurl.Polykpo2(T.global_index()), qr_E_2kpt);
    Eigen::MatrixXd MTT = compute_gram_matrix(basis_Pkpo2_T_quad, qr_E_2kpt);
    Eigen::MatrixXd pT = m_excurl.hhoOperators(T.global_index()).potential;
    
    return pT.transpose() * MTT * pT;
  }else{
    size_t iT = E.get_cells()[0]->global_index();
    size_t dim_T = m_excurl.dimensionCell(iT);
    return Eigen::MatrixXd::Zero(dim_T, dim_T);
  }
}

//------------------------------------------------------------------------------

void ReissnerMindlin::_assemble_local_jump_stab(
                                              size_t iE,
                                              const Eigen::MatrixXd & locJ,
                                              std::list<Eigen::Triplet<double> > & triplets_jump_stab,
                                              std::list<Eigen::Triplet<double> > & triplets_bdryMat
                                              )
{
  const Edge & E = *m_ddrcore.mesh().edge(iE);

  // Different set of DOFs for interior or boundary edges
  if (!E.is_boundary()){
    // Create the vector of global DOF indices: DOFs of each element around E
    Cell & T1 = *E.get_cells()[0];
    Cell & T2 = *E.get_cells()[1];
    size_t dim_T1 = m_excurl.dimensionCell(T1.global_index());
    size_t dim_T2 = m_excurl.dimensionCell(T2.global_index());
    std::vector<size_t> I_T1 = m_excurl.globalDOFIndices(T1);
    std::vector<size_t> I_T2 = m_excurl.globalDOFIndices(T2);
    std::vector<size_t> I_T(dim_T1+dim_T2);
    auto it_I_T = std::copy(I_T1.begin(), I_T1.end(), I_T.begin());
    std::copy(I_T2.begin(), I_T2.end(), it_I_T);
    
    // Assemble
    for (size_t i = 0; i < dim_T1 + dim_T2; i++) {
      int ukn_i = m_DOFtoUKN(I_T[i]);
      if (ukn_i>=0){
        for (size_t j = 0; j < dim_T1 + dim_T2; j++) {
          int ukn_j = m_DOFtoUKN(I_T[j]);
          if (ukn_j>=0){
            triplets_jump_stab.emplace_back(ukn_i, ukn_j, locJ(i,j));
          }else{
            triplets_bdryMat.emplace_back(ukn_i, I_T[j], locJ(i,j));
          }
        }
      } // for j
    } // for i
  }else{
    Cell & T = *E.get_cells()[0];
    size_t dim_T = m_excurl.dimensionCell(T.global_index());
    std::vector<size_t> I_T = m_excurl.globalDOFIndices(T);
    
    // Assemble
    for (size_t i = 0; i < dim_T; i++) {
      int ukn_i = m_DOFtoUKN(I_T[i]);
      if (ukn_i>=0){
        for (size_t j = 0; j < dim_T; j++) {
          int ukn_j = m_DOFtoUKN(I_T[j]);
          if (ukn_j>=0){
            triplets_jump_stab.emplace_back(ukn_i, ukn_j, locJ(i,j));
          }else{
            triplets_bdryMat.emplace_back(ukn_i, I_T[j], locJ(i,j));
          }
        }
      } // for j
    } // for i
  
  }
  
}

//------------------------------------------------------------------------------
double ReissnerMindlin::computeNorm(const Eigen::VectorXd & v ) const
{
  Eigen::VectorXd local_sqnorms = Eigen::VectorXd::Zero(m_ddrcore.mesh().n_cells());

  // EXcurl correspond to the first components of v, XGrad to the last ones
  Eigen::VectorXd v_curl = v.head(m_excurl.dimension());
  Eigen::VectorXd v_grad = v.tail(m_xgrad.dimension());
  
  std::function<void(size_t, size_t)> compute_local_squarednorms
    = [this, &v_curl, &v_grad, &local_sqnorms](size_t start, size_t end)->void
    {
      for (size_t iT = start; iT < end; iT++){
        Cell & T = *m_ddrcore.mesh().cell(iT);
        Eigen::VectorXd v_curl_T = m_excurl.restrict(T, v_curl);
        Eigen::VectorXd v_grad_T = m_xgrad.restrict(T, v_grad);

        QuadratureRule quad_2k_T = generate_quadrature_rule(T, 2*m_ddrcore.degree() );
        
        // Contribution of GsT and sT
        auto basis_Pk2x2_T_quad = evaluate_quad<Function>::compute(m_excurl.Polyk2x2(iT), quad_2k_T);
        Eigen::MatrixXd mass_Pk2x2_T = compute_gram_matrix(basis_Pk2x2_T_quad, quad_2k_T);
        Eigen::VectorXd GsT_v_curl_T = m_excurl.hhoOperators(iT).symmetric_gradient * v_curl_T;
        double tmp = GsT_v_curl_T.transpose() * mass_Pk2x2_T * GsT_v_curl_T;
        tmp += v_curl_T.transpose() * m_excurl.hhoOperators(iT).stabilisation * v_curl_T;
        local_sqnorms(iT) += m_para.beta0 * tmp;

        // Contribution divergence
        auto basis_trace_Pk2x2_T_quad = transform_values_quad<double>(basis_Pk2x2_T_quad,
                                                                 [](const Eigen::Matrix2d &m)->double { return m.trace();});
        Eigen::MatrixXd mass_trace_Polyk2x2 = compute_gram_matrix(basis_trace_Pk2x2_T_quad, quad_2k_T);
        tmp = GsT_v_curl_T.transpose() * mass_trace_Polyk2x2 * GsT_v_curl_T;
        local_sqnorms(iT) += m_para.beta1 * tmp;

        // Contribution norme EXcurl
        Eigen::MatrixXd mass_Pk2_T = compute_gram_matrix(evaluate_quad<Function>::compute(*m_ddrcore.cellBases(iT).Polyk2, quad_2k_T), quad_2k_T);
        Eigen::MatrixXd excurl_L2product = m_excurl.computeL2Product(iT, m_stab_par, mass_Pk2_T);
        Eigen::MatrixXd excurl_L2product_grad = m_excurl.computeL2ProductGradient(iT, m_xgrad, "both", m_stab_par, mass_Pk2_T);
        tmp = v_curl_T.transpose() * excurl_L2product * v_curl_T;
        tmp += v_grad_T.transpose() * excurl_L2product_grad * v_grad_T;
        local_sqnorms(iT) += std::min(m_para.kappa, m_para.beta0) * tmp;

        // Contribution Kirchoff term (developed)
        tmp =  v_curl_T.transpose() * excurl_L2product * v_curl_T;
        tmp += v_grad_T.transpose() * excurl_L2product_grad * v_grad_T;
        Eigen::MatrixXd excurl_L2product_grad_right = m_excurl.computeL2ProductGradient(iT, m_xgrad, "right", m_stab_par, mass_Pk2_T);
        tmp -= 2. * v_curl_T.transpose() * excurl_L2product_grad_right * v_grad_T;
        local_sqnorms(iT) += (m_para.kappa / std::pow(m_para.t, 2)) * tmp;
        

        // Contribution of L2 norms
//        local_sqnorms(iT) += v_curl_T.transpose() * m_excurl.computeL2Product(iT, m_stab_par, mass_Pk2_T) * v_curl_T;
//        local_sqnorms(iT) += v_grad_T.transpose() * m_xgrad.computeL2Product(iT, m_stab_par, mass_Pk2_T) * v_grad_T;

      }
    };
  parallel_for(m_ddrcore.mesh().n_cells(), compute_local_squarednorms, m_use_threads);
  
  return std::sqrt(std::abs(local_sqnorms.sum()));
}
