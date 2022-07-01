// Authors: Daniele Di Pietro (daniele.di-pietro@umontpellier.fr) and Jerome Droniou (jerome.droniou@monash.edu)
#include <fstream>
#include <iomanip>
#include <thread>

#include "ddr-klplate.hpp"
#include <GMpoly_cell.hpp>

#include <mesh_builder.hpp>
#include <parallel_for.hpp>
#include "vtu_writer.hpp"

#include <boost/program_options.hpp>
#include <display_timer.hpp>

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
    ("degree,k", boost::program_options::value<size_t>()->default_value(3), "The polynomial degree of the sequence")
    ("pthread,p", boost::program_options::value<bool>()->default_value(true), "Use thread-based parallelism")
    ("solution,s", boost::program_options::value<int>()->default_value(0), "Select the solution")
    ("bending-modulus,D", boost::program_options::value<double>()->default_value(1.), "Bending modulus")
    ("poisson-ratio,n", boost::program_options::value<double>()->default_value(0.3), "Poisson ratio")
    ("plot", boost::program_options::value<std::string>()->default_value("displacement"), "Save plot of the displacement to the given filename")
    ("export-matrix,e", boost::program_options::value<bool>()->default_value(false), "Export matrix to Matrix Market format")
    ("iterative-solver,i", "Use iterative linear solver")
    ("stabilization-parameter,x", boost::program_options::value<double>(), "Set the stabilization parameter");;
  
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
  size_t K = (vm.count("degree") ? vm["degree"].as<size_t>() : 3);
  std::cout << FORMAT(25) << "[main] Degree" << K << std::endl;

  // Select the solution
  int solution = (vm.count("solution") ? vm["solution"].as<int>() : 0);
  KirchhoffLove::DeflectionType u;
  KirchhoffLove::GradientDeflectionType grad_u;
  KirchhoffLove::MomentTensorType hess_u;
  KirchhoffLove::GradientDeflectionType grad_Delta_u;
  KirchhoffLove::MomentTensorEdgeDerivativeType hess_u_DE;
  KirchhoffLove::ForcingTermType divdiv_hess_u;
  
  switch (solution) {
  case 0:
    std::cout << "[main] Trigonometric solution" << std::endl;    
    u = trigonometric_u;
    grad_u = trigonometric_grad_u;
    hess_u = trigonometric_hess_u;
    grad_Delta_u = trigonometric_grad_Delta_u;
    hess_u_DE = trigonometric_hess_u_DE;
    divdiv_hess_u = trigonometric_divdiv_hess_u;
    break;

  case 1:
    std::cout << "[main] Quartic solution" << std::endl;    
    u = quartic_u;
    grad_u = quartic_grad_u;
    hess_u = quartic_hess_u;
    grad_Delta_u = quartic_grad_Delta_u;
    hess_u_DE = quartic_hess_u_DE;
    divdiv_hess_u = quartic_divdiv_hess_u;
    break;

  case 2:
    std::cout << "[main] Biquartic solution" << std::endl;    
    u = biquartic_u;
    grad_u = biquartic_grad_u;
    hess_u = biquartic_hess_u;
    grad_Delta_u = biquartic_grad_Delta_u;
    hess_u_DE = biquartic_hess_u_DE;
    divdiv_hess_u = biquartic_divdiv_hess_u;
    break;

  default:
    std::cerr << "[main] ERROR: Unknown exact solution" << std::endl;
    exit(1);
  }
  double D = vm["bending-modulus"].as<double>();
  double nu = vm["poisson-ratio"].as<double>();
  KirchhoffLove::MomentTensorType sigma = 
      [&D, &nu, &hess_u](const VectorRd & x)->MatrixRd {
          return -D*( (1.-nu)*hess_u(x) + nu * hess_u(x).trace()*MatrixRd::Identity());
        };
  KirchhoffLove::MomentTensorEdgeDerivativeType sigma_DE = 
        [&D, &nu, &hess_u_DE, &grad_Delta_u](const VectorRd & x, const Edge & E)->double {
          return -D*(1.-nu)*hess_u_DE(x, E) - D*nu*grad_Delta_u(x).dot(E.normal());
        };
  KirchhoffLove::ForcingTermType f = 
        [&D, &divdiv_hess_u](const VectorRd & x)->double {
          return D*divdiv_hess_u(x);
        };
  std::cout << FORMAT(25) << "[main] Parameters D/nu: " << D << "/" << nu << std::endl;


  // Build the mesh and re-order vertices and edges to put the boundary ones first
  MeshBuilder builder = MeshBuilder(mesh_file);
  std::unique_ptr<Mesh> mesh_ptr = builder.build_the_mesh();
  std::cout << FORMAT(25) << "[main] Mesh size" << mesh_ptr->h_max() << std::endl;

  boost::timer::cpu_timer timer;
  
  // Create plates core
  timer.start();
  bool use_threads = (vm.count("pthread") ? vm["pthread"].as<bool>() : true);
  std::cout << "[main] " << (use_threads ? "Parallel execution" : "Sequential execution") << std:: endl;
  PlatesCore platescore(*mesh_ptr, K, use_threads);
  timer.stop();
  double t_wall_platescore, t_proc_platescore;
  std::tie(t_wall_platescore, t_proc_platescore) = store_times(timer, "[main] Time PlatesCore (wall/proc) ");

  // Assemble the problem
  KirchhoffLove::ConstitutiveLawType law = 
      [D, nu](const Eigen::Matrix2d tau)->Eigen::Matrix2d 
          { return 1./(D*(1.-nu)) * (tau-nu/(1.+nu)*tau.trace()*Eigen::Matrix2d::Identity()); };
  KirchhoffLove kl(platescore, law, use_threads);
  if(vm.count("stabilization-parameter")) {
    kl.stabilizationParameter() = vm["stabilization-parameter"].as<double>();
  }

  timer.start();
  std::cout << "[main] Assembling the problem" << std::endl;
  kl.assembleLinearSystem(f, u, grad_u);
  timer.stop();
  double t_wall_model, t_proc_model;
  std::tie(t_wall_model, t_proc_model) = store_times(timer, "[main] Time model (wall/proc) ");

  // Export matrix if requested  
  if (vm.count("export-matrix") && vm["export-matrix"].as<bool>()) {
    std::cout << "[main] Exporting matrix to Matrix Market format" << std::endl;
    saveMarket(kl.systemMatrix(), "A.mtx");
    saveMarket(kl.systemVector(), "b.mtx");
  }

  // Interpolate the exact solution
  std::cout << "[main] Interpolating the exact solution" << std::endl;
  Eigen::VectorXd uI = Eigen::VectorXd::Zero(kl.dimensionSpace());
  uI.head(kl.xDivDiv().dimension()) = kl.xDivDiv().interpolate(sigma, sigma_DE);
  uI.tail(kl.polykm2Th().dimension()) = kl.interpolateDeflection(u);

  // Solve the problem
  timer.start();
  Eigen::VectorXd uh_condensed;
  if (vm.count("iterative-solver")) {
    std::cout << "[main] Solving the linear system using BiCGSTAB" << std::endl;
    
    Eigen::BiCGSTAB<KirchhoffLove::SystemMatrixType, Eigen::IncompleteLUT<double> > solver;
    // solver.preconditioner().setFillfactor(2);
    solver.compute(kl.systemMatrix());
    if (solver.info() != Eigen::Success) {
      std::cerr << "[main] ERROR: Could not factorize matrix" << std::endl;
      exit(1);
    }
    uh_condensed = solver.solve(kl.systemVector());
    if (solver.info() != Eigen::Success) {
      std::cerr << "[main] ERROR: Could not solve the system" << std::endl;
      exit(1);
    }
  } else { 
#ifdef WITH_MKL
    std::cout << "[main] Solving the linear system using Pardiso" << std::endl;    
    Eigen::PardisoLU<KirchhoffLove::SystemMatrixType> solver;
#elif WITH_UMFPACK
    std::cout << "[main] Solving the linear system using Umfpack" << std::endl;    
    Eigen::UmfPackLU<KirchhoffLove::SystemMatrixType> solver;
#else
    std::cout << "[main] Solving the linear system using direct solver" << std::endl;    
    Eigen::SparseLU<KirchhoffLove::SystemMatrixType> solver;
#endif
    solver.compute(kl.systemMatrix());
    if (solver.info() != Eigen::Success) {
      std::cerr << "[main] ERROR: Could not factorize matrix" << std::endl;
    }
    uh_condensed = solver.solve(kl.systemVector());
    if (solver.info() != Eigen::Success) {
      std::cerr << "[main] ERROR: Could not solve linear system" << std::endl;
    }
  }
  // Re-create statically condensed unknowns
  Eigen::VectorXd uh = Eigen::VectorXd::Zero(kl.dimensionSpace());
  uh.head(kl.xDivDiv().dimension() - kl.nbSCDOFs()) = uh_condensed.head(kl.xDivDiv().dimension() - kl.nbSCDOFs()); 
  uh.segment(kl.xDivDiv().dimension() - kl.nbSCDOFs(), kl.nbSCDOFs()) = kl.scVector() + kl.scMatrix() * uh_condensed;
  uh.tail(kl.polykm2Th().dimension()) = uh_condensed.segment(kl.xDivDiv().dimension()-kl.nbSCDOFs(), kl.polykm2Th().dimension());

  timer.stop();
  
  double t_wall_solve, t_proc_solve;
  std::tie(t_wall_solve, t_proc_solve) = store_times(timer, "[main] Time solve (wall/proc) ");


  // Compute and print the energy error
  double error = kl.computeNorm(uh-uI);
  double norm_sol = kl.computeNorm(uI);
  std::cout << FORMAT(25) << "[main] Error" << error << std::endl;
  std::cout << FORMAT(25) << "[main] Rel. error " << error/norm_sol << std::endl;

  // --------------------------------------------------------------------------
  //                     Creates .txt file with data and results
  // --------------------------------------------------------------------------
  std::ofstream out("results.txt");
  out << "Solution: " << solution << std::endl;
  out << "Mesh: " << mesh_file << std::endl;
  out << "Degree: " << K << std::endl;
  out << "MeshSize: " << mesh_ptr->h_max() << std::endl;
  out << "NbCells: " << mesh_ptr->n_cells() << std::endl;
  out << "NbEdges: " << mesh_ptr->n_edges() << std::endl;
  out << "NbVertices: " << mesh_ptr->n_vertices() << std::endl;
  out << "DimXDivDiv: " << kl.xDivDiv().dimension() << std::endl;
  out << "DimPkm2: " << kl.polykm2Th().dimension() << std::endl;
  out << "DimSCsystem: " << kl.sizeSystem() << std::endl;
  out << "Error: " << error << std::endl;
  out << "RelError: " << error/norm_sol << std::endl;
  out << "Bending modulus (D): " << D << std::endl;
  out << "Poisson ratio (nu): " << nu << std::endl;
  out << "TwallPlatesCore: " << t_wall_platescore << std::endl;  
  out << "TprocPlatesCore: " << t_proc_platescore << std::endl;  
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
// Implementation of the KirchhoffLove class
//------------------------------------------------------------------------------

KirchhoffLove::KirchhoffLove(
                             const PlatesCore & platescore,
                             const ConstitutiveLawType & law,
                             bool use_threads,
                             std::ostream & output
                             )
  : m_platescore(platescore),
    m_law(law),
    m_use_threads(use_threads),
    m_output(output),
    m_xdivdiv(platescore, use_threads, output),
    m_Pkm2_Th(platescore.mesh(), 0, 0, PolynomialSpaceDimension<Cell>::Poly(platescore.degree() - 2)),
    m_nloc_sc_sigma(m_xdivdiv.numLocalDofsCell()),
    m_A(sizeSystem(), sizeSystem()),
    m_b(Eigen::VectorXd::Zero(sizeSystem())),
    m_sc_A(nbSCDOFs(), sizeSystem()),
    m_sc_b(Eigen::VectorXd::Zero(nbSCDOFs())),
    m_stab_par(1.)
{
  m_output << "[KirchhoffLove] Initializing" << std::endl;
  
  // To avoid performing static condensation, initialize m_nloc_sc_sigma at 0
}

//------------------------------------------------------------------------------

std::pair<Eigen::MatrixXd, Eigen::VectorXd>
KirchhoffLove::_compute_local_contribution(size_t iT,
                                           const ForcingTermType & f,
                                           const DeflectionType & u,
                                           const GradientDeflectionType & grad_u
                                           )
{
  const Cell & T = *m_platescore.mesh().cell(iT);

  size_t dim_xdivdiv_T = m_xdivdiv.dimensionCell(iT);
  size_t dim_Pkm2_T = m_Pkm2_Th.dimensionCell(iT);
  size_t dim_T = dim_xdivdiv_T + dim_Pkm2_T;

  // Local left-hand side matrix
  Eigen::MatrixXd AT = Eigen::MatrixXd::Zero(dim_T, dim_T);

  AT.topLeftCorner(dim_xdivdiv_T, dim_xdivdiv_T)
    = m_xdivdiv.computeL2Product(iT, m_law, m_stab_par);
  AT.topRightCorner(dim_xdivdiv_T, dim_Pkm2_T)
    = m_xdivdiv.cellOperators(iT).divdiv_rhs.transpose();
  AT.bottomLeftCorner(dim_Pkm2_T, dim_xdivdiv_T)
    = -m_xdivdiv.cellOperators(iT).divdiv_rhs;

  // Local right-hand side vector: loading term
  Eigen::VectorXd bT = Eigen::VectorXd::Zero(dim_T);
  QuadratureRule quad_2kp2_T = generate_quadrature_rule(T, 2 * (m_platescore.degree() + 1));
  bT.tail(dim_Pkm2_T)
    = integrate(f, evaluate_quad<Function>::compute(*m_platescore.cellBases(iT).Polykm2, quad_2kp2_T), quad_2kp2_T);

  // Local right-hand side vector: boundary conditions
  if (T.is_boundary()){
    for (size_t iE=0; iE < T.n_edges(); iE++){
      Edge & E = *T.edge(iE);
      
      double omegaTE = T.edge_orientation(iE);
      VectorRd nE = E.normal();
      VectorRd tE = E.tangent();
      
      // Edge contributions if boundary
      if (E.is_boundary()){
          
        QuadratureRule quad_E = generate_quadrature_rule(E, 2*xDivDiv().degree()+6);

        // Term (PE tau) nE.nE partial_nE u
        FType<double> dnE_u = [&grad_u, &nE](const VectorRd & x)->double { return grad_u(x).dot(nE);};
        bT.head(dim_xdivdiv_T)
          -= omegaTE * xDivDiv().extendOperator(T, E, xDivDiv().edgePotential(E)).transpose() 
          * integrate(dnE_u, evaluate_quad<Function>::compute(*xDivDiv().edgeBases(E).Polykm1, quad_E), quad_E);

        // Term D_{tau,E} u
        bT.segment(xDivDiv().localOffset(T,E) + xDivDiv().edgeBases(E).Polykm3->dimension(), xDivDiv().edgeBases(E).Polykm2->dimension())
          += omegaTE * integrate(u, evaluate_quad<Function>::compute(*xDivDiv().edgeBases(E).Polykm2, quad_E), quad_E);

        // Vertex contributions (note that contributions of boundary vertices for internal edges is zero
        //  by compensation, and thus not counted here)
        for (size_t iV=0; iV<2; iV++){
          const Vertex & V = *E.vertex(iV);
          double omegaEV = E.vertex_orientation(iV);
          double uV = u(V.coords());
          for (size_t j = 0; j < xDivDiv().SymBasisVertex().size(); j++) {
            bT(xDivDiv().localOffset(T, V)+j) -= omegaTE * omegaEV * (xDivDiv().SymBasisVertex().function(j) * nE).dot(tE) * uV;
          } // for j
        } // for iV
      } // E.is_boundary()
    } // for iE
    
  } // T.is_boundary()
  
  return std::make_pair(AT, bT);
}

//------------------------------------------------------------------------------

LocalStaticCondensation KirchhoffLove::_compute_static_condensation(const size_t & iT) const
{
  const Cell & T = *m_platescore.mesh().cell(iT);
  
  // Dimensions
  size_t dim_sigma = m_xdivdiv.dimensionCell(iT) - m_nloc_sc_sigma;    // dimension of skeletal unknowns for sigma
  size_t dim_u = m_Pkm2_Th.dimensionCell(iT);         // dimension of u
  size_t dim_dofs = dim_sigma+dim_u;      // nb of dofs in system, remaining after SC

  // Creation of permutation matrix
  Eigen::MatrixXd Perm = Eigen::MatrixXd::Zero(dim_dofs+m_nloc_sc_sigma, dim_dofs+m_nloc_sc_sigma);
  Perm.topLeftCorner(dim_sigma, dim_sigma) = Eigen::MatrixXd::Identity(dim_sigma, dim_sigma);
  Perm.block(dim_sigma+dim_u, dim_sigma, m_nloc_sc_sigma, m_nloc_sc_sigma) = Eigen::MatrixXd::Identity(m_nloc_sc_sigma, m_nloc_sc_sigma);
  Perm.block(dim_sigma, dim_sigma+m_nloc_sc_sigma, dim_u, dim_u) = Eigen::MatrixXd::Identity(dim_u, dim_u);

  // Creation of global DOFs for system
  std::vector<size_t> IT_sys(dim_dofs, 0);
  auto IT_xdivdiv = m_xdivdiv.globalDOFIndices(T);
  auto IT_Pkm2 = m_Pkm2_Th.globalDOFIndices(T);
  // put skeletal DOFs of sigma in Idofs_T
  auto it_IT_sys = std::copy(IT_xdivdiv.begin(), IT_xdivdiv.begin()+dim_sigma, IT_sys.begin()); 
  size_t offset = m_xdivdiv.dimension() - nbSCDOFs();   // nb total of skeletal DOFs for sigma (where global dofs of u will start)
  std::transform(IT_Pkm2.begin(), IT_Pkm2.end(), it_IT_sys, [&offset](const size_t & index) { return index + offset; });

  // Creation of global DOFs for SC operator      
  std::vector<size_t> IT_sc(m_nloc_sc_sigma, 0);
  std::transform(IT_xdivdiv.begin()+dim_sigma, IT_xdivdiv.end(), IT_sc.begin(), [&offset](const size_t & index) { return index - offset; });
  
  return LocalStaticCondensation(Perm, IT_sys, IT_sc);
}

//------------------------------------------------------------------------------

void KirchhoffLove::_assemble_local_contribution( 
                                                 size_t iT,
                                                 const std::pair<Eigen::MatrixXd, Eigen::VectorXd> & lsT,
                                                 std::list<Eigen::Triplet<double> > & triplets_sys,
                                                 Eigen::VectorXd & rhs_sys,
                                                 std::list<Eigen::Triplet<double> > & triplets_sc,
                                                 Eigen::VectorXd & rhs_sc
                                                 )
{
  // Get information for local static condensation
  LocalStaticCondensation locSC = _compute_static_condensation(iT);

  Eigen::MatrixXd AT_sys, AT_sc;
  Eigen::VectorXd bT_sys, bT_sc;
  std::tie(AT_sys, bT_sys, AT_sc, bT_sc) = locSC.compute(lsT);
  
  // STATICALLY CONDENSED SYSTEM
  std::vector<size_t> IT_sys = locSC.globalDOFs_sys();
  for (size_t i = 0; i < locSC.dim_sys(); i++){
    rhs_sys(IT_sys[i]) += bT_sys(i);
    for (size_t j = 0; j < locSC.dim_sys(); j++){    
      triplets_sys.emplace_back(IT_sys[i], IT_sys[j], AT_sys(i,j));
    }
  }

  // RECOVERY OPERATOR
  std::vector<size_t> IT_sc = locSC.globalDOFs_sc();
  for (size_t i = 0; i < locSC.dim_sc(); i++){
    rhs_sc(IT_sc[i]) += bT_sc(i);
    for (size_t j = 0; j < locSC.dim_sys(); j++){
      triplets_sc.emplace_back(IT_sc[i], IT_sys[j], AT_sc(i,j));
    }
  }

}

//------------------------------------------------------------------------------

void KirchhoffLove::assembleLinearSystem(
                                         const ForcingTermType & f,
                                         const DeflectionType & u,
                                         const GradientDeflectionType & grad_u
                                         )
{
  auto assemble_all = [this, f, u, grad_u](
                                   size_t start,
                                   size_t end,
                                   std::list<Eigen::Triplet<double> > * triplets_sys,
                                   Eigen::VectorXd * rhs_sys,
                                   std::list<Eigen::Triplet<double> > * triplets_sc,
                                   Eigen::VectorXd * rhs_sc
                                   )->void
  {
    for (size_t iT = start; iT < end; iT++) {
      this->_assemble_local_contribution(
                                         iT,
                                         this->_compute_local_contribution(iT, f, u, grad_u),
                                         *triplets_sys,
                                         *rhs_sys,
                                         *triplets_sc,
                                         *rhs_sc
                                         );
    } // for iT    
  };
  if (m_use_threads) {
    m_output << "[KirchhoffLove] Parallel assembly" << std::endl;
  } else {
    m_output << "[KirchhoffLove] Sequential assembly" << std::endl;
  }
  std::tie(m_A, m_b, m_sc_A, m_sc_b) = parallel_assembly_system(m_platescore.mesh().n_cells(), this->sizeSystem(), std::make_pair(this->nbSCDOFs(), this->sizeSystem()), this->nbSCDOFs(), assemble_all, m_use_threads);

}

//------------------------------------------------------------------------------

Eigen::VectorXd KirchhoffLove::interpolateDeflection(
                                                     const DeflectionType & u,
                                                     const int deg_quad
                                                     ) const
{
  Eigen::VectorXd uI = Eigen::VectorXd::Zero(m_Pkm2_Th.dimension());

  size_t dqr = (deg_quad >= 0 ? deg_quad : 2 * (m_platescore.degree() + 1));

  auto interpolate_all = [this, &uI, u, &dqr](size_t start, size_t end)->void
  {
    for (size_t iT = start; iT < end; iT++) {
      const Cell & T = *m_platescore.mesh().cell(iT);
      
      QuadratureRule quad_dqr_T = generate_quadrature_rule(T, dqr);
      
      auto basis_Pkm2_quad
        = evaluate_quad<Function>::compute(*m_platescore.cellBases(iT).Polykm2, quad_dqr_T);

      uI.segment(m_Pkm2_Th.globalOffset(T), m_Pkm2_Th.dimensionCell(iT))
        = l2_projection(u, *m_platescore.cellBases(iT).Polykm2, quad_dqr_T, basis_Pkm2_quad);
    } // for iT
  };
  parallel_for(m_platescore.mesh().n_cells(), interpolate_all, m_use_threads);

  return uI;
}

//------------------------------------------------------------------------------
double KirchhoffLove::computeNorm(const Eigen::VectorXd & v ) const
{
  Eigen::VectorXd local_sqnorms = Eigen::VectorXd::Zero(m_platescore.mesh().n_cells());

  // xDivDiv correspond to the first components of v, Pkm2_Th to the last ones
  Eigen::VectorXd sigmaI = v.head(xDivDiv().dimension());
  Eigen::VectorXd vI = v.tail(polykm2Th().dimension());
  
  std::function<void(size_t, size_t)> compute_local_squarednorms
    = [this, &sigmaI, &vI, &local_sqnorms](size_t start, size_t end)->void
    {
      for (size_t iT = start; iT < end; iT++){
        Cell & T = *m_platescore.mesh().cell(iT);
        Eigen::VectorXd sigmaT = m_xdivdiv.restrict(T, sigmaI);
        Eigen::VectorXd vT = m_Pkm2_Th.restrict(T, vI);

        // Contribution of L2 norms
        local_sqnorms(iT) += sigmaT.transpose() * m_xdivdiv.computeL2Product(iT, m_law, m_stab_par) * sigmaT;
        local_sqnorms(iT) += vT.transpose() * GramMatrix(T, *xDivDiv().cellBases(iT).Polykm2) * vT;

      }
    };
  parallel_for(m_platescore.mesh().n_cells(), compute_local_squarednorms, m_use_threads);
  
  return std::sqrt(std::abs(local_sqnorms.sum()));
}
