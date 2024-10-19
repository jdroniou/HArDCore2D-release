// Author: Jerome Droniou (jerome.droniou@monash.edu)
#include <fstream>
#include <iomanip>
#include <thread>

#include "hho-fullgradientdiff.hpp"
#include <parallel_for.hpp>
#include <GMpoly_cell.hpp>

#include <boost/program_options.hpp>
#include <display_timer.hpp>

#include "vtu_writer.hpp"

#ifdef WITH_UMFPACK
#include <Eigen/UmfPackSupport>
#endif

#ifdef WITH_MKL
#include <Eigen/PardisoSupport>
#include <mkl.h>
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
constexpr size_t max_nb_cells_for_plot = 20000;

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
    ("solution,s", boost::program_options::value<int>()->default_value(0), "Select the solution")
    ("export-matrix,e", "Export matrix to Matrix Market format")
    ("plot", boost::program_options::value<std::string>(), "Save plot of the solution to the given filename")
    ("solver", boost::program_options::value<std::string>()->default_value("PardisoLU"), "Choice of solver, not case dependent. Options are: PardisoLU, UMFPACK, PaStiXLU, PaStiXLLT, EigenLU, EigenBiCGSTAB (reverts to EigenLU if the selected solver is not available)")
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

  std::cout << "[main] Mesh file: " << mesh_file << std::endl;
  
  // Select the degree 
  size_t K = vm["degree"].as<size_t>();
  std::cout << FORMAT(25) << "[main] Degree" << K << std::endl;

  // Select the solution
  int solution = (vm.count("solution") ? vm["solution"].as<int>() : 0);
  FullGradientDiffusion::ForcingTermType f;
  FullGradientDiffusion::PermeabilityType kappa(1.);    // Default permeability constant equal to 1.
  FullGradientDiffusion::SolutionType u;
  FullGradientDiffusion::SolutionGradientType gradu;

  switch (solution) {
  case 0:
    std::cout << "[main] Linear solution" << std::endl;    
    f = linear_f;
    u = linear_u;
    kappa = linear_kappa;
    gradu = linear_gradu;
    break;

  case 1:
    std::cout << "[main] Trigonometric solution" << std::endl;    
    f = trigonometric_f;
    u = trigonometric_u;
    kappa = trigonometric_kappa;
    gradu = trigonometric_gradu;
    break;

  default:
    std::cerr << "[main] ERROR: Unknown exact solution" << std::endl;
    exit(1);
  }

  // Build the mesh and reorder edges to handle Dirichlet BCs (Dirichlet edges at the end -- pure Dirichlet for the moment)
  MeshBuilder builder = MeshBuilder(mesh_file);
  std::unique_ptr<Mesh> mesh_ptr = builder.build_the_mesh();
  std::cout << FORMAT(25) << "[main] Mesh size" << mesh_ptr->h_max() << std::endl;
  BoundaryConditions BC("D", *mesh_ptr.get());
  BC.reorder_edges("start"); 

  boost::timer::cpu_timer timer;
  // Create HHO space
  timer.start();
  bool use_threads = (vm.count("pthread") ? vm["pthread"].as<bool>() : true);
  std::cout << "[main] " << (use_threads ? "Parallel execution" : "Sequential execution") << std:: endl;
  HHOSpace hho_space(*mesh_ptr, K, use_threads);
  timer.stop();
  double t_wall_hhospace, t_proc_hhospace;
  std::tie(t_wall_hhospace, t_proc_hhospace) = store_times(timer, "[main] Time HHOSpace (wall/proc) ");

  // Assemble the problem
  timer.start();
  FullGradientDiffusion diff(hho_space, BC, use_threads);
  if(vm.count("stabilization-parameter")) {
    diff.stabilizationParameter() = vm["stabilization-parameter"].as<double>();
  }
  Eigen::VectorXd UDir = Eigen::VectorXd::Zero(diff.numDirDOFs());
  diff.assembleLinearSystem(f, kappa, u, UDir);
  timer.stop();
  double t_wall_model, t_proc_model;
  std::tie(t_wall_model, t_proc_model) = store_times(timer, "[main] Time model (wall/proc) ");

  // Export matrix if requested  
  if (vm.count("export-matrix")) {
    std::cout << "[main] Exporting matrix to Matrix Market format" << std::endl;
    saveMarket(diff.systemMatrix(), "A_fullgradientdiff.mtx");
    saveMarket(diff.systemVector(), "b_fullgradientdiff.mtx");
  }

  // Select linear solver
  std::string name_solver = vm["solver"].as<std::string>();
  LinearSolver<FullGradientDiffusion::SystemMatrixType> solver(name_solver);
  std::cout << "[main] Solving the system using " << solver.name() << std::endl;

  // Solve the problem
  timer.start();
  Eigen::VectorXd uh_solsystem = solver.compute_and_solve(diff.systemMatrix(), diff.systemVector());

  // Re-create boundary values and statically condensed unknowns
  Eigen::VectorXd uh = Eigen::VectorXd::Zero(diff.hhospace().dimension());
  uh.head(diff.numDirDOFs()) = UDir;
  uh.segment(diff.numDirDOFs(), diff.sizeSystem()) = uh_solsystem;
  uh.tail(diff.numSCDOFs()) = diff.scVector() + diff.scMatrix() * uh.head(diff.numSkeletalDOFs());

  timer.stop();
  double t_wall_solve, t_proc_solve;
  std::tie(t_wall_solve, t_proc_solve) = store_times(timer, "[main] Time solve (wall/proc) ");
    
  // Compute the error in the L2, H1 and energy norms
  Eigen::VectorXd uI = diff.hhospace().interpolate(u);
  Eigen::VectorXd eh = uh - uI;
  
  std::vector<std::pair<double,double>> list_L2_norms = diff.hhospace().computeNorms(std::vector<Eigen::VectorXd> {eh, uI});
  double L2_err = list_L2_norms[0].first / list_L2_norms[1].first;
  std::cout << "[main] Discrete L2 error " << L2_err << std::endl;
  double H1_err = list_L2_norms[0].second / list_L2_norms[1].second;
  std::cout << "[main] Discrete H1 error " << H1_err << std::endl;

  std::vector<double> list_energy_norms = diff.computeEnergyNorms(std::vector<Eigen::VectorXd> {eh, uI});
  double energy_err = list_energy_norms[0] / list_energy_norms[1];
  std::cout << "[main] Discrete energy error " << energy_err << std::endl;
  std::cout << "[main] Mesh diameter " << mesh_ptr->h_max() << std::endl;
 

  // --------------------------------------------------------------------------
  //                     Creates a .vtu file of the solution
  // --------------------------------------------------------------------------
  // Only if we do not have too many cells
  if (vm.count("plot") && mesh_ptr->n_cells() <= max_nb_cells_for_plot) {
    std::cout << "[main] Writing solution to file" << std::endl;
    std::string filename = vm["plot"].as<std::string>();
    VtuWriter plotdata(mesh_ptr.get());
		
		// Exact solution at the vertices
    Eigen::VectorXd exact_u_vertex(mesh_ptr->n_vertices());
    for (size_t iV=0; iV< mesh_ptr->n_vertices(); iV++){
      exact_u_vertex[iV] = u(mesh_ptr->vertex(iV)->coords());
    }
		
    // Approximate solution at the vertices
    Eigen::VectorXd u_vertex = hho_space.computeVertexValues(uh);
    
    // Plot
    plotdata.write_to_vtu(filename + ".vtu", u_vertex, "approximate");
    plotdata.write_to_vtu("exact-" + filename + ".vtu", exact_u_vertex, "exact");
   }

  // --------------------------------------------------------------------------
  //                     Creates .txt file with data and results
  // --------------------------------------------------------------------------
  std::ofstream out("results.txt");
  out << "Scheme: hho-fullgradientdiff" << std::endl;
  out << "Solution: " << solution << std::endl;
  out << "Mesh: " << mesh_file << std::endl;
  out << "Degree: " << K << std::endl;
  out << "StabilisationParameter: " << diff.stabilizationParameter() << std::endl;
  out << "MeshSize: " << mesh_ptr->h_max() << std::endl;
  out << "NbCells: " << mesh_ptr->n_cells() << std::endl;
  out << "NbEdges: " << mesh_ptr->n_edges() << std::endl;
  out << "DimHHOSpace: " << diff.hhospace().dimension() << std::endl;
  out << "SizeSystem: " << diff.sizeSystem() << std::endl;
  out << "L2Error: " << L2_err << std::endl;  
  out << "H1Error: " << H1_err << std::endl;  
  out << "EnergyError: " << energy_err << std::endl;  
  out << "TwallHHOSpace: " << t_wall_hhospace << std::endl;  
  out << "TprocHHOSpace: " << t_proc_hhospace << std::endl;  
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
// Diffusion model
//------------------------------------------------------------------------------

FullGradientDiffusion::FullGradientDiffusion(
                               const HHOSpace & hho_space,
                               const BoundaryConditions & BC,
                               bool use_threads,
                               std::ostream & output
                               )
  : m_hhospace(hho_space),
    m_BC(BC),
    m_use_threads(use_threads),
    m_output(output),
    m_nloc_sc(m_hhospace.numLocalDofsCell()),
    m_A(sizeSystem(), sizeSystem()),
    m_b(Eigen::VectorXd::Zero(sizeSystem())),
    m_sc_A(numSCDOFs(), sizeSystem()),
    m_sc_b(Eigen::VectorXd::Zero(numSCDOFs())),
    m_stab_par(1.)
    
{
  m_output << "[FullGradientDiffusion] Initializing" << std::endl;
  
  // To avoid performing static condensation, initialise m_nloc_sc at 0
}

//------------------------------------------------------------------------------

void FullGradientDiffusion::assembleLinearSystem(
                                          const ForcingTermType & f,
                                          const PermeabilityType & kappa,
                                          const SolutionType & u,
                                          Eigen::VectorXd & UDir
                                          )
{
  // Assemble all local contributions
  auto assemble_all = [this, f, kappa, u](
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
                                                             this->_compute_local_contribution(iT, f, kappa),
                                                             *triplets_sys,
                                                             *rhs_sys,
                                                             *triplets_sc,
                                                             *rhs_sc
                                                             );
                        } // for iT
                      };
                      
  // Assemble the matrix and rhs
  if (m_use_threads) {
    m_output << "[FullGradientDiffusion] Parallel assembly" << std::endl;
  }else{
    m_output << "[FullGradientDiffusion] Sequential assembly" << std::endl;
  }
  FullGradientDiffusion::SystemMatrixType global_matrix;  // Matrix without BCs (but after static condensation)
  Eigen::VectorXd global_rhs;     // RHS without BCs (but after static condensation)
  std::tie(global_matrix, global_rhs, m_sc_A, m_sc_b) = parallel_assembly_system(m_hhospace.mesh().n_cells(), this->numSkeletalDOFs(), std::make_pair(this->numSCDOFs(), this->numSkeletalDOFs()), this->numSCDOFs(), assemble_all, m_use_threads);
    
  // System matrix without Dirichlet boundary dofs
  m_A = global_matrix.bottomRightCorner(sizeSystem(), sizeSystem());

  // Adjust RHS for boundary conditions
  // Creation of Dirichlet vector: L2 projection of solution only on Dirichlet edges (reordered at the start)
  for (size_t idE = 0; idE < m_BC.n_dir_edges(); idE++){
    Edge* E = mesh().edge(idE);
    QuadratureRule quadE = generate_quadrature_rule(*E, 2*m_hhospace.degree()+2);
    
    auto phiE_quadE = evaluate_quad<Function>::compute(*m_hhospace.edgeBases(idE).Polyk, quadE);
    UDir.segment(idE * m_hhospace.numLocalDofsEdge(), m_hhospace.numLocalDofsEdge()) 
        = l2_projection<HHOSpace::PolyBasisEdgeType>(u, *m_hhospace.edgeBases(idE).Polyk, quadE, phiE_quadE);
  } // for idE

  m_b = global_rhs.tail(sizeSystem()) - global_matrix.bottomLeftCorner(sizeSystem(), numDirDOFs()) * UDir;
  
}

//------------------------------------------------------------------------------
  
std::pair<Eigen::MatrixXd, Eigen::VectorXd>
FullGradientDiffusion::_compute_local_contribution(
                                            size_t iT, const ForcingTermType & f,
                                            const PermeabilityType & kappa
                                            )
{
  const Cell & T = *m_hhospace.mesh().cell(iT);

  size_t dim_T = m_hhospace.dimensionCell(iT);
  Eigen::MatrixXd AT = Eigen::MatrixXd::Zero(dim_T, dim_T);
  Eigen::VectorXd lT = Eigen::VectorXd::Zero(dim_T);

  //------------------------------------------------------------------------------
  // Local matrix
  //------------------------------------------------------------------------------

////// kappa constant scalar in element for the moment

  // Mass matrix for (P^k(T))^d
  MonomialCellIntegralsType int_mono_2k = IntegrateCellMonomials(T, 2*m_hhospace.degree());  
  Eigen::MatrixXd mass_Pkd_T = GramMatrix(T, *m_hhospace.cellBases(iT).Polykd, int_mono_2k);

  double kappaT = kappa.value(T, T.center_mass());
  AT += kappaT * m_hhospace.operators(iT).gradient.transpose() * mass_Pkd_T * m_hhospace.operators(iT).gradient
       + kappaT * m_stab_par * m_hhospace.operators(iT).stabilisation;
  
  //------------------------------------------------------------------------------
  // Local vector
  //------------------------------------------------------------------------------

  QuadratureRule quad_2k_T = generate_quadrature_rule(T, 2 * m_hhospace.degree());
  lT.tail(m_hhospace.numLocalDofsCell()) = 
      integrate(f, evaluate_quad<Function>::compute(*m_hhospace.cellBases(iT).Polyk, quad_2k_T), quad_2k_T);
  
  return std::make_pair(AT, lT);
}

//------------------------------------------------------------------------------

LocalStaticCondensation FullGradientDiffusion::_compute_static_condensation(const size_t & iT) const
{
  const Cell & T = *mesh().cell(iT);

  // Dimensions
  size_t dim_dofs = m_hhospace.dimensionCell(iT) - m_nloc_sc;    // dimension of skeletal unknowns (not statically condensed)

  // Creation of permutation matrix: trivial since cell DOFs are already at the end
  Eigen::MatrixXd Perm = Eigen::MatrixXd::Identity(dim_dofs+m_nloc_sc, dim_dofs+m_nloc_sc);

  // Creation of global DOFs for system: IT_sys contains the skeletal dofs
  std::vector<size_t> IT_sys(dim_dofs,0);
  auto IT_hhospace = m_hhospace.globalDOFIndices(T);
  std::copy(IT_hhospace.begin(), IT_hhospace.begin()+dim_dofs, IT_sys.begin());
   
  // Creation of global DOFs for SC operator: IT_sc contains global cell dofs (offset to start at 0)
  std::vector<size_t> IT_sc(m_nloc_sc);
  size_t offset = m_hhospace.dimension() - numSCDOFs();     // offset where the SC unknowns start in the global system
  std::transform(IT_hhospace.begin()+dim_dofs, IT_hhospace.end(), IT_sc.begin(), [&offset](const size_t & index) { return index - offset; });
 
  return LocalStaticCondensation(Perm, IT_sys, IT_sc);
}

//------------------------------------------------------------------------------

void FullGradientDiffusion::_assemble_local_contribution(
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

std::vector<double> FullGradientDiffusion::computeEnergyNorms( const std::vector<Eigen::VectorXd> & list_dofs ) const
{
  size_t nb_vectors = list_dofs.size();
  std::vector<Eigen::VectorXd> local_sqnorms(nb_vectors, Eigen::VectorXd::Zero(m_hhospace.mesh().n_cells()));

  std::function<void(size_t, size_t)> compute_local_squarednorms
    = [this, &list_dofs, &local_sqnorms, &nb_vectors](size_t start, size_t end)->void
    {
      for (size_t iT = start; iT < end; iT++){
        Cell & T = *m_hhospace.mesh().cell(iT);

        // Mass matrix for (P^k(T))^d
	      MonomialCellIntegralsType int_mono_2k = IntegrateCellMonomials(T, 2*m_hhospace.degree());
        Eigen::MatrixXd mass_Pkd_T = GramMatrix(T, *m_hhospace.cellBases(iT).Polykd, int_mono_2k);
        Eigen::MatrixXd GT = m_hhospace.operators(iT).gradient;
        Eigen::MatrixXd ST = m_hhospace.operators(iT).stabilisation;

        for (size_t i=0; i<nb_vectors; i++){
          Eigen::VectorXd uT = m_hhospace.restrict(T, list_dofs[i]);

          // Contribution of gradients, without any weight (no permeability)
          local_sqnorms[i](iT) += uT.transpose() * (GT.transpose() * mass_Pkd_T * GT + ST) * uT;
        }
      }
    };
  parallel_for(m_hhospace.mesh().n_cells(), compute_local_squarednorms, m_use_threads);
  
  // Vector of outputs
  std::vector<double> list_norms(nb_vectors, 0.);
  for (size_t i=0; i<nb_vectors; i++){
    list_norms[i] = std::sqrt(std::abs(local_sqnorms[i].sum()));
  }
  
  return list_norms;
}


