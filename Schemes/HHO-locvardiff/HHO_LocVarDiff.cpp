// Implementation of the HHO scheme in 2D for the diffusion equation
//
//   { -div(K \grad(u)) = f,       inside Omega
//   { K \grad(u) . nTF = g,       on GammaN
//   { 								u = g,			 on GammaD
// 
//
// Author: Jerome Droniou (jerome.droniou@monash.edu)
//

#include <iostream>
#include <fstream>

#include <boost/program_options.hpp>
#include <boost/timer/timer.hpp>
#include <unsupported/Eigen/SparseExtra>    // For saveMarket function

#include "mesh.hpp"
#include "import_mesh.hpp"
#include "mesh_builder.hpp"

#include "HHO_LocVarDiff.hpp"
#include "vtu_writer.hpp"

using namespace HArDCore2D;

// Mesh filenames
const std::string mesh_dir = "../../typ2_meshes/";
std::string default_mesh = mesh_dir + "cart5x5.typ2";
const std::string default_solver_type = "bicgstab";

// Max number of cells to plot a graph
constexpr size_t max_nb_cells_for_plot = 15000;

int main(int argc, const char* argv[]) {

  // --------------------------------------------------------------------------
  //                          Program options
  // --------------------------------------------------------------------------
  boost::program_options::options_description desc("Allowed options");
  desc.add_options()
    ("help,h", "Display this help message")
    ("mesh,m", boost::program_options::value<std::string>(), "Set the mesh")
      ("edgedegree,k", boost::program_options::value<size_t>()->default_value(1), "The polynomial degree on the edges")
    ("celldegree,l", boost::program_options::value<size_t>()->default_value(1), "The polynomial degree in the cells")
    ("bc_id,b", boost::program_options::value<std::string>()->default_value("D"), "Set the boundary conditions (D=Dirichlet, N=Neumann, Mx=Mixed number x)")
    ("testcase,c", boost::program_options::value<std::vector<int>>()->multitoken(), "Set the test case (as '-c n m'; n=exact sol, m=diffusion)")
    ("plot,p", boost::program_options::value<std::string>()->default_value("solution"), "Save plot of the solution to the given filename")
    ("solver_type", boost::program_options::value<std::string>()->default_value("bicgstab"), "Defines the linear solver for the system")
    ("use_threads", boost::program_options::value<bool>()->default_value(true), "Using multithreading")
    ("export_matrix,e", boost::program_options::value<bool>()->default_value(false), "Export matrix to Matrix Market format");

  boost::program_options::variables_map vm;
  boost::program_options::store(boost::program_options::parse_command_line(argc, argv, desc), vm);
  boost::program_options::notify(vm);

  std::ostream & output = std::cout;

  // Display the help options
  if (vm.count("help")) {
    output << desc << std::endl;
    return 0;
  }

  // Select the mesh
  std::string mesh_file = (vm.count("mesh") ? vm["mesh"].as<std::string>() : default_mesh);

  // Select the degree
  size_t K = vm["edgedegree"].as<size_t>();
  int L = vm["celldegree"].as<size_t>();

  // Select solver type
  std::string solver_type = vm["solver_type"].as<std::string>();

  // Select multithreading or not
  bool use_threads = vm["use_threads"].as<bool>();

  // Checks
  if ( (abs(int(K)-L) > 1) || (K<0) || (L<-1) ){
    output << "Degrees k and l not in acceptable range (k positive and l=k-1, k, or k+1): k=" << K << ", l=" << L << "\n\n";
    output << std::abs(int(K)-L);
    exit(1);
  }


  // --------------------------------------------------------------------------
  //                        Create the HHO data structure
  // --------------------------------------------------------------------------

  // Build mesh
  MeshBuilder builder = MeshBuilder(mesh_file);
  std::unique_ptr<Mesh> mesh_ptr = builder.build_the_mesh();

  // Get the BC and re-order the edges
  std::string bc_id = vm["bc_id"].as<std::string>();
  BoundaryConditions BC(bc_id, *mesh_ptr.get());
  BC.reorder_edges();
 
  // Create the HHO structure: need basis in Poly{K+1} on cells and Poly{K} on faces
  HybridCore hho(mesh_ptr.get(), K+1, K, use_threads, output);

  // --------------------------------------------------------------------------
  //                        Create the model equation
  // --------------------------------------------------------------------------
	
  const std::vector<int> default_id_tcase = std::vector<int>{1,1};
  std::vector<int> id_tcase = (vm.count("testcase") ? vm["testcase"].as<std::vector<int>>() : default_id_tcase);
  TestCase tcase(id_tcase);

  // Diffusion tensor, source term, solution and gradient
  CellFType<MatrixRd> kappa = tcase.diff();
  size_t deg_kappa = tcase.get_deg_diff();
  CellFType<double> source = tcase.diff_source();
  FType<double> exact_solution = tcase.sol();
  CellFType<VectorRd> grad_exact_solution = tcase.grad_sol();

  // Create the model equation
  HHO_LocVarDiff model(hho, K, L, kappa, deg_kappa, source, BC, exact_solution, grad_exact_solution, solver_type, use_threads, output);

  // --------------------------------------------------------------------------
  // 											Recalling the mesh and parameters
  // --------------------------------------------------------------------------
  output << "[Scheme] Data:\n";
  output << " Boundary conditions: " << BC.name() << " /";
  output << " Test case: solution = " << id_tcase[0] << "; diffusion = " << id_tcase[1] << "\n";
  size_t found = mesh_file.find_last_of("/\\");
  std::string mesh_name = mesh_file.substr(found+1);
  std::vector<double> meshreg = mesh_ptr->regularity();
  std::cout << "Mesh = " << mesh_name << " (nb cells= " << mesh_ptr->n_cells() << ", nb edges= " << mesh_ptr->n_edges() << ", reg= " << meshreg[0] << ", " << meshreg[1] << ")\n";
  size_t nbedgedofs = model.get_ntotal_edge_dofs();
  output << "  Degrees: edge = " << K << "; cell = " << L << " [Nb edge DOFs = " << nbedgedofs << "]\n";
  output << "  Using threads: " << (use_threads ? "true" : "false") << "\n";

  // --------------------------------------------------------------------------
  //                        Assemble and solve the model problem
  // --------------------------------------------------------------------------

//  boost::timer::cpu_timer timer;
  output << "\n[Scheme] Assembling." << std::endl;
  model.assemble(hho);
  output << "  Assembly time = " << model.get_assembly_time() << "s\n";

  output << "\n[Scheme] Solving." << std::endl;
  UVector u = model.solve(hho);
  // To export the matrix (can be read in Octave/Matlab using script mmread)
  bool export_matrix = (vm.count("export_matrix") ? vm["export_matrix"].as<bool>() : false);
  if (export_matrix) {
    output << "  Exporting matrix to Matrix Market format\n" << std::endl;
    saveMarket(model.get_SysMat(), "SystemMatrix.mtx");
  }

  output << "  Solver= " << solver_type << "\n";
  output << "  Solving time = " << model.get_solving_time() << "s\n";
  output << "  Residual of the linear system = " << model.get_solving_error() << '\n';

  printf("interm times = %f | %f | %f | %f | %f | %f \n", 
	 model.get_itime(0),  model.get_itime(1),  model.get_itime(2),  model.get_itime(3), model.get_itime(4), model.get_itime(5));
  //	printf("interm times = %f | %f | %f | %f | %f \n", 
  //			model.get_itime(5), model.get_itime(6), model.get_itime(7), model.get_itime(8), model.get_itime(9));

  // --------------------------------------------------------------------------
  //                        Compute the error
  // --------------------------------------------------------------------------

  output << "\n[Scheme] Computing error." << std::endl;

  // Interpolate the exact solution
  //timer.start();
  UVector Uh = hho.interpolate(exact_solution, L, K, 2*K+2);
  output << "  Interpolant calculated"<<std::endl;
  //output << "Interpolant of exact solution computed in " << double(timer.elapsed().wall) * pow(10, -9) << " seconds.\n";

  // Compute the L2 error
  //timer.start();
  double L2error = hho.L2norm(u - Uh)/hho.L2norm(Uh);
  //output << "L2 Error = " << L2error << ", computed in " << double(timer.elapsed().wall) * pow(10, -9) << " seconds.\n";
  output << "  L2 Error = " << L2error << "\n";

  // Compute the error in discrete H1 norm and energy norm
  //timer.start();
  double H1error = hho.H1norm(u - Uh)/hho.H1norm(Uh);
  double EnergyError = model.EnergyNorm(hho, u - Uh)/model.EnergyNorm(hho, Uh);
  //output << "H1 error = " << H1error << ", computed in " << double(timer.elapsed().wall) * pow(10, -9) << " seconds.\n";
  output << "  H1 Error = " << H1error << "\n";
  output << "  Energy Error = " << EnergyError << "\n";

  // --------------------------------------------------------------------------
  //                     Creates a .vtu file of the solution
  // --------------------------------------------------------------------------
  // Only if we do not have too many cells
  if (vm.count("plot") && mesh_ptr->n_cells() <= max_nb_cells_for_plot) {
    output << "\n[Scheme: writing solution to file]" << std::endl;
    std::string filename = vm["plot"].as<std::string>() + std::string(".vtu");
    VtuWriter plotdata(mesh_ptr.get());
		
    // Approximate solution, plots using cell or edge polynomials
//    		Eigen::VectorXd approx_sol_vertex_byT = hho.VertexValues(u, "cell");
//    		plotdata.write_to_vtu("T-" + filename,approx_sol_vertex_byT);

    Eigen::VectorXd approx_sol_vertex_byF = hho.VertexValues(u, "edge");
    plotdata.write_to_vtu(filename,approx_sol_vertex_byF);

    // Exact solution
    Eigen::VectorXd exact_sol_vertex = Eigen::VectorXd::Zero(mesh_ptr->n_vertices());
    for (size_t iV=0; iV< mesh_ptr->n_vertices(); iV++){
      auto v= mesh_ptr->vertex(iV)->coords();
      exact_sol_vertex(iV) = exact_solution(v);
    }
    plotdata.write_to_vtu(std::string("exact-")+filename,exact_sol_vertex);

  }

  // --------------------------------------------------------------------------
  //                     Creates .txt file with data and results
  // --------------------------------------------------------------------------

  std::ofstream out("results.txt");
  out << "BC: " << bc_id << '\n';
  out << "Solution: " << id_tcase[0] << '\n';
  out << "Diffusion: " << id_tcase[1] << '\n';
  out << "Mesh: " << mesh_name << "\n";
  out << "EdgeDegree: " << K << "\n";
  out << "CellDegree: " << L << "\n";
  out << "Using threads: " << (use_threads ? "true" : "false") << "\n";
  out << "AssemblyTime: " << model.get_assembly_time() << "\n";
  out << "Solving time: " << model.get_solving_time() << "\n";
  out << "L2error: " << L2error << "\n";
  out << "H1error: " << H1error << "\n";
  out << "EnergyError: " << EnergyError << "\n";
  out << "MeshSize: " << mesh_ptr->h_max() << "\n";
  out << "NbCells: " << mesh_ptr->n_cells() << "\n";
  out << "NbEdges: " << mesh_ptr->n_edges() << "\n";
  out << "NbEdgeDOFs: " << nbedgedofs << "\n";
	out << "MeshReg: " << meshreg[0] << "\n";
  out << "MeshSkewness: " << meshreg[1] << "\n";
  out << std::flush;
  out.close();

}
