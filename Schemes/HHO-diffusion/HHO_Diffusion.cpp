// Implementation of the HHO scheme in 2D for the diffusion equation, with K piecewise constant
//
//   { -div(K \grad(u)) = f,       inside Omega
//   { K \grad(u) . nTF = g,       on GammaN
//   {                 u = g,       on GammaD
// 
//  At the moment: pure Dirichlet or pure Neumann
//
// Author: Jerome Droniou (jerome.droniou@monash.edu)
//

#include <iostream>
#include <fstream>

#include <boost/program_options.hpp>
#include <boost/timer/timer.hpp>

#include "mesh.hpp"
#include "import_mesh.hpp"
#include "mesh_builder.hpp"

#include "HHO_Diffusion.hpp"
#include "hybridcore.hpp"
#include "TestCase/TestCase.hpp"
#include "vtu_writer.hpp"

using namespace HArDCore2D;

// Mesh filenames
const std::string mesh_dir = "../../typ2_meshes/";
std::string default_mesh = mesh_dir + "cart5x5.typ2";
const std::string default_solver_type = "bicgstab";

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
      ("choice_basis", boost::program_options::value<std::string>()->default_value("Mon"), "Choice of basis function: 'Mon' for monomials, 'ON' for orthonormalised")
      ("BC,b", boost::program_options::value<size_t>()->default_value(0), "Set the boundary conditions (0=Dirichlet, 1=Neumann)")
      ("testcase,c", boost::program_options::value<std::vector<int>>()->multitoken(), "Set the test case (as '-c n m'; n=exact sol, m=diffusion)")
      ("plot,p", boost::program_options::value<std::string>()->default_value("solution"), "Save plot of the solution to the given filename")
      ("solver_type", boost::program_options::value<std::string>()->default_value("bicgstab"), "Defines the linear solver for the system");


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

  // Select the degree
  size_t K = vm["edgedegree"].as<size_t>();
  int L = vm["celldegree"].as<size_t>();

  // Choice basis
  std::string choice_basis = vm["choice_basis"].as<std::string>();

  // Select solver type
  std::string solver_type = (vm.count("solver_type") ? vm["solver_type"].as<std::string>() : default_solver_type);

  // Checks
  if ( (abs(int(K)-L) > 1) || (K<0) || (L<-1) ){
    std::cout << "Degrees k and l not in acceptable range (k positive and l=k-1, k, or k+1): k=" << K << ", l=" << L << "\n\n";
std::cout << std::abs(int(K)-L);
    exit(1);
  }

  // --------------------------------------------------------------------------
  //                        Create the HHO data structure
  // --------------------------------------------------------------------------

  // Read the mesh file
  MeshReaderTyp2 mesh(mesh_file);

  std::vector<std::vector<double> > vertices;
  std::vector<std::vector<size_t> > cells;
  std::vector<std::vector<double> > centers;
  if (mesh.read_mesh(vertices, cells, centers) == false) {
    std::cout << "Could not open file" << std::endl;
    return false;
  };

  // Build the mesh
  MeshBuilder builder = MeshBuilder();
  std::unique_ptr<Mesh> mesh_ptr = builder.build_the_mesh(vertices, cells);
  if (mesh_ptr.get() == NULL) {
    printf(
      "Mesh cannot be created!\n Check the input file contains \n "
      "Vertices "
      "and cells with the correct tags");
    return 0;
  } 
  // Re-index the mesh edges, to facilitate treatment of boundary conditions (the boundary edges are put at the end)
  std::vector<size_t> new_to_old(mesh_ptr->n_edges(), 0);
  for (size_t idx = 0; idx < mesh_ptr->n_i_edges(); idx++){
    new_to_old[idx] = mesh_ptr->i_edge(idx)->global_index();
  }
  for (size_t idx = 0; idx < mesh_ptr->n_b_edges(); idx++){
    new_to_old[mesh_ptr->n_i_edges()+idx] = mesh_ptr->b_edge(idx)->global_index();
  }
  mesh_ptr->renum('E', new_to_old);

  // Create the HHO structure
  HybridCore hho(mesh_ptr.get(), K, L, choice_basis);

  // --------------------------------------------------------------------------
  //                        Create the model equation
  // --------------------------------------------------------------------------

  // Select the boundary conditions
  size_t BC = vm["BC"].as<size_t>();
  
  const std::vector<int> default_id_tcase = std::vector<int>{1,1};
  std::vector<int> id_tcase = (vm.count("testcase") ? vm["testcase"].as<std::vector<int>>() : default_id_tcase);
  TestCase tcase(id_tcase);

  // Diffusion tensor
  HHO_Diffusion::tensor_function_type kappa = [&](double x, double y, Cell* cell) {
      return tcase.diff(x,y,cell);
  };

  // Source term
  HHO_Diffusion::source_function_type source = [&](double x, double y, Cell* cell) {
      return tcase.source(x,y, cell);
  };

  // Exact solution and gradient
  HHO_Diffusion::solution_function_type exact_solution = [&](double x, double y) {
      return tcase.sol(x,y);
  };
  HHO_Diffusion::grad_function_type grad_exact_solution = [&](double x, double y, Cell* cell) {
      return tcase.grad_sol(x,y,cell);
  };

  // Create the model equation
  HHO_Diffusion model(kappa, source, BC, exact_solution, grad_exact_solution, solver_type);


  // --------------------------------------------------------------------------
  //                       Recalling the mesh and parameters
  // --------------------------------------------------------------------------
  if (BC==0){
    std::cout << "\nBoundary conditions: Dirichlet\n";
  } else if (BC==1){
    std::cout << "\nBoundary conditions: Neumann\n";
  }
  std::cout << "Test case: solution = " << id_tcase[0] << "; diffusion = " << id_tcase[1] << "\n";
  size_t found = mesh_file.find_last_of("/\\");
  std::string mesh_name = mesh_file.substr(found+1);
  std::vector<double> meshreg = mesh_ptr->regularity();
  std::cout << "Mesh = " << mesh_name << " (nb cells= " << mesh_ptr->n_cells() << ", nb edges= " << mesh_ptr->n_edges() << ", reg= " << meshreg[0] << ", " << meshreg[1] << ")\n";
  size_t nbedgedofs = mesh_ptr->n_edges()*hho.nlocal_edge_dofs();
  std::cout << "Degrees: edge = " << K << "; cell = " << L << " | Nb edge DOFs = " << nbedgedofs << "\n";
  std::cout << "Type basis functions: " << choice_basis << "\n\n";

  // --------------------------------------------------------------------------
  //                        Solve the model problem
  // --------------------------------------------------------------------------

  std::cout << std::string(80, '-') << "\nAssembling and solving the problem..." << std::endl;
  boost::timer::cpu_timer timer;
  Eigen::VectorXd u = model.solve(hho);
  std::cout<<"Model is solved"<<std::endl;
  
  std::cout << "Solver= " << solver_type << ", solved problem in " << double(timer.elapsed().wall) * pow(10, -9) << " seconds\n";
  std::cout << "\tAssembly time = " << model.get_assembly_time() << "s\n";
  std::cout << "\tSolving time = " << model.get_solving_time() << "s\n";
  std::cout << "\tResidual of the linear system = " << model.get_solving_error() << '\n';

  printf("interm times = %f | %f | %f | %f | %f\n", 
      model.get_itime(0),  model.get_itime(1),  model.get_itime(2),  model.get_itime(3), model.get_itime(4));
//  printf("interm times = %f | %f | %f | %f | %f \n", 
//      model.get_itime(5), model.get_itime(6), model.get_itime(7), model.get_itime(8), model.get_itime(9));

  // --------------------------------------------------------------------------
  //                        Compute the error
  // --------------------------------------------------------------------------

  std::cout << "Computing error..." << std::endl;

  // Interpolate the exact solution
  //timer.start();
  Eigen::VectorXd Uh = hho.interpolate(exact_solution, 2*hho.K()+5);
  std::cout << "Interpolant calculated"<<std::endl;
  //std::cout << "Interpolant of exact solution computed in " << double(timer.elapsed().wall) * pow(10, -9) << " seconds.\n";

  // Compute the L2 error
  //timer.start();
  double L2error = hho.L2norm(u - Uh)/hho.L2norm(Uh);
  //std::cout << "L2 Error = " << L2error << ", computed in " << double(timer.elapsed().wall) * pow(10, -9) << " seconds.\n";
  std::cout << "L2 Error = " << L2error << "\n";

  // Compute the error in discrete H1 norm and energy norm
  //timer.start();
  double H1error = hho.H1norm(u - Uh)/hho.H1norm(Uh);
  double EnergyError = model.EnergyNorm(hho,u - Uh)/model.EnergyNorm(hho,Uh);
  //std::cout << "H1 error = " << H1error << ", computed in " << double(timer.elapsed().wall) * pow(10, -9) << " seconds.\n";
  std::cout << "H1 Error = " << H1error << "\n";
  std::cout << "Energy Error = " << EnergyError << "\n";

  // Compute the Linf error of the edge coefficients
  //timer.start();
//  double erroredges = hho.Linf_edge(u - Uh)/hho.Linf_edge(Uh);
//  std::cout << "Linf edges = " << erroredges << "\n";

  // Compute the integral
//  double average = hho.integral(u);
//  double average_exact_sol = hho.integrate_over_domain([&](auto x, auto y) {
//    return exact_solution(x,y);
//    });
//  std::cout << "Integral of solution = " << average << " (should be " << average_exact_sol << ")\n"; 

  // --------------------------------------------------------------------------
  //                     Creates a .vtu file of the solution
  // --------------------------------------------------------------------------

  if (vm.count("plot")) {
    std::string filename = vm["plot"].as<std::string>() + std::string(".vtu");
    VtuWriter plotdata(mesh_ptr.get());
    
    // Approximate solution, plots using cell or edge polynomials
//    Eigen::VectorXd approx_sol_vertex_byT = hho.VertexValues(u, "cell");
//    plotdata.write_to_vtu("T-" + filename,approx_sol_vertex_byT,1);

    Eigen::VectorXd approx_sol_vertex_byE = hho.VertexValues(u, "edge");
    plotdata.write_to_vtu(filename,approx_sol_vertex_byE,1);

    // Exact solution
    Eigen::VectorXd exact_sol_vertex = Eigen::VectorXd::Zero(mesh_ptr->n_vertices());
    for (size_t iV=0; iV< mesh_ptr->n_vertices(); iV++){
      auto v= mesh_ptr->vertex(iV)->coords();
      exact_sol_vertex(iV) = exact_solution(v(0),v(1));
    }
    plotdata.write_to_vtu(std::string("exact-")+filename,exact_sol_vertex,1);

    std::cout << "\nApproximate solution in '" << filename << "', exact solution in '" << std::string("exact-")+filename << "'\n";

}

    // --------------------------------------------------------------------------
    //                     Creates .txt file with data and results
    // --------------------------------------------------------------------------

    std::ofstream out("results.txt");
    out << "BC: " << BC << '\n';
    out << "Solution: " << id_tcase[0] << '\n';
    out << "Diffusion: " << id_tcase[1] << '\n';
    out << "Mesh: " << mesh_name << "\n";
    out << "EdgeDegree: " << K << "\n";
    out << "CellDegree: " << L << "\n";
    out << "Choice basis: " << choice_basis << "\n";
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
    out << "lambda from TestCase: " << tcase.get_lambda() << "\n";
    out << std::flush;
    out.close();


}
