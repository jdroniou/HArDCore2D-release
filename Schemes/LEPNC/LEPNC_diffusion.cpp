// Implementation of the Locally Enriched Non-Conforming scheme on generic polygonal meshes for the diffusion problem
//
//   { - div(K \grad(zeta(u))) = f,       inside Omega
//   { K \grad(zeta(u)) . nTF = g,       on GammaN
//   { 								u = g,			 on GammaD
// 
// At the moment, only pure Neumann or pure Dirichlet
//
// Author: Jerome Droniou (jerome.droniou@monash.edu)
//

#include <iostream>

#include <boost/program_options.hpp>
#include <boost/timer/timer.hpp>

#include "mesh.hpp"
#include "import_mesh.hpp"
#include "mesh_builder.hpp"

#include "LEPNC_diffusion.hpp"
#include "TestCase/TestCase.hpp"
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
      ("bc_id,b", boost::program_options::value<std::string>()->default_value("D"), "Set the boundary conditions (D=Dirichlet, N=Neumann, Mx=Mixed number x)")
      ("testcase,c", boost::program_options::value<std::vector<int>>()->multitoken(), "Set the test case (as '-c n m'; n=exact sol, m=diffusion)")
      ("plot,p", boost::program_options::value<std::string>()->default_value("solution"), "Save plot of the solution to the given filename")
      ("solver_type", boost::program_options::value<std::string>()->default_value("bicgstab"), "Defines the linear solver for the system");

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

	// Select solver type
	std::string solver_type = (vm.count("solver_type") ? vm["solver_type"].as<std::string>() : default_solver_type);


	// --------------------------------------------------------------------------
  //                        Create the NC data structure
  // --------------------------------------------------------------------------

  // Read the mesh file
	MeshReaderTyp2 mesh(mesh_file);

	std::vector<std::vector<double> > vertices;
	std::vector<std::vector<size_t> > cells;
	std::vector<std::vector<double> > centers;
	if (mesh.read_mesh(vertices, cells, centers) == false) {
		output << "Could not open file" << std::endl;
		return 0;
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
	// Re-index the mesh edges, to facilitate treatment of boundary conditions (Dirichlet boundary edges are put at the end)
  // Get boundary conditions and re-order edges
  std::string bc_id = vm["bc_id"].as<std::string>();
  BoundaryConditions BC(bc_id, *mesh_ptr.get());
  BC.reorder_edges(); 

  // Create the NC structure
  LEPNCCore nc(mesh_ptr.get(), 1, 1, output);

  // --------------------------------------------------------------------------
  //                        Create the model equation
  // --------------------------------------------------------------------------

	if (BC.name()=="Neumann"){
		output << "\n\n@@@@@@@@@@@@@@@@@@@@@@@@@\n          WARNING: Neumann only works if there is no reaction term!! \n\n@@@@@@@@@@@@@@@@@@@@@@@@@\n\n";
	}
	
	const std::vector<int> default_id_tcase = std::vector<int>{1,1};
	std::vector<int> id_tcase = (vm.count("testcase") ? vm["testcase"].as<std::vector<int>>() : default_id_tcase);
	TestCase tcase(id_tcase);

  // Diffusion tensor
   LEPNC_diffusion::tensor_function_type kappa = [&](double x, double y, Cell* cell) {
      VectorRd p = VectorRd(x,y);
			return tcase.get_diffusion().value(p,cell);
  };
	size_t deg_kappa = tcase.get_diffusion().degree;

  // Exact solution and gradient
  LEPNC_diffusion::solution_function_type exact_solution = [&](double x, double y) {
      VectorRd p = VectorRd(x,y);
			return tcase.get_solution().value(p);
  };
  LEPNC_diffusion::grad_function_type grad_exact_solution = [&](double x, double y, Cell* cell) {
      VectorRd p = VectorRd(x,y);
			return tcase.get_solution().gradient(p,cell);
  };

	
  // Source term
  LEPNC_diffusion::source_function_type source_term = [&](double x, double y, Cell* cell) {
    VectorRd p = VectorRd(x,y);
		return tcase.diff_source()(p, cell);
  };

  // Create the model equation
  LEPNC_diffusion model(nc, kappa, deg_kappa, source_term, BC, exact_solution, grad_exact_solution, solver_type);

  // --------------------------------------------------------------------------
	// 											Recalling the mesh and parameters
  // --------------------------------------------------------------------------
  output << " Boundary conditions: " << BC.name() << " /";
	output << "Test case: solution = " << id_tcase[0] << "; diffusion = " << id_tcase[1] << "\n";
	size_t found = mesh_file.find_last_of("/\\");
	std::string mesh_name = mesh_file.substr(found+1);
  std::vector<double> meshreg = mesh_ptr->regularity();
  output << "Mesh = " << mesh_name << " (nb cells= " << mesh_ptr->n_cells() << ", nb edges= " << mesh_ptr->n_edges() << ", reg= " << meshreg[0] << ", " << meshreg[1] << ")\n";
	size_t nbedgedofs = mesh_ptr->n_edges();
	output << "Nb edge DOFs = " << nbedgedofs << "\n\n";

  // --------------------------------------------------------------------------
  //                        Solve the model problem
  // --------------------------------------------------------------------------

  output << std::string(80, '-') << "\nAssembling and solving the problem..." << std::endl;
  boost::timer::cpu_timer timer;
  Eigen::VectorXd Xh = model.solve();
  output<<"Model is solved"<<std::endl;
  

  // --------------------------------------------------------------------------
  //                        Compute the error
  // --------------------------------------------------------------------------

  output << "Computing error..." << std::endl;

  // Interpolate the exact solution
  //timer.start();
  Eigen::VectorXd Uh_ml = nc.nc_interpolate_ml(exact_solution, 10);
  Eigen::VectorXd Uh_moments = nc.nc_interpolate_moments(exact_solution, 10);
  output << "Interpolant calculated"<<std::endl;

  // Compute the mass-lumped L2 error
  double L2error = nc.nc_L2norm(Xh - Uh_moments)/nc.nc_L2norm(Uh_moments);
	output << "L2 Error = " << L2error << "\n";

  // Compute the error in discrete H1 norm and energy norm
  //timer.start();
  double H1error = nc.nc_H1norm(Xh - Uh_moments)/nc.nc_H1norm(Uh_moments);
  double EnergyError = model.EnergyNorm(Xh - Uh_moments)/model.EnergyNorm(Uh_moments);
	output << "H1 Error = " << H1error << "\n";
	output << "Energy Error = " << EnergyError << "\n";

	// --------------------------------------------------------------------------
  //                     Creates .vtu files
  // --------------------------------------------------------------------------
  // Only if we do not have too many cells
  if (vm.count("plot") && mesh_ptr->n_cells() <= max_nb_cells_for_plot) {
		std::string filename = vm["plot"].as<std::string>() + std::string(".vtu");
		VtuWriter plotdata(mesh_ptr.get());
		
		// Approximate solution
		Eigen::VectorXd approx_sol_vertex = nc.nc_VertexValues(Xh, 0.5);
		plotdata.write_to_vtu(filename,approx_sol_vertex,1);

		// Exact solution
		Eigen::VectorXd exact_sol_vertex = Eigen::VectorXd::Zero(mesh_ptr->n_vertices());
		for (size_t iV=0; iV< mesh_ptr->n_vertices(); iV++){
			auto v= mesh_ptr->vertex(iV)->coords();
			exact_sol_vertex(iV) = exact_solution(v(0),v(1));
		}
		plotdata.write_to_vtu(std::string("exact-")+filename,exact_sol_vertex,1);

	}

//		// --------------------------------------------------------------------------
//	  //                     Creates .txt file with data and results
//	  // --------------------------------------------------------------------------

    std::ofstream out("results.txt");
    out << "BC: " << bc_id << '\n';
    out << "Solution: " << id_tcase[0] << '\n';
    out << "Diffusion: " << id_tcase[1] << '\n';
		out << "Mesh: " << mesh_name << "\n";
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
		out << "MeshAnisotropy: " << meshreg[1] << "\n";
		out << std::flush;
    out.close();


}
