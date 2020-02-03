// Implementation of the Non-conforming scheme on generic polygonal meshes for the Stefan and PME models
//
//   { u - div(K \grad(zeta(u))) = f,       inside Omega
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

#include "BPNC_StefanPME.hpp"
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
      ("testcaseNL", boost::program_options::value<int>()->default_value(1), "Set the nonlinearity zeta (as '--testcaseNL n')")
      ("powerPME", boost::program_options::value<double>()->default_value(1.0), "Set the power of the nonlinearity, only useful when the PME nonlinearity is chosen (as '--powerPME m')")
      ("weight", boost::program_options::value<double>()->default_value(0), "Value between 0 and 1, set the proportion of weight of mass-lumping on the edge unknowns (as '--weight w')")
      ("source", boost::program_options::value<int>()->default_value(1), "1= use correct source for given testcase, 0= put zero source (as '--source s')")
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
	// Re-index the mesh edges, to facilitate treatment of boundary conditions (Dirichlet boundary edges are put at the end)
  // Get boundary conditions and re-order edges
  std::string bc_id = vm["bc_id"].as<std::string>();
  BoundaryConditions BC(bc_id, *mesh_ptr.get());
  BC.reorder_edges(); 

  // Create the NC structure
  BPNCCore nc(mesh_ptr.get(), 1, 1, output);

  // --------------------------------------------------------------------------
  //                        Create the model equation
  // --------------------------------------------------------------------------

	if (BC.name()=="Neumann"){
		output << "\n\n@@@@@@@@@@@@@@@@@@@@@@@@@\n          WARNING: Neumann only works if there is no reaction term!! \n\n@@@@@@@@@@@@@@@@@@@@@@@@@\n\n";
	}
	
	const std::vector<int> default_id_tcase = std::vector<int>{1,1};
	std::vector<int> id_tcase = (vm.count("testcase") ? vm["testcase"].as<std::vector<int>>() : default_id_tcase);
	TestCaseStefanPME tcase(id_tcase);

	int id_tcaseNL = vm["testcaseNL"].as<int>();
	int mPME = vm["powerPME"].as<double>();
	TestCaseNonLinearity tcaseNL(id_tcaseNL, mPME);

	double weight = vm["weight"].as<double>();
	double source = vm["source"].as<int>();

	// Nonlinearity zeta(u)
	TestCaseNonLinearity::nonlinearity_function_type zeta = [&tcaseNL](double s, std::string type){
		return tcaseNL.nonlinearity(s, type);
	};

  // Diffusion tensor
  BPNC_StefanPME::tensor_function_type kappa = [&](double x, double y, Cell* cell) {
			return tcase.diff(x,y,cell);
  };
	size_t deg_kappa = tcase.get_deg_diff();

  // Exact solution and gradient
  BPNC_StefanPME::solution_function_type exact_solution = [&](double x, double y) {
			return tcase.sol(x,y);
  };
  BPNC_StefanPME::grad_function_type grad_exact_solution = [&](double x, double y, Cell* cell) {
			return tcase.grad_sol(x,y,cell);
  };
	
  // Source term, adjusted for reaction term
  BPNC_StefanPME::source_function_type source_term = [&](double x, double y, Cell* cell) {
		// Two options: zero source term, or exact given by exact solution (latter is only valid for diffusion=Id)
		double val = 0;
		if (source == 1){
			val = exact_solution(x, y) - zeta(exact_solution(x, y), "hess") * (kappa(x, y, cell)*grad_exact_solution(x, y, cell)).dot(grad_exact_solution(x, y, cell))	+ zeta(exact_solution(x, y), "der")*tcase.source(x,y, cell);
		}
		return val;
  };

  // Create the model equation
  BPNC_StefanPME StefanPME_model(nc, kappa, deg_kappa, source_term, BC, exact_solution, grad_exact_solution, zeta, weight, solver_type);

  // --------------------------------------------------------------------------
	// 											Recalling the mesh and parameters
  // --------------------------------------------------------------------------
  output << " Boundary conditions: " << BC.name() << " /";
	output << "Test case: solution = " << id_tcase[0] << "; diffusion = " << id_tcase[1] << "; nonlinearity = " << id_tcaseNL << " (power PME: " << mPME << "); weight for edges: " << weight << "\n";
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
  Eigen::VectorXd Xh = StefanPME_model.solve();
  output<<"Model is solved"<<std::endl;
  

  // --------------------------------------------------------------------------
  //                        Compute the error
  // --------------------------------------------------------------------------

  output << "Computing error..." << std::endl;

	// Compute vector for zeta of approximate solution
	Eigen::VectorXd zetau = StefanPME_model.apply_nonlinearity(Xh, "fct");

  // Interpolate the exact solution
  //timer.start();
  Eigen::VectorXd Uh = nc.nc_interpolate_ml(exact_solution, 10);
	std::function<double(double,double)> exact_zetau = [&zeta, &exact_solution](double x, double y)->double {
			return zeta(exact_solution(x, y), "fct");
		};
  Eigen::VectorXd ZetaUh = nc.nc_interpolate_ml(exact_zetau, 10);
  output << "Interpolant calculated"<<std::endl;

  // Compute the mass-lumped L2 error
  double L2error_ml = StefanPME_model.L2_MassLumped(Xh - Uh)/StefanPME_model.L2_MassLumped(Uh);
	output << "L2 Error ML = " << L2error_ml << "\n";

  // Compute the error in discrete H1 norm and energy norm
  //timer.start();
  double H1error = nc.nc_H1norm(zetau - ZetaUh)/nc.nc_H1norm(ZetaUh);
  double EnergyError = StefanPME_model.EnergyNorm(zetau - ZetaUh)/StefanPME_model.EnergyNorm(ZetaUh);
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
		// In the cases weight=0 and weight=1, some components of Xh actually correspond to zeta(u), not u,
		// but nc_VertexValues then does not use these components
		Eigen::VectorXd approx_sol_vertex = nc.nc_VertexValues(Xh, weight);
		plotdata.write_to_vtu(filename,approx_sol_vertex,1);

		// Zeta of approximate solution
		Eigen::VectorXd zeta_approx_sol_vertex = nc.nc_VertexValues(zetau, weight);
		std::string filename_zeta = std::string("zeta-") + vm["plot"].as<std::string>() + std::string(".vtu");
		plotdata.write_to_vtu(filename_zeta,zeta_approx_sol_vertex,1);

		// Exact solution
		Eigen::VectorXd exact_sol_vertex = Eigen::VectorXd::Zero(mesh_ptr->n_vertices());
		for (size_t iV=0; iV< mesh_ptr->n_vertices(); iV++){
			auto v= mesh_ptr->vertex(iV)->coords();
			exact_sol_vertex(iV) = exact_solution(v(0),v(1));
		}
		plotdata.write_to_vtu(std::string("exact-")+filename,exact_sol_vertex,1);

		// Zeta of exact solution
		Eigen::VectorXd zeta_exact_sol_vertex = Eigen::VectorXd::Zero(mesh_ptr->n_vertices());
		for (size_t iV=0; iV< mesh_ptr->n_vertices(); iV++){
			auto v= mesh_ptr->vertex(iV)->coords();
			zeta_exact_sol_vertex(iV) = zeta(exact_solution(v(0),v(1)), "fct");
		}
		plotdata.write_to_vtu(std::string("zeta-exact-")+filename,zeta_exact_sol_vertex,1);

	}

//		// --------------------------------------------------------------------------
//	  //                     Creates .txt file with data and results
//	  // --------------------------------------------------------------------------

    std::ofstream out("results.txt");
    out << "BC: " << bc_id << '\n';
    out << "Solution: " << id_tcase[0] << '\n';
    out << "Diffusion: " << id_tcase[1] << '\n';
    out << "Nonlinearity: " << id_tcaseNL << '\n';
    out << "Power PME: " << mPME << '\n';
    out << "Weight: " << weight << '\n';
		out << "Mesh: " << mesh_name << "\n";
		out << "AssemblyTime: " << StefanPME_model.get_assembly_time() << "\n";
		out << "Solving time: " << StefanPME_model.get_solving_time() << "\n";
		out << "L2error: " << L2error_ml << "\n";
		out << "L2error_ml: " << L2error_ml << "\n";
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
