// Implementation of the Non-conforming scheme on generic polygonal meshes for the Stefan and PME models
//
//   { d_t u - div(K \grad(zeta(u))) = f,       inside Omega
//   { K \grad(zeta(u)) . nTF = g,       on GammaN
//   { 								u = g,			 on GammaD
//   { u(0) = u_ini
//
// At the moment, only pure Neumann or pure Dirichlet
//
// Author: Jerome Droniou (jerome.droniou@monash.edu)
//

#include <iostream>

#include <boost/program_options.hpp>
#include <boost/timer/timer.hpp>

#include "mesh.hpp"
//#include "import_mesh.hpp"
//#include "mesh_builder.hpp"

#include "LEPNC_StefanPME_transient.hpp"
#include "TestCase/TestCaseTransient.hpp"
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
      ("FinalTps", boost::program_options::value<double>()->default_value(1.0), "Set the final time")
      ("dt", boost::program_options::value<double>()->default_value(0.1), "Set the time step")
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

  // Select the mesh, final time and time step
	std::string mesh_file = (vm.count("mesh") ? vm["mesh"].as<std::string>() : default_mesh);
	double FinalTps = vm["FinalTps"].as<double>();
	double dt = vm["dt"].as<double>();

	// Select solver type
	std::string solver_type = (vm.count("solver_type") ? vm["solver_type"].as<std::string>() : default_solver_type);


	// --------------------------------------------------------------------------
  //                        Create the NC data structure
  // --------------------------------------------------------------------------

  // Read the mesh file
//	MeshReaderTyp2 mesh(mesh_file);

//	std::vector<std::vector<double> > vertices;
//	std::vector<std::vector<size_t> > cells;
//	std::vector<std::vector<double> > centers;
//	if (mesh.read_mesh(vertices, cells, centers) == false) {
//		output << "Could not open file" << std::endl;
//		return 0;
//	};

//	// Build the mesh
//	MeshBuilder builder = MeshBuilder();
//	std::unique_ptr<Mesh> mesh_ptr = builder.build_the_mesh(vertices, cells);
//	if (mesh_ptr.get() == NULL) {
//		printf(
//		  "Mesh cannot be created!\n Check the input file contains \n "
//		  "Vertices "
//		  "and cells with the correct tags");
//		return 0;
//	} 

	MeshBuilder builder(mesh_file);
	std::unique_ptr<Mesh> mesh_ptr = builder.build_the_mesh();
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
	
	int id_tcaseNL = vm["testcaseNL"].as<int>();
	double mPME = vm["powerPME"].as<double>();
	TestCaseNonLinearity tcaseNL(id_tcaseNL, mPME);

	const std::vector<int> default_id_tcase = std::vector<int>{1,1};
	std::vector<int> id_tcase = (vm.count("testcase") ? vm["testcase"].as<std::vector<int>>() : default_id_tcase);
	TestCaseTransient tcase(id_tcase, mPME);

	double weight = vm["weight"].as<double>();
	int source = vm["source"].as<int>();

	// Nonlinearity zeta(u)
	TestCaseNonLinearity::nonlinearity_function_type zeta = [&tcaseNL](double s, std::string type){
		return tcaseNL.nonlinearity(s, type);
	};

  // Diffusion tensor
  LEPNC_StefanPME_Transient::tensor_function_type kappa = [&](const double x, const double y, const Cell* cell) {
			return tcase.diff(x,y,cell);
  };
	size_t deg_kappa = tcase.get_deg_diff();

  // Exact solution, gradient and time derivative
  LEPNC_StefanPME_Transient::solution_function_type exact_solution = [&](const double t, const VectorRd p) {
			return tcase.solution()(t,p);
  };
  LEPNC_StefanPME_Transient::grad_function_type grad_exact_solution = [&](const double t, const VectorRd p, const Cell* cell) {
			return tcase.grad_solution()(t,p,cell);
  };
  LEPNC_StefanPME_Transient::solution_function_type time_der_exact_solution = [&](double t, VectorRd p) {
			return tcase.delt_solution()(t,p);
  };
	
  // Source term, adjusted for reaction term
  LEPNC_StefanPME_Transient::source_function_type source_term = [&](const double t, const VectorRd p, const Cell* cell) {
		// Two options: zero source term, or exact given by exact solution (latter is only valid for diffusion=Id)
		double val = 0;
		if (source == 1){
			val = time_der_exact_solution(t, p) - zeta(exact_solution(t, p), "hess") * (kappa(p.x(), p.y(), cell)*grad_exact_solution(t, p, cell)).dot(grad_exact_solution(t, p, cell)) + zeta(exact_solution(t, p), "der")*tcase.minus_div_diff_grad(t, p.x(), p.y(), cell);
		}
		return val;
  };

  // Create the model equation
  LEPNC_StefanPME_Transient StefanPME_model(nc, kappa, deg_kappa, source_term, BC, exact_solution, grad_exact_solution, zeta, weight, solver_type);

  // --------------------------------------------------------------------------
	// 											Recalling the mesh and parameters
  // --------------------------------------------------------------------------
  output << " Boundary conditions: " << BC.name() << " /";
	output << "Test case: solution = " << id_tcase[0] << "; diffusion = " << id_tcase[1] << "; nonlinearity = " << id_tcaseNL << " (power PME: " << mPME << "); weight for edges: " << weight << "\n";
	size_t found = mesh_file.find_last_of("/\\");
	std::string mesh_name = mesh_file.substr(found+1);
  std::vector<double> meshreg = mesh_ptr->regularity();
  output << "Mesh = " << mesh_name << " (nb cells= " << mesh_ptr->n_cells() << ", nb edges= " << mesh_ptr->n_edges() << ", reg= " << meshreg[0] << ", " << meshreg[1] << ")\n";
  output << "dt = " << dt << "\n";
	size_t nbedgedofs = mesh_ptr->n_i_edges();
	output << "Nb edge DOFs = " << nbedgedofs << "\n\n";

  // --------------------------------------------------------------------------
  //                        Time iterations
  // --------------------------------------------------------------------------

  output << std::string(80, '-') << "\nTime stepping..." << std::endl;
  boost::timer::cpu_timer timer;
  size_t NbTimeSteps = FinalTps/dt;
  
  // Initial solution
  std::function<double(double, double)> initial_condition = [&](double x, double y)->double{
      return exact_solution(0.0, VectorRd(x,y));
    };
  Eigen::VectorXd Xn = nc.nc_interpolate_ml(initial_condition, 10);
  double ave_newton = 0;
  for (size_t idt = 1; idt <= NbTimeSteps; idt++){
    output << "Time step " << idt << "/" << NbTimeSteps << std::endl;
    Eigen::VectorXd Xnp1 = StefanPME_model.iterate(dt*idt, dt, Xn);
    Xn = Xnp1;
    ave_newton += StefanPME_model.get_nb_newton();
  } 
  ave_newton /= NbTimeSteps;
  
  output<<"Time stepping finished (averaged number of Newton iterations: "<< ave_newton << ")" << std::endl;
  

  // --------------------------------------------------------------------------
  //                        Compute the error
  // --------------------------------------------------------------------------

  output << "Computing error..." << std::endl;

	// Compute vector for zeta of approximate solution
	Eigen::VectorXd zetau = StefanPME_model.apply_nonlinearity(Xn, "fct");

  // Interpolate the exact solution
  //timer.start();
  std::function<double(double,double)> solution_FinalTps = [&](double x, double y)->double{
    return exact_solution(FinalTps, VectorRd(x,y));
  };
  Eigen::VectorXd Uh = nc.nc_interpolate_ml(solution_FinalTps, 10);
	std::function<double(double,double)> exact_zetau = [&zeta, &solution_FinalTps](double x, double y)->double {
			return zeta(solution_FinalTps(x, y), "fct");
		};
  Eigen::VectorXd ZetaUh = nc.nc_interpolate_ml(exact_zetau, 10);
  output << "Interpolant calculated" << std::endl;

  // Compute the mass-lumped L2 and L^{m+1} errors
  double L2error = StefanPME_model.L2_MassLumped(Xn - Uh)/StefanPME_model.L2_MassLumped(Uh);
	output << "L2 Error ML = " << L2error << "\n";
  double Lmp1error = StefanPME_model.Lp_MassLumped(Xn - Uh, mPME+1)/StefanPME_model.Lp_MassLumped(Uh, mPME+1);
	output << "Lm+1 Error ML (m+1=" << mPME+1.0 << ")= " << Lmp1error << "\n";

  // Compute the error in discrete H1 norm and energy norm
  //timer.start();
  double H1error = nc.nc_H1norm(zetau - ZetaUh)/nc.nc_H1norm(ZetaUh);
  double EnergyError = StefanPME_model.EnergyNorm(zetau - ZetaUh)/StefanPME_model.EnergyNorm(ZetaUh);
	output << "H1 Error = " << H1error << "\n";
	output << "Energy Error = " << EnergyError << "\n";

  // Percentage (area/total_area) of negative values
  double neg_area = 0;
  double total_area = 0;
  double min_value = 0;
  double max_value = 0;
  double neg_mass = 0;
  double total_mass = 0;
  for (size_t iT=0; iT < mesh_ptr->n_cells(); iT++){
    Cell* T = mesh_ptr->cell(iT);
    total_area += T->measure();
    Eigen::VectorXd XT = nc.nc_restr(Xn, iT);
    Eigen::MatrixXd MassLoc = StefanPME_model.get_MassT(iT);
    total_mass += (MassLoc * (XT.array().abs()).matrix()).sum();
    for (size_t i=0; i < size_t(XT.size()); i++){
      max_value = std::max(XT(i), max_value);
      if (XT(i) < 0) {
        neg_area += MassLoc(i,i);
        min_value = std::min(XT(i), min_value);
        neg_mass += MassLoc(i,i) * std::abs(XT(i));
      }
    }
  }
  double fraction_neg_values = neg_area/total_area;
  double fraction_neg_mass = neg_mass/total_mass;
  output << "Fraction negative values: " << fraction_neg_values << "; min/max values: " << min_value << "/" << max_value << "\n";
  output << "Fraction negative mass: " << fraction_neg_mass << "\n";
  

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
		Eigen::VectorXd approx_sol_vertex = nc.nc_VertexValues(Xn, weight);
		plotdata.write_to_vtu(filename,approx_sol_vertex);

		// Zeta of approximate solution
		Eigen::VectorXd zeta_approx_sol_vertex = nc.nc_VertexValues(zetau, weight);
		std::string filename_zeta = std::string("zeta-") + vm["plot"].as<std::string>() + std::string(".vtu");
		plotdata.write_to_vtu(filename_zeta,zeta_approx_sol_vertex);

		// Exact solution
		Eigen::VectorXd exact_sol_vertex = Eigen::VectorXd::Zero(mesh_ptr->n_vertices());
		for (size_t iV=0; iV< mesh_ptr->n_vertices(); iV++){
			auto v= mesh_ptr->vertex(iV)->coords();
			exact_sol_vertex(iV) = exact_solution(FinalTps,v);
		}
		plotdata.write_to_vtu(std::string("exact-")+filename,exact_sol_vertex);

		// Zeta of exact solution
		Eigen::VectorXd zeta_exact_sol_vertex = Eigen::VectorXd::Zero(mesh_ptr->n_vertices());
		for (size_t iV=0; iV< mesh_ptr->n_vertices(); iV++){
			auto v= mesh_ptr->vertex(iV)->coords();
			zeta_exact_sol_vertex(iV) = zeta(exact_solution(FinalTps,v), "fct");
		}
		plotdata.write_to_vtu(std::string("zeta-exact-")+filename,zeta_exact_sol_vertex);

	}

		// --------------------------------------------------------------------------
	  //                     Creates .txt file with data and results
	  // --------------------------------------------------------------------------

    std::ofstream out("results.txt");
    out << "BC: " << bc_id << '\n';
    out << "Solution: " << id_tcase[0] << '\n';
    out << "Diffusion: " << id_tcase[1] << '\n';
    out << "Nonlinearity: " << id_tcaseNL << '\n';
    out << "mPME: " << mPME << '\n';
    out << "Weight: " << weight << '\n';
		out << "Mesh: " << mesh_name << "\n";
		out << "TimeStep: " << dt << "\n";
		out << "AssemblyTime: " << StefanPME_model.get_assembly_time() << "\n";
		out << "Solving time: " << StefanPME_model.get_solving_time() << "\n";
    out << "AveNewton: " << ave_newton << "\n";
		out << "L2error: " << L2error << "\n";
		out << "Lmp1error: " << Lmp1error << "\n";
		out << "H1error: " << H1error << "\n";
		out << "EnergyError: " << EnergyError << "\n";
		out << "FractionNegValues: " << fraction_neg_values << "\n";
		out << "FractionNegMass: " << fraction_neg_mass << "\n";
		out << "MinValue: " << min_value << "\n";
		out << "MaxValue: " << max_value << "\n";
		out << "MeshSize: " << mesh_ptr->h_max() << "\n";
		out << "NbCells: " << mesh_ptr->n_cells() << "\n";
		out << "NbEdges: " << mesh_ptr->n_edges() << "\n";
		out << "NbEdgeDOFs: " << nbedgedofs << "\n";
		out << "MeshReg: " << meshreg[0] << "\n";
		out << "MeshAnisotropy: " << meshreg[1] << "\n";
		out << std::flush;
    out.close();


}
