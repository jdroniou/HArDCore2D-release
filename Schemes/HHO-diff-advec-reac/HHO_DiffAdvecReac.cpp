#include "HHO_DiffAdvecReac.hpp"

int main(const int argc, const char *argv[])
{
    // Set options
    if (!program_options(argc, argv, mesh_name, bc_id, id_tcase, L, K, plot_file, use_threads, export_matrix))
    {
        return 0;
    }

    // Build mesh and reorder edges
    const std::string mesh_file = mesh_dir + mesh_name + ".typ2";
    MeshBuilder builder = MeshBuilder(mesh_file);
    std::unique_ptr<Mesh> mesh_ptr = builder.build_the_mesh();

    BoundaryConditions BC(bc_id, *mesh_ptr.get());
    BC.reorder_edges();

    // Load the test case
    TestCase tcase(id_tcase);

    // Get flags to speed up assembly time in the case of zero or constant advection or reaction
    bool advec_zero = (id_tcase[2] == 0 ? true : false);
    bool reac_zero = (id_tcase[3] == 0 ? true : false);
    bool reac_const = tcase.is_reac_const();
    bool div_advec_zero = tcase.is_div_advec_zero();
    bool div_advec_const = tcase.is_div_advec_const();

    // Exact solution
    const FType<double> exact_solution = tcase.sol();
    const CellFType<VectorRd> grad_exact_solution = tcase.grad_sol();

    // Diffusion
    const CellFType<MatrixRd> diff = tcase.diff();

    // Advection
    const CellFType<VectorRd> advec = tcase.advec();
    const CellFType<double> div_advec = tcase.div_advec();

    // Reaction
    const CellFType<double> reac = tcase.reac();

    // Source term
    const CellFType<double> source = ((advec_zero && reac_zero) ? tcase.diff_source() : tcase.diff_advec_reac_source());

    // Load hybrid core
    HybridCore hho = HybridCore(mesh_ptr.get(), K + 1, K, use_threads);

    // Load model
    std::cout << "\n[Scheme] Loading model\n";
    HHO2D model(hho, L, K, use_threads);

    // Print the data
    std::cout << "\n[Scheme] Data:\n";
    std::cout << "     Boundary conditions: " << BC.name() << "\n";
    std::cout << "     Test case: solution = " << id_tcase[0] << "; diffusion = " << id_tcase[1] << "; advection = " << id_tcase[2] << "; reaction = " << id_tcase[3] << "\n";
    std::cout << "     Mesh = " << mesh_name << ", nb cells= " << mesh_ptr->n_cells() << ", nb edges = " << mesh_ptr->n_edges() << "\n";
    std::cout << "     Degrees: edge = " << K << "; cell = " << L << "\n";
    std::cout << "     Using threads = " << (use_threads ? "true" : "false") << "; Export matrix = " << (export_matrix ? "true" : "false") << "\n";

    const size_t n_local_edge_dofs = DimPoly<Edge>(K);
    const size_t n_local_highorder_dofs = DimPoly<Cell>(K + 1);
    const size_t n_local_cell_dofs = DimPoly<Cell>(L);

    const size_t d = mesh_ptr->dim();

    // Calculate the global operator
    const MatrixFType AT = [&](Cell *cell, ElementQuad &elquad) -> Eigen::MatrixXd {
        const size_t n_local_edges = cell->n_edges();
        const size_t n_local_dofs = n_local_cell_dofs + n_local_edges * n_local_edge_dofs;

        // HybridCore::PolyCellBasisType phiT = hho.CellBasis(cell->global_index());

        // Declare RHS of G_T and P_T matrices
        Eigen::MatrixXd BG = Eigen::MatrixXd::Zero(n_local_cell_dofs, n_local_dofs);
        Eigen::MatrixXd BP = Eigen::MatrixXd::Zero(n_local_highorder_dofs, n_local_dofs);

        // Get quadratures
        QuadratureRule quadT = elquad.get_quadT();
        BasisQuad<double> phiT_quadT = elquad.get_phiT_quadT();
        BasisQuad<VectorRd> dphiT_quadT = elquad.get_dphiT_quadT();

        // cell mass matrix
        Eigen::MatrixXd MTT = compute_gram_matrix(phiT_quadT, phiT_quadT, quadT, n_local_cell_dofs, n_local_highorder_dofs, "sym");

        // Get L_T matrix to impose average condition on P_T later
        Eigen::VectorXd LT = (MTT.row(0)).transpose();
        Eigen::MatrixXd LTtLT = LT * (LT.transpose());

        // Calculate stiffness matrix
        FType<MatrixRd> diff_on_cell = [&](const VectorRd x) -> MatrixRd { return diff(cell->center_mass(), cell); };
        Eigen::MatrixXd ST = compute_weighted_gram_matrix(diff_on_cell, dphiT_quadT, dphiT_quadT, quadT, "sym");
        double scalT = ST.norm() / LTtLT.norm();

        // face-face and face-cell mass matrices required for the stabilisation
        std::vector<Eigen::MatrixXd> MFF(n_local_edges, Eigen::MatrixXd::Zero(n_local_edge_dofs, n_local_edge_dofs));
        std::vector<Eigen::MatrixXd> MFT(n_local_edges, Eigen::MatrixXd::Zero(n_local_edge_dofs, n_local_highorder_dofs));

        if (!advec_zero)
        {
            // If nonzero advection, set cell - cell term of RHS of G_T
            FType<VectorRd> advec_on_cell = [&](const VectorRd x) -> VectorRd { return advec(x, cell); };
            BG.topLeftCorner(n_local_cell_dofs, n_local_cell_dofs) = compute_weighted_gram_matrix(advec_on_cell, phiT_quadT, dphiT_quadT, quadT, n_local_cell_dofs, n_local_cell_dofs);
        }

        // Set cell - cell term of RHS of P_T
        BP.topLeftCorner(n_local_highorder_dofs, n_local_cell_dofs) = ST.topLeftCorner(n_local_highorder_dofs, n_local_cell_dofs) + scalT * LTtLT.topLeftCorner(n_local_highorder_dofs, n_local_cell_dofs);

        // Loop over edges to fet face terms of the RHS of PT and GT
        for (size_t iTF = 0; iTF < n_local_edges; iTF++)
        {
            const size_t offset_F = n_local_cell_dofs + iTF * n_local_edge_dofs;
            const VectorRd &nTF = cell->edge_normal(iTF);

            QuadratureRule quadF = elquad.get_quadF(iTF);
            BasisQuad<double> phiF_quadF = elquad.get_phiF_quadF(iTF);
            BasisQuad<double> phiT_quadF = elquad.get_phiT_quadF(iTF);
            BasisQuad<VectorRd> dphiT_quadF = elquad.get_dphiT_quadF(iTF);

            MFF[iTF] = compute_gram_matrix(phiF_quadF, phiF_quadF, quadF, "sym");
            MFT[iTF] = compute_gram_matrix(phiF_quadF, phiT_quadF, quadF);

            if (!advec_zero)
            {
                FType<double> advec_dot_nTF = [&](const VectorRd x) -> double { return advec(x, cell).dot(nTF); };
                BG.topLeftCorner(n_local_cell_dofs, n_local_cell_dofs) -= compute_weighted_gram_matrix(advec_dot_nTF, phiT_quadF, phiT_quadF, quadF, n_local_cell_dofs, n_local_cell_dofs, "sym");
                BG.block(0, offset_F, n_local_cell_dofs, n_local_edge_dofs) = compute_weighted_gram_matrix(advec_dot_nTF, phiT_quadF, phiF_quadF, quadF, n_local_cell_dofs, n_local_edge_dofs);
            }

            FType<VectorRd> diff_nTF = [&](const VectorRd x) -> VectorRd { return diff(x, cell) * nTF; };
            BP.topLeftCorner(n_local_highorder_dofs, n_local_cell_dofs) -= compute_weighted_gram_matrix(diff_nTF, dphiT_quadF, phiT_quadF, quadF, n_local_highorder_dofs, n_local_cell_dofs);
            BP.block(0, offset_F, n_local_highorder_dofs, n_local_edge_dofs) = compute_weighted_gram_matrix(diff_nTF, dphiT_quadF, phiF_quadF, quadF);

            // for(size_t i = 1; i < n_local_highorder_dofs; i++) {
            //     FType<double> diff_dphiT_quadF_dot_nTF = [&](const VectorRd x) -> double {
            //         return nTF.dot(diff(cell->center_mass(), cell) * phiT.gradient(i, x));
            //     };
            //     BP.row(i).segment(0, n_local_cell_dofs) -= integrate(diff_dphiT_quadF_dot_nTF, phiT_quadF, quadF, n_local_cell_dofs).transpose();
            //     BP.row(i).segment(offset_F, n_local_edge_dofs) += integrate(diff_dphiT_quadF_dot_nTF, phiF_quadF, quadF).transpose();
            // }
        }

        // Set PT and the diffusion consistency term
        Eigen::MatrixXd PT = (ST + scalT * LTtLT).ldlt().solve(BP);
        Eigen::MatrixXd ATDIFF = PT.transpose() * ST * PT;

        // Operator to transform an element from the local space of unknowns to an element on a cell
        Eigen::MatrixXd ID_CELL = Eigen::MatrixXd::Zero(n_local_cell_dofs, n_local_dofs);
        ID_CELL.topLeftCorner(n_local_cell_dofs, n_local_cell_dofs) = Eigen::MatrixXd::Identity(n_local_cell_dofs, n_local_cell_dofs);

        Eigen::MatrixXd ATADVEC = Eigen::MatrixXd::Zero(n_local_dofs, n_local_dofs);
        Eigen::MatrixXd ATREAC = Eigen::MatrixXd::Zero(n_local_cell_dofs, n_local_cell_dofs);

        // If nonzero advection, set the advection consistency term
        if (!advec_zero)
        {
            Eigen::MatrixXd GT = MTT.topLeftCorner(n_local_cell_dofs, n_local_cell_dofs).ldlt().solve(BG);
            Eigen::MatrixXd ADVEC = ID_CELL.transpose() * MTT.topLeftCorner(n_local_cell_dofs, n_local_cell_dofs) * GT;
            ATADVEC = 0.5 * (ADVEC - ADVEC.transpose());
            if (!div_advec_zero)
            {
                if (div_advec_const)
                {
                    ATADVEC.topLeftCorner(n_local_cell_dofs, n_local_cell_dofs) += 0.5 * div_advec(cell->center_mass(), cell) * MTT.topLeftCorner(n_local_cell_dofs, n_local_cell_dofs);
                }
                else
                {
                    FType<double> div_advec_on_cell = [&](const VectorRd x) -> double { return div_advec(x, cell); };
                    ATADVEC.topLeftCorner(n_local_cell_dofs, n_local_cell_dofs) += 0.5 * compute_weighted_gram_matrix(div_advec_on_cell, phiT_quadT, phiT_quadT, quadT, n_local_cell_dofs, n_local_cell_dofs, "sym");
                }
            }
        }

        // If nonzero reaction, set the reaction consistency term
        if (!reac_zero)
        {
            if (reac_const)
            {
                ATREAC = reac(cell->center_mass(), cell) * MTT.topLeftCorner(n_local_cell_dofs, n_local_cell_dofs);
            }
            else
            {
                FType<double> reac_on_cell = [&](const VectorRd x) -> double { return reac(x, cell); };
                ATREAC = compute_weighted_gram_matrix(reac_on_cell, phiT_quadT, phiT_quadT, quadT, n_local_cell_dofs, n_local_cell_dofs, "sym");
            }
        }

        Eigen::MatrixXd ATCONS = ATDIFF;
        if (!advec_zero)
        {
            ATCONS += ATADVEC;
        }
        if (!reac_zero)
        {
            ATCONS.topLeftCorner(n_local_cell_dofs, n_local_cell_dofs) += ATREAC;
        }

        // Cell difference operator
        Eigen::MatrixXd DT = MTT.topLeftCorner(n_local_cell_dofs, n_local_cell_dofs).inverse() * MTT * PT - ID_CELL;

        // Declare diffusion and advection stabilisation
        Eigen::MatrixXd DIFFSTAB = Eigen::MatrixXd::Zero(n_local_dofs, n_local_dofs);
        Eigen::MatrixXd ADVECSTAB = Eigen::MatrixXd::Zero(n_local_dofs, n_local_dofs);

        for (size_t iTF = 0; iTF < n_local_edges; iTF++)
        {
            const size_t offset_F = n_local_cell_dofs + iTF * n_local_edge_dofs;
            // Operator to transform an element from the local space of unknowns to an element on a edge
            Eigen::MatrixXd ID_EDGE = Eigen::MatrixXd::Zero(n_local_edge_dofs, n_local_dofs);
            ID_EDGE.block(0, offset_F, n_local_edge_dofs, n_local_edge_dofs) = Eigen::MatrixXd::Identity(n_local_edge_dofs, n_local_edge_dofs);

            const VectorRd xF = cell->edge(iTF)->center_mass();

            const VectorRd &nTF = cell->edge_normal(iTF);
            // diffTF = k_Tn.n
            const double diffTF = (diff(xF, cell) * nTF).dot(nTF);
            const double hf = cell->edge(iTF)->diam();

            Eigen::MatrixXd MFF_INV_MFT = MFF[iTF].ldlt().solve(MFT[iTF]);

            if (!advec_zero)
            {
                // advecTF = |b.n|
                const double advecTF = std::abs(advec(xF, cell).dot(nTF));

                // Advection stabilisation
                Eigen::MatrixXd ID_MINUS_MFF_INV_MFT = ID_EDGE - MFF_INV_MFT.topLeftCorner(n_local_edge_dofs, n_local_cell_dofs) * ID_CELL;
                ADVECSTAB += advecTF * ID_MINUS_MFF_INV_MFT.transpose() * MFF[iTF] * ID_MINUS_MFF_INV_MFT;
            }

            // Diffusion stabilisation
            Eigen::MatrixXd DTF_MINUS_MFF_INV_MFT_DT = MFF_INV_MFT * PT - ID_EDGE - MFF_INV_MFT.topLeftCorner(n_local_edge_dofs, n_local_cell_dofs) * DT;
            DIFFSTAB += (diffTF / hf) * DTF_MINUS_MFF_INV_MFT_DT.transpose() * MFF[iTF] * DTF_MINUS_MFF_INV_MFT_DT;
        }

        Eigen::MatrixXd ATSTAB = d * DIFFSTAB;
        if (!advec_zero)
        {
            ATSTAB += 0.5 * ADVECSTAB;
        }
        // return the consistency term plus the stabilisation term
        return ATCONS + ATSTAB;
    };

    // diff * grad for Neumann Bcs
    CellFType<VectorRd> diff_grad_sol = [&](const VectorRd x, const Cell *cell) -> VectorRd {
        return diff(x, cell) * grad_exact_solution(x, cell);
    };

    // Set standard load vector
    const VectorFType bT = model.standard_load_vector(source, diff_grad_sol, BC);

    // Set operators and vectors
    model.set_global_operator(AT);
    model.set_load_vector(bT);
    model.set_dirichlet(exact_solution, BC.n_dir_edges());

    // Assemble model
    std::cout << "\n[Scheme] Assembling\n";
    model.assemble();
    std::cout << "     Assembly time = " << model.get_assembly_time() << "s\n";

    // Solve model
    std::cout << "\n[Scheme] Solving\n";
    // Once neumann_solve is working
    // UVector sol = (bc_id == "N" && id_tcase[3] == 0 ? model.neumann_solve() : model.solve());
    UVector sol = model.solve();
    std::cout << "     Solving time = " << model.get_solving_time() << "s\n";
    std::cout << "     Residual of the linear system = " << model.get_solving_error() << "\n";

    if (export_matrix)
    {
        std::cout << "\n[Scheme] Exporting matrix to Matrix Market format\n";
        saveMarket(model.get_SysMat(), "SystemMatrix.mtx");
    }

    // Compute errors
    std::cout << "\n[Scheme] Computing error\n";

    std::cout << "     Interpolant calculated\n";
    UVector Uh = hho.interpolate(exact_solution, L, K, 2 * K + 2);

    double L2error = hho.L2norm(sol - Uh) / hho.L2norm(Uh);
    std::cout << "     L2 Error = " << L2error << "\n";

    double H1error = hho.H1norm(sol - Uh) / hho.H1norm(Uh);
    std::cout << "     H1 Error = " << H1error << "\n";

    double energyerror = model.energy_norm(sol - Uh) / model.energy_norm(Uh);
    std::cout << "     Energy Error = " << energyerror << "\n";

    std::cout << "\n[Scheme] Writing solution to file\n";

    // Write to results file
    std::ofstream out("results.txt");
    out << "BC: " << BC.name() << '\n';
    out << "Solution: " << id_tcase[0] << '\n';
    out << "Diffusion: " << id_tcase[1] << '\n';
    out << "Advection: " << id_tcase[2] << '\n';
    out << "Reaction: " << id_tcase[3] << '\n';
    out << "Mesh: " << mesh_name << "\n";
    out << "EdgeDegree: " << K << "\n";
    out << "CellDegree: " << L << "\n";
    out << "L2Error: " << L2error << "\n";
    out << "H1Error: " << H1error << "\n";
    out << "EnergyError: " << energyerror << "\n";
    out << "MeshSize: " << mesh_ptr->h_max() << "\n";
    out << "NbCells: " << mesh_ptr->n_cells() << "\n";
    out << "NbEdges: " << mesh_ptr->n_edges() << "\n";
    out << "MeshReg: " << mesh_ptr->regularity()[0] << "\n";
    out << "MeshSkew: " << mesh_ptr->regularity()[1] << "\n";
    out << std::flush;
    out.close();

    // Plot numerical and exact solutions
    model.plot(plot_file, sol, exact_solution);

    return 0;
}

bool program_options(int argc, const char *argv[], std::string &mesh_name, std::string &bc_id, std::vector<int> &id_tcase, size_t &L, size_t &K, std::string &plot_file, bool &use_threads, bool &export_matrix)
{
    namespace po = boost::program_options;

    // Program options
    po::options_description desc("Allowed options");
    desc.add_options()("help,h", "Produce help message")("mesh,m", po::value<std::string>(), "Set the mesh")("bc,b", po::value<std::string>(), "Set the boundary conditions (D=Dirichlet, N=Neumann, M0=Mixed)")("testcase,c", po::value<std::vector<int>>()->multitoken(), "Set the test case (sol diff advec reac)")("celldegree,l", po::value<size_t>(), "Set the degree of the cell polynomials")("edgedegree,k", po::value<size_t>(), "Set the degree of the edge polynomials")("plot,p", po::value<std::string>(), "Plot to file")("use_threads,u", po::value<bool>(), "Using multithreading")("export_matrix,e", po::value<bool>(), "Export matrix to Matrix Market format");

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    // Display help message
    if (vm.count("help"))
    {
        std::cout << desc << "\n";
        return false;
    }

    // Get mesh file
    mesh_name = (vm.count("mesh") ? vm["mesh"].as<std::string>() : "hexa1_3");

    // Get bc
    bc_id = (vm.count("bc") ? vm["bc"].as<std::string>() : "D");

    // Get test case
    std::vector<int> default_tcase = {1, 1, 1, 1};
    id_tcase = (vm.count("testcase") ? vm["testcase"].as<std::vector<int>>() : default_tcase);

    // Get polynomial degrees
    L = (vm.count("celldegree") ? vm["celldegree"].as<size_t>() : 0);
    K = (vm.count("edgedegree") ? vm["edgedegree"].as<size_t>() : 0);

    // Check compatible edge and cell degrees
    if ((std::abs(int(K) - int(L)) > 1) || (K < 0) || (L < 0))
    {
        std::cout << "Degrees k and l not in acceptable range (k,l non-negative, |k-l| <= 1): k=" << K << ", l=" << L << "\n";
        return false;
    }

    // Get plot file
    plot_file = (vm.count("plot") ? vm["plot"].as<std::string>() : "");

    // Get use_threads
    use_threads = (vm.count("use_threads") ? vm["use_threads"].as<bool>() : true);

    // Get export_matrix
    export_matrix = (vm.count("export_matrix") ? vm["export_matrix"].as<bool>() : false);

    return true;
}