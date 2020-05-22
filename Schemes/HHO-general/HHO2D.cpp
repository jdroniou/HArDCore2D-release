#include "HHO2D.hpp"

HHO2D::HHO2D(HybridCore &hho, const size_t L, const size_t K, const bool use_threads, size_t doeT, size_t doeF)
    : m_hho(hho),
      m_L(L),
      m_K(K),
      m_use_threads(use_threads),
      m_doeT(doeT),
      m_doeF(doeF)
{   
    // Resize GlobMat and ScBeMat for static condesation
    GlobMat.resize(n_total_edge_dofs, n_total_edge_dofs);
    ScBeMat.resize(n_total_cell_dofs, n_total_edge_dofs);

    // Resize AT for energy norm
    AT.resize(n_cells);

    // If default given, set quadrature degrees
    if (m_doeT == 0)
        m_doeT = m_L + m_K + 1;
    if (m_doeF == 0)
        m_doeF = 2 * m_K + 1;
}

// Set the global operator
void HHO2D::set_global_operator(const MatrixFType &AT)
{
    global_operator = [&](Cell *cell, ElementQuad &elquad) -> Eigen::MatrixXd {
        return AT(cell, elquad);
    };
}

// Set the load vector
void HHO2D::set_load_vector(const VectorFType &bT)
{
    load_vector = [&](Cell *cell, ElementQuad &elquad) -> Eigen::VectorXd {
        return bT(cell, elquad);
    };
}

// Plot exact and numerical solutions to plot_file
void HHO2D::plot(const std::string plot_file, const UVector &sol, const FType<double> &exact_sol)
{
    // Plot if plot_file given and cells less than 5000
    if (plot_file != "" && n_cells < 5000)
    {
        VectorRd coord;

        std::vector<Vertex *> vertices = mesh_ptr->get_vertices();
        size_t n_verts = mesh_ptr->n_vertices();
        Eigen::VectorXd interpolant(n_verts);

        // Calculate interpolant of exact solution
        for (size_t i = 0; i < n_verts; i++)
        {
            coord = vertices[i]->coords();
            interpolant(i) = exact_sol(coord);
        }

        // Mesh mesh = *mesh_ptr;
        VtuWriter plotdata(mesh_ptr);        

        // Calculate cell and face values on the vertices
        Eigen::VectorXd approx_sol_vertex_byF = m_hho.VertexValues(sol, "edge");
        Eigen::VectorXd approx_sol_vertex_byT = m_hho.VertexValues(sol, "cell");       

        // Plot solutions. T == cell, F == face, E == exact
        plotdata.write_to_vtu("T-" + plot_file + ".vtu", approx_sol_vertex_byF);
        plotdata.write_to_vtu("F-" + plot_file + ".vtu", approx_sol_vertex_byT);
        plotdata.write_to_vtu("E-" + plot_file + ".vtu", interpolant);
    }
}

// The standard load vector with no Neumann BCs
VectorFType HHO2D::standard_load_vector(const CellFType<double> &source)
{
    VectorFType bT = [&](Cell *cell, ElementQuad &elquad) -> Eigen::VectorXd {
        // Get the source term on the cell
        FType<double> source_on_cell = [&](const VectorRd x) -> double { return source(x, cell); };

        // Get number of local edges and edge_dofs
        const size_t n_local_edges = cell->n_edges();
        const size_t n_local_dofs = n_local_cell_dofs + n_local_edges * n_local_edge_dofs;

        // get quadrature rule and basis function on quadrature nodes
        QuadratureRule quadT = elquad.get_quadT();
        BasisQuad<double> phiT_quadT = elquad.get_phiT_quadT();

        Eigen::VectorXd source_vec = Eigen::VectorXd::Zero(n_local_dofs);

        // Integrate source term with basis functions
        source_vec.head(n_local_cell_dofs) = integrate(source_on_cell, phiT_quadT, quadT, n_local_cell_dofs);
        return source_vec;
    };

    return bT;
}

// The standard load vector when exact solution is known
VectorFType HHO2D::standard_load_vector(const CellFType<double> &source, const CellFType<VectorRd> &f, const BoundaryConditions &BC)
{
    VectorFType bT = [&](Cell *cell, ElementQuad &elquad) -> Eigen::VectorXd {
        // Get the source term on the cell
        FType<double> source_on_cell = [&](const VectorRd x) -> double { return source(x, cell); };

        // Get number of local edges and edge_dofs
        const size_t n_local_edges = cell->n_edges();
        const size_t n_local_dofs = n_local_cell_dofs + n_local_edges * n_local_edge_dofs;

        // get quadrature rule and basis function on quadrature nodes
        QuadratureRule quadT = elquad.get_quadT();
        BasisQuad<double> phiT_quadT = elquad.get_phiT_quadT();

        Eigen::VectorXd source_vec = Eigen::VectorXd::Zero(n_local_dofs);

        // Integrate source term with basis functions
        source_vec.head(n_local_cell_dofs) = integrate(source_on_cell, phiT_quadT, quadT, n_local_cell_dofs);

        // If cell is boundary, might have Neumann edges
        if (cell->is_boundary())
        {   
            // loop over edges and check for Neumann edges
            for (size_t iTF = 0; iTF < n_local_edges; iTF++)
            {
                Edge *edge = cell->edge(iTF);
                if (BC.type(*edge) == "neu")
                {
                    const size_t offset_F = n_local_cell_dofs + iTF * n_local_edge_dofs;
                    const size_t iF = edge->global_index();
                    const VectorRd &nTF = cell->edge_normal(iTF);

                    // Need higher order quadrature
                    QuadratureRule quadF = generate_quadrature_rule(*edge, 2 * m_K + 2);
                    BasisQuad<double> phiF_quadF = evaluate_quad<Function>::compute(m_hho.EdgeBasis(iF), quadF);

                    // Get Neumann BC
                    FType<double> f_dot_nTF = [&](const VectorRd x) -> double {
                        return nTF.dot(f(x, cell));
                    };

                    // Integrate Neumann BC condition with basis function and add to source vector
                    source_vec.segment(offset_F, n_local_edge_dofs) += integrate(f_dot_nTF, phiF_quadF, quadF);
                }
            }
        }
        return source_vec;
    };

    return bT;
}

// The standard load vector when exact solution is unknown but Neumann BCs given
VectorFType HHO2D::standard_load_vector(const CellFType<double> &source, const FType<double> &f, const BoundaryConditions &BC)
{
    VectorFType bT = [&](Cell *cell, ElementQuad &elquad) -> Eigen::VectorXd {
        // Get the source term on the cell
        FType<double> source_on_cell = [&](const VectorRd x) -> double { return source(x, cell); };

        // Get number of local edges and edge_dofs
        const size_t n_local_edges = cell->n_edges();
        const size_t n_local_dofs = n_local_cell_dofs + n_local_edges * n_local_edge_dofs;

        // get quadrature rule and basis function on quadrature nodes
        QuadratureRule quadT = elquad.get_quadT();
        BasisQuad<double> phiT_quadT = elquad.get_phiT_quadT();

        Eigen::VectorXd source_vec = Eigen::VectorXd::Zero(n_local_dofs);

        // Integrate source term with basis functions
        source_vec.head(n_local_cell_dofs) = integrate(source_on_cell, phiT_quadT, quadT, n_local_cell_dofs);

        // If cell is boundary, might have Neumann edges
        if (cell->is_boundary())
        {   
            // loop over edges and check for Neumann edges
            for (size_t iTF = 0; iTF < n_local_edges; iTF++)
            {
                Edge *edge = cell->edge(iTF);
                if (BC.type(*edge) == "neu")
                {
                    const size_t offset_F = n_local_cell_dofs + iTF * n_local_edge_dofs;
                    const size_t iF = edge->global_index();

                    // Need higher order quadrature
                    QuadratureRule quadF = generate_quadrature_rule(*edge, 2 * m_K + 2);
                    BasisQuad<double> phiF_quadF = evaluate_quad<Function>::compute(m_hho.EdgeBasis(iF), quadF);

                    // Integrate Neumann BC condition with basis function and add to source vector
                    source_vec.segment(offset_F, n_local_edge_dofs) += integrate(f, phiF_quadF, quadF);
                }
            }
        }
        return source_vec;
    };

    return bT;
}

// Set the Dirichlet BC
void HHO2D::set_dirichlet(const FType<double> &f, const size_t n_dir_edges)
{
    UDir = Eigen::VectorXd::Zero(n_local_edge_dofs * n_dir_edges);
    for (size_t idF = 0; idF < n_dir_edges; idF++)
    {
        Edge *edge = mesh_ptr->edge(n_edges - n_dir_edges + idF);
        const size_t iF = edge->global_index();
        QuadratureRule quadF = generate_quadrature_rule(*edge, 2 * m_K + 2);
        BasisQuad<double> phiF_quadF = evaluate_quad<Function>::compute(m_hho.EdgeBasis(iF), quadF);
        UDir.segment(idF * n_local_edge_dofs, n_local_edge_dofs) = l2_projection<HybridCore::PolyEdgeBasisType>(f, m_hho.EdgeBasis(iF), quadF, phiF_quadF);
    }
}

// Set zero Dirichlet BC
void HHO2D::set_dirichlet(const size_t n_dir_edges)
{
    UDir = Eigen::VectorXd::Zero(n_local_edge_dofs * n_dir_edges);
}

void HHO2D::assemble()
{   
    // Start the timer
    boost::timer::cpu_timer timer;
    timer.start();

    // Set up triplets for sparse matrix initialisation
    std::vector<Eigen::Triplet<double>> triplets_GlobMat;
    std::vector<Eigen::Triplet<double>> triplets_ScBe;

    // set up vectors of local matrices to emplace into global matrices
    std::vector<Eigen::MatrixXd> invATT_ATF(n_cells);
    std::vector<Eigen::VectorXd> cell_source(n_cells);
    std::vector<Eigen::MatrixXd> MatF(n_cells);

    // mesh cells
    std::vector<Cell *> cells = mesh_ptr->get_cells();

    // Construct the local matrices using multithreading is use_threads is true
    std::function<void(size_t, size_t)> construct_all_local_contributions = [&](size_t start, size_t end) -> void {
        for (size_t iT = start; iT < end; iT++)
        {
            Cell *cell = cells[iT];
            const size_t n_local_edges = cell->n_edges();
            const size_t n_edge_dofs = n_local_edges * n_local_edge_dofs;

            ElementQuad elquad(m_hho, iT, m_doeT, m_doeF);

            Eigen::VectorXd bT = load_vector(cell, elquad);
            AT[iT] = global_operator(cell, elquad);

            // Get each contribution of AT
            Eigen::MatrixXd ATT = AT[iT].topLeftCorner(n_local_cell_dofs, n_local_cell_dofs);
            Eigen::MatrixXd ATF = AT[iT].topRightCorner(n_local_cell_dofs, n_edge_dofs);
            Eigen::MatrixXd AFT = AT[iT].bottomLeftCorner(n_edge_dofs, n_local_cell_dofs);
            Eigen::MatrixXd AFF = AT[iT].bottomRightCorner(n_edge_dofs, n_edge_dofs);

            Eigen::PartialPivLU<Eigen::MatrixXd> invATT;
            invATT.compute(ATT);
            
            Eigen::VectorXd invATT_bTcell = invATT.solve(bT.head(n_local_cell_dofs));

            // Store the local matrices
            invATT_ATF[iT] = invATT.solve(ATF);
            MatF[iT] = AFF - AFT * invATT_ATF[iT];
            cell_source[iT] = bT.tail(n_edge_dofs) - AFT * invATT_bTcell;

            ScRHS.segment(iT * n_local_cell_dofs, n_local_cell_dofs) = invATT_bTcell;
        }
    };

    // Running the local constructions in parallel
    parallel_for(n_cells, construct_all_local_contributions, m_use_threads);

    // Set the local matrices into the correct positions of the global matrices
    for (size_t iT = 0; iT < n_cells; iT++)
    {
        Cell *cell = cells[iT];
        const size_t n_local_edges = cell->n_edges();
        for (size_t iTF = 0; iTF < n_local_edges; iTF++)
        {
            const size_t iF = cell->edge(iTF)->global_index();
            for (size_t ik = 0; ik < n_local_edge_dofs; ik++)
            {
                const size_t iLocal = iTF * n_local_edge_dofs + ik;
                const size_t iGlobal = iF * n_local_edge_dofs + ik;
                GlobRHS(iGlobal) += cell_source[iT](iLocal);
                for (size_t i = 0; i < n_local_cell_dofs; i++)
                {
                    triplets_ScBe.emplace_back(iT * n_local_cell_dofs + i, iGlobal, invATT_ATF[iT](i, iLocal));
                }
                for (size_t jTF = 0; jTF < n_local_edges; jTF++)
                {
                    const size_t jF = cell->edge(jTF)->global_index();
                    for (size_t jk = 0; jk < n_local_edge_dofs; jk++)
                    {
                        const size_t jLocal = jTF * n_local_edge_dofs + jk;
                        const size_t jGlobal = jF * n_local_edge_dofs + jk;
                        triplets_GlobMat.emplace_back(iGlobal, jGlobal, MatF[iT](iLocal, jLocal));
                    }
                }
            }
        }
    }

    // Construct the global matrices
    GlobMat.setFromTriplets(std::begin(triplets_GlobMat), std::end(triplets_GlobMat));
    ScBeMat.setFromTriplets(std::begin(triplets_ScBe), std::end(triplets_ScBe));

    // Get assembly time
    assembly_time = timer.elapsed().wall;
}

UVector HHO2D::solve()
{
    /*
        * GlobMat = AFF - AFT * INV_ATT * ATF
        * GlobRHS = BF - AFT * INV_ATT_BT
        * 
        * ScBeMat = INV_ATT * ATF 
        * ScRHS   = INV_ATT_BT
    */

    // Start the timer
    boost::timer::cpu_timer timer;
    timer.start();

    const size_t n_fixed_dofs = UDir.size();
    const size_t n_unknowns = n_total_edge_dofs - n_fixed_dofs;

    SysMat = GlobMat.topLeftCorner(n_unknowns, n_unknowns);

    // (D == Dirichlet, I == Internal/Neumann, T == Cell, F == Face)
    //
    //                = BI - AIT * INV_ATT * BT        - (AID - AIT * INV_ATT * ATD * UD)
    Eigen::VectorXd B = GlobRHS.segment(0, n_unknowns) - GlobMat.topRightCorner(n_unknowns, n_fixed_dofs) * UDir;

    // Solve the statically condesed system using BiCGSTAB
    Eigen::BiCGSTAB<Eigen::SparseMatrix<double>> solver;
    solver.compute(SysMat);

    // xF = INV( AII - AIT * INV_ATT * ATI ) * ( BI - AIT * INV_ATT * BT - (AID - AIT * INV_ATT * ATD * UD) )
    Eigen::VectorXd xF = solver.solve(B);

    // Print solver iterations and estimated error
    std::cout << "     [solver] #iterations: " << solver.iterations() << ", estimated error: " << solver.error() << std::endl;

    // Find the residual of the system
    solving_error = (SysMat * xF - B).norm();

    // Set each part of the solution
    Eigen::VectorXd Xh = Eigen::VectorXd::Zero(n_total_cell_dofs + n_total_edge_dofs);
    Xh.tail(n_fixed_dofs) = UDir;
    Xh.segment(n_total_cell_dofs, n_unknowns) = xF;

    //                         = INV_ATT_BT - INV_ATT * ATF * UF
    Xh.head(n_total_cell_dofs) = ScRHS - ScBeMat * Xh.tail(n_total_edge_dofs);

    // Get solving time
    solving_time = timer.elapsed().wall;

    return UVector(Xh, *mesh_ptr, m_L, m_K);
}

// UVector HHO2D::neumann_solve()
// {

// }

double HHO2D::energy_norm(const UVector Xh)
{
    double value = 0.0;

    for (size_t iT = 0; iT < n_cells; iT++)
    {
        Eigen::VectorXd XTF = Xh.restr(iT);
        // a(v, v)
        value += XTF.transpose() * AT[iT] * XTF;
    }

    // a(v, v)^(1/2)
    return sqrt(value);
}