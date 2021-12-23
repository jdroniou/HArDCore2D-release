#include "DirectedGraph.hpp"
#include <string>
#include <boost/timer/timer.hpp>

//const std::string mesh_root = "/home/liam/HHO/liamyemm/HArDCore2D-liam/typ2_meshes/";
const std::string mesh_root = "../../typ2_meshes/";

int main(int argc, char *argv[])
{
    std::string mesh_file = argv[1];
    std::size_t iterations = std::stoi(argv[2]);
    
    std::cout.unsetf ( std::ios::floatfield );
    std::cout.precision(4);

    //  Open mesh file

    std::ifstream inMesh;
    inMesh.open(mesh_root + mesh_file + ".typ2");
    if (!inMesh)
    {
        std::cerr << "Unable to open mesh file\n";
        exit(1);
    }

    std::cout << "Reading mesh file...\n";

    //  Read mesh file

    std::string ignore_line;
    std::size_t n_verts;
    std::size_t n_cells;
    std::vector<Eigen::VectorXd> vertices;

    std::getline(inMesh, ignore_line);
    inMesh >> n_verts;
    std::size_t count = 0;
    double xcoord, ycoord;
    Eigen::VectorXd coord(2);

    while (count < n_verts && inMesh >> xcoord >> ycoord)
    {
        coord << xcoord, ycoord;
        vertices.push_back(coord);
        count++;
    }

    std::cout << "Read " + std::to_string(vertices.size()) + "/" + std::to_string(n_verts) + " vertices\n";

    std::getline(inMesh, ignore_line);
    std::getline(inMesh, ignore_line);
    inMesh >> n_cells;

    DirectedGraph graph;

    std::size_t cell_count = 0;
    std::size_t n_local_cell_nodes;

    while (cell_count < n_cells && inMesh >> n_local_cell_nodes)
    {
        std::size_t node_count = 0;
        std::vector<size_t> cell_nodes(n_local_cell_nodes);
        while (node_count < n_local_cell_nodes && inMesh >> cell_nodes[node_count])
        {
            node_count++;
        }

        DirectedCell cell(cell_count + 1);
        for (std::size_t iV = 0; iV < n_local_cell_nodes - 1; iV++)
        {
            DirectedEdge edge(cell_nodes[iV], cell_nodes[iV + 1]);
            cell.add_edge(edge);
        }
        DirectedEdge edge(cell_nodes[n_local_cell_nodes - 1], cell_nodes[0]);
        cell.add_edge(edge);
        graph.add_cell(cell);
        cell_count++;
    }

    std::cout << "Read " + std::to_string(graph.G.size()) + "/" + std::to_string(n_cells) + " cells\n";

    inMesh.close();

    if (!graph.test_graph())
    {
        exit(1);
    }

    std::cout << "\nCoarsening graph...\n";
    graph.randomise();
    for (size_t i = 0; i < iterations; i++)
    {
    	boost::timer::cpu_timer timer;
    	timer.start();        
        graph.coarsen();
        double time = double(timer.elapsed().wall) * pow(10, -9);
		std::cout << std::setprecision(4) << "Iteration " + std::to_string(i+1) + " completed in " << time << "s. " + std::to_string(graph.G.size()) + " cells remaining\n";
		graph.order();
    }
    std::cout << "Coarsening complete\n";

    if (!graph.test_graph())
    {
        exit(1);
    }   
    
    // Plot the mesh
    // std::ofstream plotout(mesh_root + "agglomerated/" + mesh_file + ".coarse." + std::to_string(iterations) + ".dat");
    // graph.plotfile(&plotout, vertices);
    // plotout.close();    
    
    std::string partition = graph.get_partition();
    
    NodeArray cell_node_array = graph.graph_to_array();

    // Remove unused vertices
    for (size_t iV = n_verts - 1; iV > 0; --iV)
    {
        if (!(cell_node_array.node_exists(iV + 1)))
        {
            vertices.erase(vertices.begin() + iV);
            cell_node_array.renum_nodes(iV + 1);
            n_verts--;
        }
    }

    std::cout << "\nWriting mesh to file\n";

    std::ofstream out(mesh_root + "agglomerated/" + mesh_file + ".coarse." + 			std::to_string(iterations) + ".typ2");

    out << "Vertices"
        << "\n";
    out << n_verts << "\n";
    for (size_t iV = 0; iV < n_verts; iV++)
    {
        out << std::setprecision(16) << vertices[iV](0);
        out << " ";
        out << std::setprecision(16) << vertices[iV](1);
        out << "\n";
    }
    
    cell_node_array.print(&out);
    out << partition;
    out.close();

    return 0;
}
