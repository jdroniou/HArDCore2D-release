// Classes to assist with mesh handling and coarsening

#include <iostream>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <functional>
#include <random>
#include <fstream>
#include <Eigen/Dense>

#ifndef _DIRECTEDGRAPH_HPP
#define _DIRECTEDGRAPH_HPP

/*!	
 * @defgroup DirectedGraph 
 * @brief The DirectedGraph class as well as classes that make it up and the NodeArray class
 */

/*!
 * \addtogroup DirectedGraph
 * @{
 */

// ----------------------------------------------------------------------------
//                            Free functions
// ----------------------------------------------------------------------------
 
/// Type definition for a matrix of unsigned ints
typedef std::vector<std::vector<std::size_t>> Array;

// ----------------------------------------------------------------------------
//                            DirectedEdge class definition
// ----------------------------------------------------------------------------

/** The DirectedEdge class is the defining element of a DirectedCell. It takes its two end points as parameters and
	contains functions to assist with comparing edges.
 */

class DirectedEdge
{
public:
	/// Node number of the edge tail
    std::size_t a;
    
    /// Node number of the edge head
    std::size_t b;
	
	/// Null constructor
    DirectedEdge();
    
    /// Destructor
    ~DirectedEdge();
    
    /**@brief Class Constructor: Initialises the edge with the two nodes that form the edge
     **/
    DirectedEdge(size_t, size_t);
	
	/// Returns the edge consisting of the same nodes, pointing in the opposite direction
    DirectedEdge anti_edge();
};

/// Boolean operation to test if two edges are equal
bool operator==(const DirectedEdge &, const DirectedEdge &);

// ----------------------------------------------------------------------------
//                            DirectedCell class definition
// ----------------------------------------------------------------------------

/** The DirectedCell class is the defining element of a DirectedGraph. It is defined by a vector of DirectedEdge's,
	and a vector of unsigned int's describing the fine cells that make up the cell upon coarsening. The class contains functions
	related to cell manipulation.
 **/

class DirectedCell
{
public:
	/// The edges the cell consists of
    std::vector<DirectedEdge> T;
    
    /// The cell ID's of the cell's that have formed this cell
    std::vector<std::size_t> part;
    
    /// Null Constructor
    DirectedCell();
	
	/// Constructor
    DirectedCell(std::size_t);
    
    /// Destructor
    ~DirectedCell();
    
    /**@brief Class Constructor: Initialises the cell with a vector of edges and the cell ID
     **/
    DirectedCell(
    	std::vector<DirectedEdge> &, ///< A reference to the edges that the cell consists of
    	std::size_t ///< Initial cell ID
    );
	
	/// Method to append an edge to the end of the cell
    void add_edge(
    	DirectedEdge & ///< Edge to be added
    );
	
	/// Remove edge at a given position
    void remove_edge(
    	std::size_t ///< Position of edge to be removed
    );
	
	/** Appends the vector of edges and partition vector of a cell to the
	  *	edges and partition of this cell
	 **/
    void append_cell(
    	DirectedCell & ///< Cell to be added
    );
    
    /// Checks for repeated or lone nodes
    bool check_nodes();
	
	/// Tests if the edges of the cell are ordered
    bool is_ordered();
	
	/// Method to remove all edge - antiedge pairs
    void remove_duplicate_edges();
	
	/// Orders the cell edges - returns vector of edges to be removed
    std::vector<DirectedEdge> order_edges();
    
    /// Tests if cell has an edge
    bool has_edge(
    	DirectedEdge ///< Edge to check for
    );
	
	/// Prints all the edges of the cell
    void print(); // Useful for debugging
};

// ----------------------------------------------------------------------------
//                            NodeArray class definition
// ----------------------------------------------------------------------------

/** The NodeArray class is generated from the directed graph, and combined with the vertex coordinates
  *	fully describes a mesh.
 **/
 
class NodeArray
{
public:
	/// The cell-node array
    Array A;
	
	/// Null constructor
    NodeArray();
    
    /// Default constructor
    ~NodeArray();
    
    /**@brief Class Constructor: Initialises the NodeArray with the cell-node array
     **/
    NodeArray(
    	Array & ///< The cell-node array
    );
    
    /// Test if a given node exists in the cell-node array
    bool node_exists(
    	std::size_t ///< Node to test for
    );
    
    /// Subtracts one from all nodes greater than a given node
    void renum_nodes(
    	std::size_t ///< Node given
    );
    
    /// Prints the cell node array to an out file stream
    void print(
    	std::ofstream *, ///< Pointer to the file stream to print to 
    	std::size_t width = 6 ///< Optional parameter to format output
    );
};

class DirectedGraph
{
public:
	/// Vector of cells that form the graph
    std::vector<DirectedCell> G;
    
    /// Null constructor
    DirectedGraph();
    
    /// Default constructor
    ~DirectedGraph();
    
    /**@brief Class Constructor: Initialises the DirectedGraph with the cells
     **/
    DirectedGraph(
    	std::vector<DirectedCell> & ///< Initial cells
    );
	
	/// Appends a cell to the end of the graph
    void add_cell(
    	DirectedCell & ///< Cell to be added
    );

	/// Removes cell from a given position
    void remove_cell(
    	std::size_t ///< Cell to be removed
    );
	
	/// Test if graph has duplicate edges or unordered cells
    bool test_graph();    
	
	/// Randomise order of cells in graph
    void randomise();
	
	/// Order cells in graph by the number of edges in each cell
    void order();
	
	/// Coarsen the graph by merging cells
    void coarsen();
	
	/// Returns the cell-node array the graph corresponds to
    NodeArray graph_to_array();
    
    /// Returns list of fine cells that each coarse cell consists of
    std::string get_partition();
    
    /// Outputs to a .dat file for gnuplot to read
    void plotfile(
    	std::ofstream *, ///< Pointer to the file stream to print to 
    	std::vector<Eigen::VectorXd> & ///< Vertices
    );

private:
	/* Checks if there exists a shared edge, and checks if there exists
	   a shared node at a non-edge interface i.e. they touch at a corner.
	   Must be true, false to return true. Also stores a vector of shared 
	   edges to pass to merge_cells - which deletes them.
	 */
	bool can_merge(std::size_t, std::size_t, std::vector<size_t> &);
	
	/* Appends second cell to the first. Removes second cell. Removes 
	   duplicate edges of first cell then orders edges. Handles
	   non-simply connected case.
	 */
    void merge_cells(std::size_t, std::size_t, std::vector<size_t> &);    
    
    // Random boolean generator for coarsen()
    std::knuth_b rand_engine;
    bool random_bool(double);
};

//@}

#endif
