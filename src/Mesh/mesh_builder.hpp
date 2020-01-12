// Class to build the mesh data after having read the mesh file
//
// Author: Jerome Droniou (jerome.droniou@monash.edu)
//

#ifndef MESH_BUILDER_HPP
#define MESH_BUILDER_HPP
#include "mesh.hpp"
#include "cell.hpp"
#include "edge.hpp"
#include "vertex.hpp"
#include <string>
#include <vector>
#include <Eigen/Dense>
#include <memory>

#include <stdlib.h>     /* exit, EXIT_FAILURE */

namespace HArDCore2D {

/*!
*  @addtogroup Mesh
* @{
*/

// ----------------------------------------------------------------------------
//                            Class definition
// ----------------------------------------------------------------------------

/// The MeshBuilder class provides build tools to create a full mesh with all connectivities
class MeshBuilder {
public:
    /**
    * Constructor for MeshBuilder.
    */
    MeshBuilder();
    /**
    * Build a Mesh from vertices and cells
    *
    * @param vertices vector containing the coordinates of the vertices ordered
    *using the global ordering. Note that indexes start at 0 in c++
    * @param cells vector containing vectors with  the global indexes of
    *vertices making up a cell
    *
    * @return a pointer to the mesh that was build
    */
    std::unique_ptr<Mesh> build_the_mesh(std::vector<std::vector<double> > vertices,
                           std::vector<std::vector<size_t> > cells);  ///< construct the connectivity in the mesh

private:
    void build_boundary(Mesh* mesh);  ///< identifies boundary cells and vertices, and compute lists of boundary cells, edges and vertices

};

/*@}*/
}
#endif /* MESH_BUILDER_HPP */

