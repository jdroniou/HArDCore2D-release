// Creates a .vtu file of the solution, to be visualised by paraview
//
// Author: Jerome Droniou (jerome.droniou@monash.edu)
//

/*
*
*  This library was developed around HHO methods, although some parts of it have a more
* general purpose. If you use this code or part of it in a scientific publication, 
* please mention the following book as a reference for the underlying principles
* of HHO schemes:
*
* The Hybrid High-Order Method for Polytopal Meshes: Design, Analysis, and Applications.
* D. A. Di Pietro and J. Droniou. 2019, 516p. 
* url: https://hal.archives-ouvertes.fr/hal-02151813.
*
*/



#ifndef VTU_WRITER_HPP
#define VTU_WRITER_HPP
#include <mesh.hpp>
#include <cstring>
#include <vector>
#include <Eigen/Dense>
using Eigen::Vector2d;

/*!  
*  @defgroup Plot 
* @brief Classes providing tools to create vtu files to visualise solutions
*/

namespace HArDCore2D {

/*!
*  @addtogroup Plot
* @{
*/
class VtuWriter {
public:
    /**
    * @brief Constructor for mesh writer 
    *
    * @param mesh pointer to the mesh
    */
    VtuWriter(Mesh* mesh); 

    /// Writes the vtu file  
    bool write_to_vtu(
          std::string file_name, ///< name of file to write to
          Eigen::VectorXd sol_vertex,   ///< vector of values of the solution at the mesh vertices
          bool dimen = true  ///< elevation of the plot: planar if dimen=0, elevated if dimen=1
          );

    /// overloaded writer for the mesh alone
    bool write_to_vtu(std::string file_name); 

private:
    void write_vertices(FILE* pFile);/// <\brief add vertices coords to the vtk file 
    void write_vertices(FILE* pFile, Eigen::VectorXd data);/// <\brief add vertices coords to the vtk file for 3d graphs
    void write_header(FILE* pFile); /// <\brief add header to the file
    void write_cells(FILE* pFile); /// <\brief add the cell data - note we just use a polygon form vtk type 7
    /**
    * @brief     *
    * @param pFile
    * @param alldata
    * @param names
    */
   
    void write_point_scalar_property(FILE* pFile, Eigen::VectorXd data,
                                     std::string name); /// <\brief writes all of the point data to the vtk file
    void write_footer(FILE* pFile); /// <\brief add the footer to the vtk file

    Mesh* _mesh;

};

//@}
}  // namespace HArDCore2D

#endif /* VTU_WRITER_HPP */
