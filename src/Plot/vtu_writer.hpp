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
*  D. A. Di Pietro and J. Droniou. Modeling, Simulation and Applications, vol. 19. 
*  Springer International Publishing, 2020, xxxi + 525p. doi: 10.1007/978-3-030-37203-3. 
*  url: https://hal.archives-ouvertes.fr/hal-02151813.
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
    VtuWriter(const Mesh* mesh); 

    /// Writes the vtu file  
    bool write_to_vtu(
          const std::string & filename, ///< name of file to write to
          const std::vector<Eigen::VectorXd> & sol_vrtx,   ///< Each element in the std::vector is a vector of values of a function to plot at the mesh vertices
  				const std::vector<std::string> & sol_names = {}     ///< each string is the name of one function
          );

    /// Overloaded writer for the mesh alone
    bool write_to_vtu(std::string file_name); 

    /// Overload to simplify the call when only one solution is involved
    inline bool write_to_vtu(
				  const std::string & filename, ///< name of file to write to
				  const Eigen::VectorXd & sol_vrtx,  ///< values of the solution at the mesh vertices
				  const std::string & sol_name = "solution"  ///< name of the solution
				  )
		  {
		    return write_to_vtu(filename, std::vector<Eigen::VectorXd> {sol_vrtx}, std::vector<std::string> {sol_name});
		  }

private:
    void write_vertices(FILE* pFile);/// <\brief add vertices coords to the vtk file 
    void write_vertices(FILE* pFile, const Eigen::VectorXd & data);/// <\brief add vertices coords to the vtk file for 3d graphs
    void write_header(FILE* pFile); /// <\brief add header to the file
    void write_cells(FILE* pFile); /// <\brief add the cell data - note we just use a polygon form vtk type 7
    /**
    * @brief     *
    * @param pFile
    * @param alldata
    * @param names
    */
   
    void write_solution(FILE* pFile, const Eigen::VectorXd & sol, const std::string & name); /// <\brief writes a solution to the vtu file
    void write_footer(FILE* pFile); /// <\brief add the footer to the vtu file

    const Mesh* _mesh;

};

//@}
}  // namespace HArDCore2D

#endif /* VTU_WRITER_HPP */
