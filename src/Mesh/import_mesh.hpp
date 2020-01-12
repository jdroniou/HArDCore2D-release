// Class to read a typ2 mesh file
//
// Author: Jerome Droniou (jerome.droniou@monash.edu)
//

#ifndef IMPORT_MESH_HPP
#define IMPORT_MESH_HPP

#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
namespace HArDCore2D {
/**
* Flag for the importer to know what to do with the text file
*/
enum Flag {
    header_,
    nv_,
    vertices_,
    nc_,
    cells_,
    centers_
};

/*!
*  @addtogroup Mesh
* @{
*/

// ----------------------------------------------------------------------------
//                            Class definition
// ----------------------------------------------------------------------------

/// The MeshReaderTyp2 class provides functions to read a .typ2 mesh file
class MeshReaderTyp2 {
public:
    /**
    * Constructor for mesh reader
    *
    * @param file_name name of the file name, needs to include the full path
    */
    MeshReaderTyp2(std::string file_name);  ///< class to read the cells and vertices in a .typ2 file
    /**
    * Reads the file into the specified containers
    *
    * @param vertices reference to a vector to hold the vertices coordinates
    * @param cells reference to a vector to hold the cell indexes
    * @param centers reference to a vector to hold the cell centers coordinates
    */
    bool read_mesh(std::vector<std::vector<double> >& vertices,
                   std::vector<std::vector<size_t> >& cells,
                   std::vector<std::vector<double> >& centers);  ///< reads the .typ2 file and fills in cells, vertices and centers

private:
    std::string _file_name;  ///< name of the file being read
};
/*@}*/

};     // end namespace HArDCore2D
#endif /* IMPORT_MESH_HPP */
