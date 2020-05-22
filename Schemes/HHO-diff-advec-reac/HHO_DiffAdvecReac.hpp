// Implementation of the HHO scheme in 2D for the diffusion equation, with K piecewise constant, \beta \in H^1(\Omega)^d and \mu piecewise continuous
//
//   { -div(K \grad(u) + \beta u) + \mu u = f,       inside Omega
//   {                   K \grad(u) . nTF = g,       on GammaN
//   {                                  u = g,       on GammaD
//
//  Pure Neumann not yet working for nonzero reaction

/*!
 * @defgroup HHO_DiffAdvecReac
 * @brief HHO scheme for a diffusion-advection-reaction scheme.
 */

/** The HHO_DiffAdvecReac code uses a Hybrid High-Order method to solve the equation -div(K \grad(u) + \beta u) + \mu u = f 
 *  with K piecewise constant, \beta \in H^1(\Omega)^d and \mu piecewise continuous. 
 **/

#include "HHO-general/HHO2D.hpp"
#include "boost/program_options.hpp"
#include <unsupported/Eigen/SparseExtra>

const std::string mesh_dir = "../../typ2_meshes/";

/*!
 * @addtogroup HHO_DiffAdvecReac
 * @{
 */

/// A method to set all program paramters from options passed to the main function. Returns true if options are valid, false otherwise.
bool program_options(
    const int,          ///< The number of options passed to the main function
    const char *[],     ///< An array of all options passed to the main function
    std::string &,      ///< The mesh name
    std::string &,      ///< The boundary condition ID
    std::vector<int> &, ///< The test case ID
    size_t &,           ///< Cell polynomial degree
    size_t &,           ///< Edge polynomial degree
    std::string &,      ///< File to plot to
    bool &,             ///< Option to use multi-threading
    bool &              ///< Option to export the system matrix
);

// Parameters that get passed to program options
std::string mesh_name, bc_id, plot_file;
std::vector<int> id_tcase;
size_t L, K;
bool use_threads, export_matrix;

// @}