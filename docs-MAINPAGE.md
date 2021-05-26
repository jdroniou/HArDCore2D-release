HArD::Core (sources: https://github.com/jdroniou/HArDCore) provides a suite of C++ tools to implement numerical schemes whose unknowns are polynomials in the cells, on the edges, and on the faces. The focus is on dealing on generic polytopal meshes. This document addresses the 2D version of HArD::Core. The tools are similar in 2D and 3D, and we refer to the main page of documentation of the 3D version for a more thorough introduction to the library. We only present here the specific way of loading a 2D mesh structure (which differs from the 3D version), and a brief description of the schemes available in this 2D library.

Transferring a scheme's implementation from 3D to 2D or vice-versa is very straightforward, provided that the scheme's mathematical definition does not depend on the dimension and that the generic types provided in `basis.hpp` are used; see README file of the HArD::Core github depository https://github.com/jdroniou/HArDCore.

\tableofcontents

* [Loading a 2D mesh](#mesh) -- How to load a mesh.
* [Schemes](#schemes) -- The list of schemes currently implemented in HArD::Core2D, and scripts to run them.

<a name="mesh">
\section loading_mesh Loading a mesh
</a>

HArDCore2D currently reads meshes in the `typ2` format designed for the <a href="https://www.i2m.univ-amu.fr/fvca5/benchmark/index.html">FVCA5 Benchmark</a>. A short documentation describing this format is provided in the `typ2_meshes` directory (see README.pdf therein). Several meshes can also be found in this directory.

A mesh file must be read using an instance of the `MeshReaderTyp2` class, and then built using `MeshBuilder`.  A working example is given below (assuming the executable will be in `build/Schemes` for example).

\code{.cpp}
#include "mesh.hpp"
#include "import_mesh.hpp"
#include "mesh_builder.hpp"

using namespace HArDCore2D;

int main() {

  // Build mesh
  MeshBuilder builder = MeshBuilder(mesh_file);
  std::unique_ptr<Mesh> mesh_ptr = builder.build_the_mesh();

  // Get the BC and re-order the edges (useful to set BC for hybrid schemes)
  std::string bc_id = vm["bc_id"].as<std::string>();
  BoundaryConditions BC(bc_id, *mesh_ptr.get());
  BC.reorder_edges();

	// Create an HybridCore instance
  HybridCore hho(mesh_ptr.get(), K+1, K, use_threads, output);
}
\endcode

<i>Note</i>: the `typ2` format allows for meshes with very generic polygonal cells, including non-convex cells.
However, the builder assumes that each cell is star-shaped with respect to the isobarycenter of its vertices -- otherwise, the calculation of the center of mass may be incorrect. Similarly, the quadrature rules assume that each cell is star-shaped with respect to its center of mass.



<a name="schemes">
\section schemes Schemes
</a>

The following schemes are currently available in HArD::Core2D. The Hybrid High-Order schemes follow the implementation principles described in Appendix B of the book available at https://hal.archives-ouvertes.fr/hal-02151813.

 - [HHO_diffusion](@ref HHO_Diffusion): Hybrid High-Order (HHO) for \f$-\mathrm{div}(K\nabla u)=f\f$, for Dirichlet, Neumann or mixed boundary conditions, with \f$K\f$ a diffusion tensor that is piecewise constant on the mesh.

 - [HHO_locvardiff](@ref HHO_LocVarDiff): HHO for \f$-\mathrm{div}(K\nabla u)=f\f$, for Dirichlet, Neumann or mixed boundary conditions, with \f$K\f$ a diffusion tensor that can vary in each cell.

 - [HHO_diffadvecreac](@ref HHO_DiffAdvecReac): Hybrid High-Order (HHO) for \f$-\mathrm{div}(K\nabla u+\beta u)+\mu u=f\f$, for Dirichlet or mixed boundary conditions, with \f$K\f$ a diffusion tensor that is piecewise constant on the mesh.

 - [LEPNC_diffusion](@ref HArDCore2D::LEPNC_diffusion), in module [LEPNC](@ref LEPNC): Locally Enriched Polytopal Non-Conforming (LEPNC) method for the pure diffusion problem \f$-\mathrm{div}(K\nabla\zeta(u))=f\f$.

 - [LEPNC_StefanPME](@ref HArDCore2D::LEPNC_StefanPME), in module [LEPNC](@ref LEPNC): Locally Enriched Polytopal Non-Conforming (LEPNC) method for the stationnary Stefan/PME problem \f$u-\mathrm{div}(K\nabla\zeta(u))=f\f$.

 - [LEPNC_StefanPME_Transient](@ref HArDCore2D::LEPNC_StefanPME_Transient), in module [LEPNC](@ref LEPNC): LEPNC for the transient Stefan/PME problem \f$\partial_t u-\mathrm{div}(K\nabla\zeta(u))=f\f$.

 - [HMM_StefanPME_Transient](@ref HArDCore2D::HMM_StefanPME_Transient), in module [HMM](@ref HMM): Hybrid Mimetic Mixed (HMM) method for the transient Stefan/PME problem \f$\partial_t u-\mathrm{div}(K\nabla\zeta(u))=f\f$.

 - [DDR_rmplate](@ref DDR_rmplate): Discrete de Rham (DDR) scheme for the Reissner-Mindlin plate bending problem.

The directory `runs` contains BASH to run series of tests on families of meshes. The files `data.sh` describe the parameters of the test cases (polynomial degrees, boundary conditions, mesh families, etc.). The script produces results in the `output` directory, including a pdf file `rate.pdf` describing the rates of convergence in various energy norms.

To run the scripts as they are, you will need `pdflatex` and (for the LEPNC and HMM schemes) a FORTRAN compiler (adjust the `Makefile` to your compiler).




