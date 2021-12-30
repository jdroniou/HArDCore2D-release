#ifndef PLATESCORE_HPP
#define PLATESCORE_HPP

#include <memory>
#include <iostream>

#include <basis.hpp>
#include <polynomialspacedimension.hpp>

/*!	
 * @defgroup PlatesCore 
 * @brief Classes providing tools to the DDR sequence for plates
 */
 
namespace HArDCore2D
{

  /*!
   *	\addtogroup PlatesCore
   * @{
   */

  /// Construct all polynomial spaces for the plates sequence
  class PlatesCore
  {
  public:
    // Types for element bases
    typedef Family<MonomialScalarBasisCell> PolyBasisCellType;
    typedef TensorizedVectorFamily<PolyBasisCellType, dimspace> Poly2BasisCellType;
    typedef Family<MatrixFamily<PolyBasisCellType, dimspace>> PolySymBasisCellType;
    typedef Family<HessianBasis<ShiftedBasis<MonomialScalarBasisCell> > > HolyBasisCellType;
    typedef Family<HolyComplBasisCell> HolyComplBasisCellType;
    
    // Types for edge bases
    typedef Family<MonomialScalarBasisEdge> PolyBasisEdgeType;

    /// Structure to store element bases
    struct CellBases
    {
      /// Geometric support
      typedef Cell GeometricSupport;

      std::unique_ptr<PolyBasisCellType> Polykp1;
      std::unique_ptr<PolyBasisCellType> Polykm2;
      std::unique_ptr<Poly2BasisCellType> Poly2km1;
      std::unique_ptr<PolySymBasisCellType> PolySymkm1;
      std::unique_ptr<HolyBasisCellType> Holykm4;
      std::unique_ptr<HolyComplBasisCellType> HolyComplkm1;
    };

    /// Structure to store edge bases
    struct EdgeBases
    {
      /// Geometric support
      typedef Edge GeometricSupport;

      std::unique_ptr<PolyBasisEdgeType> Polykm1;
      std::unique_ptr<PolyBasisEdgeType> Polykm2;
      std::unique_ptr<PolyBasisEdgeType> Polykm3;
    };

    /// Constructor
    PlatesCore(const Mesh & mesh, size_t K, bool use_threads = true, std::ostream & output = std::cout);

    /// Return a const reference to the mesh
    const Mesh & mesh() const
    {
      return m_mesh;
    }

    /// Return the polynomial degree
    const size_t & degree() const
    {
      return m_K;
    }
    
    /// Return cell bases for element iT
    inline const CellBases & cellBases(size_t iT) const
    {
      // Make sure that the basis has been created
      assert( m_cell_bases[iT] );
      return *m_cell_bases[iT].get();
    }

    /// Return edge bases for edge iE
    inline const EdgeBases & edgeBases(size_t iE) const
    {
      // Make sure that the basis has been created
      assert( m_edge_bases[iE] );
      return *m_edge_bases[iE].get();
    }
  private:
    /// Compute the bases on an element T
    CellBases _construct_cell_bases(size_t iT);

    /// Compute the bases on an edge E
    EdgeBases _construct_edge_bases(size_t iE);
    
    // Pointer to the mesh
    const Mesh & m_mesh;
    // Degree
    const size_t m_K;
    // Output stream
    std::ostream & m_output;
    
    // Cell bases
    std::vector<std::unique_ptr<CellBases> > m_cell_bases;
    // Edge bases
    std::vector<std::unique_ptr<EdgeBases> > m_edge_bases;
  };
} // namespace HArDCore2D

#endif
