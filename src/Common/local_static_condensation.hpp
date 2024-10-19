// Structure to store information for and perform local static condensation

#ifndef LOCAL_STATIC_CONDENSATION_HPP
#define LOCAL_STATIC_CONDENSATION_HPP

#include <Eigen/Dense>

namespace HArDCore2D
{

  /*!
   *	\addtogroup Common
   * @{
   */
   
  /// Structure to store information for, and perform, local static condensation
  /** This structure is used to perform, from a matrix-vector representing local contributions from a cell,
  the static condensation and produce what is necessary to fill in the global matrix-vector of the system and
  the recovery operator.
    - Perm: permutation matrix that puts all the statically condensed DOFs at the end
    - globalDOFs_gl: list, in the order (at the start) they appear after the permutation, of the global DOFs for the unknowns that are globally coupled (not statically condensed) in the system
    - globalDOFs_sc: list, in the order (at the end) they appear after the permutation, of the global DOFs for the unknowns that are statically condensed
  */
  struct LocalStaticCondensation
  {
    /// Constructor
    LocalStaticCondensation(
                    const Eigen::MatrixXd & Perm,         ///< Permutation that moves the DOFs to condense at the end
                    const std::vector<size_t> & globalDOFs_gl, ///< List of global DOFs that are globally coupled (not statically condensed)
                    const std::vector<size_t> & globalDOFs_sc  ///< List of global DOFs that are statically condensed
                    )
         : m_Perm(Perm),
           m_globalDOFs_gl(globalDOFs_gl),
           m_globalDOFs_sc(globalDOFs_sc),
           m_dim_gl(globalDOFs_gl.size()),
           m_dim_sc(globalDOFs_sc.size())
           {
            // do nothing
           };

    /// Compute the local static condensation
    /** From a local pair of Matrix-Vector product, constructed with all the DOFs, returns a 4-tuple; the first two elements are the matrix-vector of the locally statically condensed system, the last two elements are the matrix-vector of the local recovery operator */
    std::tuple<Eigen::MatrixXd, Eigen::VectorXd, Eigen::MatrixXd, Eigen::VectorXd>
      compute(const std::pair<Eigen::MatrixXd, Eigen::VectorXd> & lsT)
    {
      // We permute the DOFs that are statically condensed
      Eigen::MatrixXd AT = m_Perm * lsT.first * m_Perm.transpose();
      Eigen::VectorXd bT = m_Perm * lsT.second;
  
      Eigen::MatrixXd AT_gl, AT_sc;
      Eigen::VectorXd bT_gl, bT_sc;
  
      if (m_dim_sc > 0){
        // Extract 4 blocks of AT, bT for static condensation
        Eigen::MatrixXd A11 = AT.topLeftCorner(m_dim_gl, m_dim_gl);  
        Eigen::MatrixXd A12 = AT.topRightCorner(m_dim_gl, m_dim_sc);
        Eigen::MatrixXd A21 = AT.bottomLeftCorner(m_dim_sc, m_dim_gl);
        Eigen::MatrixXd A22 = AT.bottomRightCorner(m_dim_sc, m_dim_sc);
        Eigen::VectorXd b1 = bT.head(m_dim_gl);
        Eigen::VectorXd b2 = bT.tail(m_dim_sc);
        // Create condensed system (AT_reduced, bT_reduced) and SC recovery operator (RT, cT)
        // Eigen::PartialPivLU<Eigen::MatrixXd> A22pivlu(A22);
        Eigen::FullPivLU<Eigen::MatrixXd> A22pivlu(A22);
        if( !A22pivlu.isInvertible() ) {
          std::cout << "[LocalStaticCondensation] Found non invertible local matrix" << std::endl;
          exit(1);
        } // if
        Eigen::MatrixXd A22inv_A21 = A22pivlu.solve(A21);
        Eigen::VectorXd A22inv_b2 = A22pivlu.solve(b2);
        AT_gl = A11 - A12 * A22inv_A21;
        bT_gl = b1 - A12 * A22inv_b2;
        AT_sc = -A22inv_A21;
        bT_sc = A22inv_b2;
      }else{
        // No static condensation
        AT_gl = AT;
        bT_gl = bT;
        AT_sc = Eigen::MatrixXd::Zero(0,0);
        bT_sc = Eigen::VectorXd::Zero(0);
      }
      
      return std::make_tuple(AT_gl, bT_gl, AT_sc, bT_sc);
    };

    /// Returns global DOFs that are not statically condensend
    inline std::vector<size_t> globalDOFs_gl() { return m_globalDOFs_gl;};
    /// Returns the number of DOFs that are not statically condensed
    inline size_t dim_gl() { return m_dim_gl;};
    
    ///// Functions only for backward compatibility, will be removed at some point.
    inline std::vector<size_t> globalDOFs_sys() { return m_globalDOFs_gl;};
    inline size_t dim_sys() { return m_dim_gl;};
    
    
    /// Returns global DOFs that are statically condensend
    inline std::vector<size_t> globalDOFs_sc() { return m_globalDOFs_sc;};
    /// Returns the number of DOFs that are statically condensed
    inline size_t dim_sc() { return m_dim_sc;};

    // Members
    Eigen::MatrixXd m_Perm;     
    std::vector<size_t> m_globalDOFs_gl;
    std::vector<size_t> m_globalDOFs_sc;
    size_t m_dim_gl;
    size_t m_dim_sc;
  };


  //@}

} // end of namespace HArDCore2D

#endif
