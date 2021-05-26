// Free functions to help handle BC and DOFs/indices
//
// Author: Jerome Droniou (jerome.droniou@monash.edu)
//

/*
*
*	This library was developed around HHO methods, although some parts of it have a more
* general purpose. If you use this code or part of it in a scientific publication, 
* please mention the following book as a reference for the underlying principles
* of HHO schemes:
*
* The Hybrid High-Order Method for Polytopal Meshes: Design, Analysis, and Applications.
* D. A. Di Pietro and J. Droniou. 2019, 516p. 
* url: https://hal.archives-ouvertes.fr/hal-02151813.
*
*/

#include <BoundaryConditions/BoundaryConditions.hpp>
#include <ddrspace.hpp>
#include <mesh.hpp>

using namespace HArDCore2D;

// ----------------------------------------------------------------------------
//                            Class definition
// ----------------------------------------------------------------------------
 
  /*!
   *	\addtogroup BoundaryConditions
   * @{
   */

/// Adds BC labels do ddrspace DOFs. The default label is 0; we leave it 0 for internal DOF, -1 for Neumann DOF and 1 for Dirichlet DOF
void setBCLabels(const BoundaryConditions &BC, DDRSpace &ddrspace) {
  
  for (Vertex *v : ddrspace.mesh().get_vertices()){
    int label = 0;
    if (BC.type(*v)=="dir"){
      label = 1;
    }else if (BC.type(*v)=="neu"){
      label = -1;
    }
    size_t offset_v = ddrspace.globalOffset(*v);
    for (size_t i=offset_v; i < offset_v + ddrspace.numLocalDofsVertex(); i++){
      ddrspace.setLabelDOF(i, label);
    }          
  }

  for (Edge *e : ddrspace.mesh().get_edges()){
    int label = 0;
    if (BC.type(*e)=="dir"){
      label = 1;
    }else if (BC.type(*e)=="neu"){
      label = -1;
    }
    size_t offset_e = ddrspace.globalOffset(*e);
    for (size_t i=offset_e; i < offset_e + ddrspace.numLocalDofsEdge(); i++){
      ddrspace.setLabelDOF(i, label);
    }          
  }

}


/// Function to offset and index i according to a vector c0,c1,...,c2n of increasing numbers
/**  The indices between [c0,c1),...[c2n-1,c2n) are removed (return -1), and the other are offset according to this removal
*/
int offsetIndex(const std::vector<size_t> &c, const int &i){
  assert( c.size() >= 2 && c.size()%2 == 0);
  
  int val = -1;
  int idx = -1;
  int offset = 0;
  // Look for index in c such that c[idx]>i and add offsets (every other pair of indices)
  do {
    idx++;
    if (idx%2 == 1){
      // c.size() even so no risk of having idx=c.size() here
      offset += c[idx]-c[idx-1];
    } 
  }
  while (idx<int(c.size()) && int(c[idx])<=i);
  
  if (idx%2 == 0){
    val = i - offset;
  }
  
  return val;
  
}

/// Create a map from DOFs 0..N-1 to values obtained by cutting the DOFs corresponding to c (as per offsetIndex).
Eigen::ArrayXi create_mapDOF(const std::vector<size_t> &c, const size_t N){
  assert( c.size() >= 2 && c.size()%2 == 0);
  
  Eigen::ArrayXi dofs = Eigen::ArrayXi::LinSpaced(N, 0, N-1);
  Eigen::ArrayXi map = -Eigen::ArrayXi::Ones(N);
  
  size_t offset = 0;
  map.head(c[0]) = dofs.head(c[0]);
  for (size_t i = 1; i < c.size()-1; i += 2){
    offset += c[i] - c[i-1];
    map.segment(c[i], c[i+1]-c[i]) = dofs.segment(c[i], c[i+1]-c[i]) - offset;
  }
  offset += c[c.size()-1]-c[c.size()-2];
  map.tail(N-c[c.size()-1]) = dofs.tail(N-c[c.size()-1]) - offset;
  
  return map;

}


/// Replace sections of vector V by values from vector Z into vector V; the sections are determined by 'sec'. The vectors are any (identical) types of Eigen::Vector
/** Returns a vector of same length as V, whose segments starting at sec[i].first and with length sec[i].second
  are replaced using the values (successively) of Z.
  Can be used to:
  
   1) create a map from DOFs to unknowns: starting from a vector V of -1 of size nb of DOFs and a linearly spaced vector Z of size
    the number of final unknowns, replaces the -1 by the unknowns. The final vector maps -1 for DOFs that are eliminated (not in
    the 'sec' positions) and linearly lists unknown numbers for the other DOFs

   2) insert, in a vector already containing values for Dirichlet DOFs, values calculated by solving a system on the other DOFs.
*/
template<typename VecType>
VecType replaceSectionsVector(const VecType &V, const VecType &Z, const std::vector<std::pair<size_t,size_t>> &sec){
  assert ( sec[sec.size()-1].first + sec[sec.size()-1].second <= size_t(V.rows()) );
  
  VecType val = V;
  size_t posZ = 0;
  for (size_t i=0; i < sec.size(); i++){
    assert( posZ + sec[i].second <= size_t(Z.rows()) );
    
    val.segment(sec[i].first, sec[i].second) = Z.segment(posZ, sec[i].second);
    posZ += sec[i].second;
  }
  return val;
}

//@}

