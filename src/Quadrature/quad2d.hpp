// Creates quadrature rule in a cell
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


#ifndef QUAD2D_HPP
#define QUAD2D_HPP
#include <vector>
#include <triangle_dunavant_rule.hpp>
#include <assert.h>
#include <string>


/*!  
*  @defgroup Quadratures 
* @brief Classes providing quadratures on edges and in cells
*/

namespace HArDCore2D {

// \addtogroup Quadratures
//@{

/**
* @brief Wrapper for dunavant quadrature rules
*/
class QuadRuleTriangle {
    static constexpr size_t max_doe = 20;

public:
    /**
    * @brief Default constructor
    *
    * @param doe degrees of exactness (e.g. how many points for approximating
    *integral
    * @param warn
    */
    QuadRuleTriangle(size_t doe, bool warn);
    ~QuadRuleTriangle();

    size_t nq();             /// <\brief returns number of points
    double xq(size_t i);  /// <\brief
    double yq(size_t i);  /// <\brief
    double wq(size_t i);  /// <\brief get the weight for a given point
    void setup(double xV[], double yV[]);  /// <\brief setup the quad rule given
                                           /// vertex coords

private:
    size_t _npts;
    double* _xy;
    double* _w;
    double* _xyphys;
    double area;
};

//@}
}
#endif /* QUAD2D_HPP */
