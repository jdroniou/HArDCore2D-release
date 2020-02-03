// Class to Provides various test cases (diffusion, exact solution, and their derivatives)
//
//
// Author: Jerome Droniou (jerome.droniou@monash.edu)
//

#include "TestCaseStefanPME.hpp"
#include "cell.hpp"
#include <memory>
#include <string>
#include <vector>
#include <iostream>

#include <Eigen/Dense>

using namespace HArDCore2D;

// ----------------------------------------------------------------------------
//                          Implementation
// ----------------------------------------------------------------------------

// Class - we set up the anisotropy when we create the class
TestCaseStefanPME::TestCaseStefanPME(const std::vector<int> iTCS)
  : m_iTCS(iTCS),
    _deg_diff(0),
    _lambda(pow(10,6)),
    m_tcase({std::max(m_iTCS[0],1), m_iTCS[1]}) {
      validate();
      if (m_iTCS[1]==2 || m_iTCS[1]==4){
        _deg_diff = 2;
      }
//      if (m_iTCS[0]>0){
//        TestCase m_tcase(m_iTCS);
//      }
  }

/////////////////// SOLUTION ///////////////////////////////:

// Solution
double TestCaseStefanPME::sol(const double x, const double y){
  double u = 0;
  if (m_iTCS[0]>0){
    u = m_tcase.sol(x, y);
  }else{
    switch(m_iTCS[0]){
      /// iTCS[0]=-1, \f$u=\cosh(\frac{x+y}{\sqrt{2}}-\gamma\f$ if \f$\frac{x+y}{\sqrt{2}}>\gamma\f$.
      case -1: {
               double s=(x+y)/sqrt2;
               if (s>gamma){             
                 u = cosh(s - gamma);
               }
               break;
              }
      /// iTCS[0]=-2, \f$u=(\frac{x+y}{\sqrt{2}}-\frac12)^3\f$.
      case -2: {
               double s=(x+y)/sqrt2;
               u = std::pow(s - 0.5, 3);
               break;
              }
      /// iTCS[0]=-3, \f$u=max(\rho - r^2, 0)\f$ where \f$r^2=(x-.5)^2+(y-.5)^2\f$.
      case -3: u = std::max(rho-std::pow(radial(x,y), 2), 0.0);
               break;
      default: break;
    }
  }
  return u;
}

// Gradient of the solution
Eigen::Vector2d TestCaseStefanPME::grad_sol(const double x, const double y, const Cell* cell){
  Eigen::Vector2d G = Eigen::Vector2d::Zero();
  if (m_iTCS[0]>0){
    G = m_tcase.grad_sol(x, y, cell);
  }else{
    switch(m_iTCS[0]){
      case -1: {
                double s=(x+y)/sqrt2;
                if (s>gamma){             
                  G(0) = 1.0/sqrt2 * sinh(s - gamma);
                  G(1) = 1.0/sqrt2 * sinh(s - gamma);
                }
                break;
              }
      case -2: {
                double s=(x+y)/sqrt2;
                G(0) = 3.0 * std::pow(s-0.5, 2)/sqrt2;
                G(1) = 3.0 * std::pow(s-0.5, 2)/sqrt2;
                break;
              }
      case -3: {
                bool inside = (rho>std::pow(radial(x,y), 2));
                G(0) = (inside ? -2*translate(x,y).x() : 0);
                G(1) = (inside ? -2*translate(x,y).y() : 0);
                break;
              }
      default: break;
    }
  }
  return G;
}

// Hessian of the solution
Eigen::Matrix2d TestCaseStefanPME::hess_sol(const double x, const double y, const Cell* cell){
  Eigen::MatrixXd H = Eigen::Matrix2d::Zero();
  if (m_iTCS[0]>0){
    H = m_tcase.hess_sol(x, y, cell);
  }else{
    switch(m_iTCS[0]){
      case -1: {
                double s=(x+y)/sqrt2;
                if (s>gamma){     
                  H.row(0) << 0.5*cosh(s-gamma), 0.5*cosh(s-gamma);
                  H.row(1) << 0.5*cosh(s-gamma), 0.5*cosh(s-gamma);
                }
                break;
              }
      case -2: {
                double s=(x+y)/sqrt2;
                H.row(0) << 3.0 * (s-0.5), 3.0*(s-0.5);
                H.row(1) << 3.0 * (s-0.5), 3.0*(s-0.5);
                break;
              }
      case -3: {
                bool inside = (rho>std::pow(radial(x,y), 2));
                if (inside){
                  H = -2 * Eigen::Matrix2d::Identity();
                }
                break;
              }
      default: break;
    }
  }
  return H;
}

//////////////////////////// DIFFUSION /////////////////////////////

// Diffusion matrix
Eigen::Matrix2d TestCaseStefanPME::diff(const double x, const double y, const Cell* cell){
  Eigen::Matrix2d K = Eigen::Matrix2d::Identity();
  switch(m_iTCS[1]){
      /// iTCS[1]=1: Diff = Id
    case 1: break;    
      /// iTCS[1]=2: Diff = \f$\left[\begin{array}{cc}y^2+1 & -xy\\  -xy & x^2+1\end{array}\right]\f$
    case 2: K.row(0) << pow(y,2)+1, -x*y;        
            K.row(1) << -x*y , pow(x,2)+1;
            break;
      /// iTCS[1]=3: Diff=\f$\left[\begin{array}{cc}\lambda & 0\\ 0 & 1\end{array}\right]\f$ if \f$y<1/2\f$, Diff=Id if \f$y\ge 1/2\f$. Only valid with iTCS[0]=2 (\f$\partial_x u\f$ must vanish along \f$y=1/2\f$). 
    case 3: if (cell->center_mass().y()<0.5){      
              K.row(0) << _lambda, 0;
            }else{
              K.row(0) << 1, 0;
            }
            K.row(1) << 0, 1;
            break;
      /// iTCS[1]=4: rotating diffusion. Diff=\f$\left[\begin{array}{cc}\epsilon \bar{x}^2 +  \bar{y}^2 & (\epsilon-1)\bar{x}\bar{y}\\ (\epsilon-1)\bar{x}\bar{y} & \bar{x}^2+\epsilon \bar{y}^2\end{array}\right]\f$, where \f$\bar{x}=x+0.1\f$ and \f$\bar{y}=y+0.1\f$.
    case 4: {
            double barx = x+0.1;
            double bary = y+0.1;
            K.row(0) << eps*barx*barx + bary*bary, (eps-1)*barx*bary;
            K.row(1) << (eps-1)*barx*bary, barx*barx + eps*bary*bary;
            break;
            }
      /// iTCS[1]=5: Diff=\f$\left[\begin{array}{cc}\lambda & 0\\ 0 & 1\end{array}\right]\f$. 
    case 5: K.row(0) << _lambda, 0;
            K.row(1) << 0, 1;
            break;
    default: break;
  }
  return K;
}

// Divergence by row of the diffusion matrix
Eigen::Vector2d TestCaseStefanPME::div_diff(const double x, const double y, const Cell* cell){
  Eigen::Vector2d divK = Eigen::Vector2d::Zero();
  switch(m_iTCS[1]){
    case 1: break;
    case 2: divK(0) = -x;
            divK(1) = -y;
            break;
    case 3: break;
    case 4: {
            double barx = x+0.1;
            double bary = y+0.1;
            divK(0) = 2*eps*barx + (eps-1)*barx;
            divK(1) = (eps-1)*bary + 2*eps*bary;
            break;
            }
    case 5: break;
    default: break;
  }
  return divK;
}

///////////////////////////// SOURCE TERM ///////////////////////////

// Source term
double TestCaseStefanPME::source(const double x, const double y, const Cell* cell){

  Eigen::Matrix2d AHu = diff(x,y,cell) * hess_sol(x,y,cell);

  return -AHu.trace() - div_diff(x,y,cell).dot(grad_sol(x,y,cell));
}



///////////////////////////// VALIDATION ////////////////////////////

void TestCaseStefanPME::validate(){
  
  if (m_iTCS[0]<-3 || m_iTCS[0]==0 || m_iTCS[1]>5 || (m_iTCS[1]==3 && m_iTCS[0] !=2)){
    std::cout << "Incorrect choice of test cases: iTCS= " << m_iTCS[0] << ", " << m_iTCS[1] << "\n";
    exit(EXIT_FAILURE);
  }  

}

