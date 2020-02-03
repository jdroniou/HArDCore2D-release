// Class to Provides various test cases (diffusion, exact solution, and their derivatives)
//
//
// Author: Jerome Droniou (jerome.droniou@monash.edu)
//

#include "TestCase.hpp"
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
TestCase::TestCase(const std::vector<int> iTC)
  : m_iTC(iTC),
    _deg_diff(0),
    _lambda(pow(10,6)) {
    validate();
    if (m_iTC[1]==2 || m_iTC[1]==4){
      _deg_diff = 2;
    }
  }

/////////////////// SOLUTION ///////////////////////////////:

// Solution
double TestCase::sol(const double x, const double y){
  double u = 0;
  switch(m_iTC[0]){
      /// iTC[0]=1: \f$u(x,y)=sin(\pi x)  sin(\pi y)\f$
    case 1: u = sin(pi*x) * sin(pi*y);
            break;
      /// iTC[0]=2: \f$u(x,y)=cos(\pi x)  cos(\pi y)\f$
    case 2: u = cos(pi*x) * cos(pi*y);    
            break;
      /// iTC[0]=3: \f$u(x,y)= x\f$
    case 3: u = x;    
            break;
      /// iTC[0]=4: \f$u(x,y)= y\f$
    case 4: u = y;    
            break;
      /// iTC[0]=5: \f$u(x,y)= x^2 + y^2\f$
    case 5: u = pow(x, 2) + pow(y, 2); 
            break;
      /// iTC[0]=6: \f$u(x,y)= e^x \sin(pi x) \cos(\pi x)\f$
    case 6: u = exp(x) * sin(pi*x) * sin(pi*y); 
            break;
    default: break;
  }
  return u;
}

// Gradient of the solution
Eigen::Vector2d TestCase::grad_sol(const double x, const double y, const Cell* cell){
  Eigen::Vector2d G = Eigen::Vector2d::Zero();
  switch(m_iTC[0]){
    case 1: G(0) = pi * cos(pi*x) * sin(pi*y);
            G(1) = pi * sin(pi*x) * cos(pi*y);
            break;

    case 2: G(0) = -pi * sin(pi*x) * cos(pi*y);
            G(1) = -pi * cos(pi*x) * sin(pi*y);
            break;

    case 3: G(0) = 1;
            G(1) = 0;
            break;

    case 4: G(0) = 0;
            G(1) = 1;
            break;

    case 5: G(0) = 2*x;
            G(1) = 2*y;
            break;

    case 6: G(0) = exp(x) * sin(pi*x) * sin(pi*y) + exp(x) * pi * cos(pi*x) * sin(pi*y) ;
            G(1) = exp(x) * pi * sin(pi*x) * cos(pi*y);
            break;

    default: break;
  }
  return G;
}

// Hessian of the solution
Eigen::Matrix2d TestCase::hess_sol(const double x, const double y, const Cell* cell){
  Eigen::MatrixXd H = Eigen::Matrix2d::Zero();
  switch(m_iTC[0]){

    case 1: H.row(0) << - pi*pi*sin(pi*x)*sin(pi*y), pi*pi*cos(pi*x)*cos(pi*y);
            H.row(1) <<  pi*pi*cos(pi*x)*cos(pi*y), -pi*pi*sin(pi*x)*sin(pi*y);
            break;

    case 2: H.row(0) << - pi*pi*cos(pi*x)*cos(pi*y), pi*pi*sin(pi*x)*sin(pi*y);
            H.row(1) <<  pi*pi*sin(pi*x)*sin(pi*y), -pi*pi*cos(pi*x)*cos(pi*y);
            break;

    case 3: break;
    case 4: break;
    case 5: H.row(0) << 2, 0;
            H.row(1) << 0, 2;
            break;
    case 6: H.row(0) << exp(x)*sin(pi*x)*sin(pi*y) + 2*pi*exp(x)*cos(pi*x)*sin(pi*y) - pi*pi*exp(x)*sin(pi*x)*sin(pi*y), pi*exp(x)*sin(pi*x)*cos(pi*y)+pi*pi*exp(x)*cos(pi*x)*cos(pi*y);
            H.row(1) << pi*exp(x)*sin(pi*x)*cos(pi*y)+pi*pi*exp(x)*cos(pi*x)*cos(pi*y), -pi*pi*exp(x)*sin(pi*x)*sin(pi*y);
            break;

    default: break;
  }
  return H;
}

//////////////////////////// DIFFUSION /////////////////////////////

// Diffusion matrix
Eigen::Matrix2d TestCase::diff(const double x, const double y, const Cell* cell){
  Eigen::Matrix2d K = Eigen::Matrix2d::Identity();
  switch(m_iTC[1]){
      /// iTC[1]=1: Diff = Id
    case 1: break;    
      /// iTC[1]=2: Diff = \f$\left[\begin{array}{cc}y^2+1 & -xy\\  -xy & x^2+1\end{array}\right]\f$
    case 2: K.row(0) << pow(y,2)+1, -x*y;        
            K.row(1) << -x*y , pow(x,2)+1;
            break;
      /// iTC[1]=3: Diff=\f$\left[\begin{array}{cc}\lambda & 0\\ 0 & 1\end{array}\right]\f$ if \f$y<1/2\f$, Diff=Id if \f$y\ge 1/2\f$. Only valid with iTC[0]=2 (\f$\partial_x u\f$ must vanish along \f$y=1/2\f$). 
    case 3: if (cell->center_mass().y()<0.5){      
              K.row(0) << _lambda, 0;
            }else{
              K.row(0) << 1, 0;
            }
            K.row(1) << 0, 1;
            break;
      /// iTC[1]=4: rotating diffusion. Diff=\f$\left[\begin{array}{cc}\epsilon \bar{x}^2 +  \bar{y}^2 & (\epsilon-1)\bar{x}\bar{y}\\ (\epsilon-1)\bar{x}\bar{y} & \bar{x}^2+\epsilon \bar{y}^2\end{array}\right]\f$, where \f$\bar{x}=x+0.1\f$ and \f$\bar{y}=y+0.1\f$.
    case 4: {
            double barx = x+0.1;
            double bary = y+0.1;
            K.row(0) << eps*barx*barx + bary*bary, (eps-1)*barx*bary;
            K.row(1) << (eps-1)*barx*bary, barx*barx + eps*bary*bary;
            break;
            }
      /// iTC[1]=5: Diff=\f$\left[\begin{array}{cc}\lambda & 0\\ 0 & 1\end{array}\right]\f$. 
    case 5: K.row(0) << _lambda, 0;
            K.row(1) << 0, 1;
            break;
    default: break;
  }
  return K;
}

// Divergence by row of the diffusion matrix
Eigen::Vector2d TestCase::div_diff(const double x, const double y, const Cell* cell){
  Eigen::Vector2d divK = Eigen::Vector2d::Zero();
  switch(m_iTC[1]){
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
double TestCase::source(const double x, const double y, const Cell* cell){

  Eigen::Matrix2d AHu = diff(x,y,cell) * hess_sol(x,y,cell);

  return -AHu.trace() - div_diff(x,y,cell).dot(grad_sol(x,y,cell));
}



///////////////////////////// VALIDATION ////////////////////////////

void TestCase::validate(){
  
  if (m_iTC[0]>6 || m_iTC[1]>5 || (m_iTC[1]==3 && m_iTC[0] !=2)){
    std::cout << "Incorrect choice of test cases: iTC= " << m_iTC[0] << ", " << m_iTC[1] << "\n";
    exit(EXIT_FAILURE);
  }  

}

