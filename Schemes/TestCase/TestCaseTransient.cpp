// Class to Provides various test cases (diffusion, exact solution, and their derivatives)
//
//
// Author: Jerome Droniou (jerome.droniou@monash.edu)
//

#include "TestCaseTransient.hpp"
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
TestCaseTransient::TestCaseTransient(const std::vector<int> iTC, const double mPME)
  : m_iTC(iTC),
    _deg_diff(0),
    _lambda(pow(10,6)),
    m_mPME(mPME),
    m_alpha(1/mPME),
    m_rho(m_alpha/2.0),
    m_gamma(m_alpha*(mPME-1.0)/(4.0*mPME)) {
    validate();
    if (m_iTC[1]==2 || m_iTC[1]==4){
      _deg_diff = 2;
    }
  }

/////////////////// SOLUTION ///////////////////////////////:

// Solution
std::function<double(const double&, const VectorRd&)> TestCaseTransient::solution(){
  std::function<double(const double&, const VectorRd&)> u = [](double t, VectorRd p)->double {return 0.0;};

  switch(m_iTC[0]){
      /// iTC[0]=1: \f$u(x,y)=cos(t) sin(\pi x)  sin(\pi y)\f$
    case 1: 
        u = [this](double t, VectorRd p) -> double {
            return cos(t) * sin(pi*p.x()) * sin(pi*p.y());
          };
        break;

      /// iTC[0]=2: \f$u(x,y)=cos(t) cos(\pi x)  cos(\pi y)\f$
    case 2: 
        u = [this](double t, VectorRd p) -> double {
            return cos(t) * cos(pi*p.x()) * cos(pi*p.y());    
          };
        break;

      /// iTC[0]=3: \f$u(x,y)= (1-t) x\f$
    case 3: 
        u = [](double t, VectorRd p) -> double {
            return (2.0-t) * p.x();
          };
        break;

      /// iTC[0]=4: Barenblatt solution (with offset)
    case 4:
        u = [this](double t, VectorRd p) -> double {
            double tps = t+m_t0;
            double normz2 = (p-VectorRd(.5,.5)).squaredNorm();
            double val = m_CB-m_gamma*normz2*std::pow(tps, -2.0*m_rho);
            val = std::pow(std::max(0.0, val), 1.0/(m_mPME-1.0));
            return val * std::pow(tps, -m_alpha);
        };
        break;

    default: break;
  }
  return u;
}

// Gradient of the solution
std::function<VectorRd(const double&, const VectorRd&, const Cell*)> TestCaseTransient::grad_solution(){
  
  std::function<VectorRd(const double&, const VectorRd&, const Cell*)> grad = [](const double t, const VectorRd p, const Cell* cell)->VectorRd {return VectorRd::Zero();};
  switch(m_iTC[0]){
    case 1: 
      grad = [this](const double t, const VectorRd p, const Cell* cell)->VectorRd {
          VectorRd G = VectorRd::Zero();
          G(0) = cos(t) * pi * cos(pi*p.x()) * sin(pi*p.y());
          G(1) = cos(t) * pi * sin(pi*p.x()) * cos(pi*p.y());
          return G;
        };
        break;

    case 2: 
      grad = [this](const double t, const VectorRd p, const Cell* cell)->VectorRd {
          VectorRd G = VectorRd::Zero();
          G(0) = -cos(t) * pi * sin(pi*p.x()) * cos(pi*p.y());
          G(1) = -cos(t) * pi * cos(pi*p.x()) * sin(pi*p.y());
          return G;
        };
        break;

    case 3: 
      grad = [](const double t, const VectorRd p, const Cell* cell)->VectorRd {
          VectorRd G = VectorRd::Zero();
          G(0) = 2.0-t;
          G(1) = 0;
          return G;
        };
        break;

    default: break;
  }
  return grad;
}

// Hessian of the solution
std::function<Eigen::Matrix2d(const double&, const VectorRd&, const Cell*)> TestCaseTransient::hess_solution(){

  std::function<Eigen::Matrix2d(const double&, const VectorRd&, const Cell*)> Hess = [](const double t, const VectorRd p, const Cell* cell)->Eigen::Matrix2d { return Eigen::Matrix2d::Zero();};
  switch(m_iTC[0]){

    case 1: 
      Hess = [this](const double t, const VectorRd p, const Cell* cell)->Eigen::Matrix2d {
            Eigen::Matrix2d H = Eigen::Matrix2d::Zero();
            H.row(0) << - cos(t) * pi*pi*sin(pi*p.x())*sin(pi*p.y()), cos(t) * pi*pi*cos(pi*p.x())*cos(pi*p.y());
            H.row(1) <<  cos(t) * pi*pi*cos(pi*p.x())*cos(pi*p.y()), -cos(t) * pi*pi*sin(pi*p.x())*sin(pi*p.y());
            return H;
          };
          break;

    case 2:       
      Hess = [this](const double t, const VectorRd p, const Cell* cell)->Eigen::Matrix2d {
            Eigen::Matrix2d H = Eigen::Matrix2d::Zero();
            H.row(0) << - cos(t) * pi*pi*cos(pi*p.x())*cos(pi*p.y()), cos(t) * pi*pi*sin(pi*p.x())*sin(pi*p.y());
            H.row(1) <<  cos(t) * pi*pi*sin(pi*p.x())*sin(pi*p.y()), -cos(t) * pi*pi*cos(pi*p.x())*cos(pi*p.y());
            return H;
          };
          break;

    case 3: break;

    default: break;
  }
  return Hess;
}

// Time derivative of the solution
std::function<double(const double &, const VectorRd &)> TestCaseTransient::delt_solution(){

  std::function<double(const double &, const VectorRd &)> time_der = [](double t, VectorRd p)->double { return 0.0;};

  switch(m_iTC[0]){

    case 1:         
      time_der = [this](double t, VectorRd p) -> double {
        return -sin(t) * sin(pi*p.x()) * sin(pi*p.y());
      };
      break;

    case 2: 
      time_der = [this](double t, VectorRd p) -> double {
        return -sin(t) * cos(pi*p.x()) * cos(pi*p.y());   
      };
      break;

    case 3: 
      time_der = [](double t, VectorRd p) -> double {
        return -p.x();
      };
      break;

    default: break;
  }
  return time_der;
}

//////////////////////////// DIFFUSION /////////////////////////////

// Diffusion matrix
Eigen::Matrix2d TestCaseTransient::diff(const double x, const double y, const Cell* cell){
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
Eigen::Vector2d TestCaseTransient::div_diff(const double x, const double y, const Cell* cell){
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
double TestCaseTransient::div_diff_grad(const double t, const double x, const double y, const Cell* cell){

  Eigen::Matrix2d AHu = diff(x,y,cell) * hess_solution()(t,VectorRd(x,y),cell);

  return  -AHu.trace() - div_diff(x,y,cell).dot(grad_solution()(t,VectorRd(x,y),cell));
}



///////////////////////////// VALIDATION ////////////////////////////

void TestCaseTransient::validate(){
  
  if (m_iTC[0]>4 || m_iTC[1]>5 || (m_iTC[1]==3 && m_iTC[0] !=2)){
    std::cout << "Incorrect choice of test cases: iTC= " << m_iTC[0] << ", " << m_iTC[1] << "\n";
    exit(EXIT_FAILURE);
  }  

}

