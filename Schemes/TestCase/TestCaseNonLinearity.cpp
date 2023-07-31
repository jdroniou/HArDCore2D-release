// Class to provide a non-linear function
//
//
// Author: Jerome Droniou (jerome.droniou@monash.edu)
//

#include "TestCaseNonLinearity.hpp"
#include <string>
#include <iostream>
#include <cmath>

using namespace HArDCore2D;

// ----------------------------------------------------------------------------
//                          Implementation
// ----------------------------------------------------------------------------

// Class
TestCaseNonLinearity::TestCaseNonLinearity(const int iTCNL, const double m)
	: m_iTCNL(iTCNL),
    m_m(m) { }

/////////////////// NONLINEARITY ///////////////////////////////:

double TestCaseNonLinearity::nonlinearity(double s, std::string type) {
		double val = 0;

		switch(m_iTCNL){
		/// iTCR[0]=-1: \f$f(u)=0\f$
		case -1: 
				break;
		/// iTCR[0]=1: \f$f(u)=u\f$
		case 1: 
				if (type == "fct"){
					val = s;
				}else if (type == "nlin" || type == "der"){
					val = 1;
				}else if (type == "hess"){
					val = 0;
				}else{
					std::cout << "Type nonlinear function unknown: " << type << "\n";
					exit(EXIT_FAILURE);
				}
				break;
		/// iTCR[0]=2: \f$f(u)=|u|^{m-1}u = |u|^m sign(u)\f$
		case 2: 
				if (type == "fct"){
					val = std::pow(std::abs(s), m_m) * sign(s);
				}else if (type == "nlin"){
					val = std::pow(std::abs(s), m_m-1);
				}else if (type == "der"){
					val = m_m * std::pow(std::abs(s), m_m-1);
				}else if (type == "hess"){
          if (s != 0){
            val = m_m * (m_m-1) * std::pow(std::abs(s), m_m-2) * sign(s);
          }
				}else{
					std::cout << "Type nonlinear function unknown: " << type << "\n";
					exit(EXIT_FAILURE);
				}
				break;
		/// iTCR[0]=3: \f$f(u)=max(u,0)^2\f$
		case 3: 
				if (type == "fct"){
					val = std::pow(std::max(s,double(0)), 2);
				}else if (type == "nlin"){
					val = std::max(s,double(0));
				}else if (type == "der"){
					val = 2*std::max(s,double(0));
				}else if (type == "hess"){
					val = 0;
				}else{
					std::cout << "Type nonlinear function unknown: " << type << "\n";
					exit(EXIT_FAILURE);
				}
				break;
		/// iTCR[0]=4: \f$f(u)=(u+th)^+-th + eps u\f$
		case 4: 
				{
					double eps = 0;
					double th = 0;
					if (type == "fct"){
						val = std::max(s+th,double(0))-th + eps*s;
					}else if (type == "nlin" || type == "der"){
						if (s+th<0){
							val = 0 + eps;
						}else{
							val = 1 + eps;
						}
					}else if (type == "hess"){
						val = 0;
					}else{
						std::cout << "Type nonlinear function unknown: " << type << "\n";
						exit(EXIT_FAILURE);
					}
					break;
				}
		/// iTCR[0]=5: \f$f(u)=u\f$ truncated at 0 and th
		case 5: 
				{
					double th = 2;
					if (type == "fct"){
						val = std::min(th, std::max(s,double(0)));
					}else if (type == "nlin" || type == "der"){
						if (0<=s && s<=th){
							val = 1;
						}else{
							val = 0;
						}
					}else if (type == "hess"){
						val = 0;
					}else{
						std::cout << "Type nonlinear function unknown: " << type << "\n";
						exit(EXIT_FAILURE);
					}
					break;
				}
    /// iTCR[0]=6: Classical Stefan \f$f(s)=s\f$ if \f$s<0\f$, \f$f(s)=0\f$ if \f$0\le s\le 1\f$, \f$f(s)=s-1\f$ if \f$s>1\f$.
		case 6: 
        {
          if (type == "fct"){
            val = std::min(s,double(0)) + std::max(s-1.0,double(0));
          }else if (type == "nlin" || type == "der"){
            val = ( (s<0) ? 1 : 0 ) + ( (s>1) ? 1 : 0 );
          }else if (type == "hess"){
						val = 0;
					}else{
						std::cout << "Type nonlinear function unknown: " << type << "\n";
						exit(EXIT_FAILURE);
					}
  				break;
        }
    // Default unknown
		default: 
		  std::cout << "[TestCaseNonLinear] Test case " + std::to_string(m_iTCNL) + " unknown" << std::endl;
      exit(1);
    break;
	}
	
	return val;

}

