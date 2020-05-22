// Class to provide various test cases (diffusion, exact solution, and their derivatives)
//
//
// Author: Jerome Droniou (jerome.droniou@monash.edu)
//

#include "TestCase.hpp"

using namespace HArDCore2D;

// ----------------------------------------------------------------------------
//                          Implementation
// ----------------------------------------------------------------------------

// Class - we set up the anisotropy when we create the class
TestCase::TestCase(std::vector<int> iTC)
    : m_iTC(iTC),
      _deg_diff(0),
      _lambda(pow(10, 2))
{
     if (m_iTC.size() == 2)
     {
          // If only two entries passed to TestCase, force advec and reac to be zero
          m_iTC.push_back(0);
          m_iTC.push_back(0);
     }
     validate();
     if (m_iTC[1] == 2 || m_iTC[1] == 4)
     {
          _deg_diff = 2;
     }
}

/////////////////// SOLUTION ///////////////////////////////:

// Solution
FType<double> TestCase::sol()
{
     FType<double> u = [&](const VectorRd x) -> double {
          return 0;
     };

     switch (m_iTC[0])
     {
          /// iTC[0]=1: \f$u(x,y)=sin(\pi x)  sin(\pi y)\f$
     case 1:
          u = [&](const VectorRd x) -> double {
               return sin(pi * x(0)) * sin(pi * x(1));
          };
          break;
          /// iTC[0]=2: \f$u(x,y)=cos(\pi x)  cos(\pi y)\f$
     case 2:
          u = [&](const VectorRd x) -> double {
               return cos(pi * x(0)) * cos(pi * x(1));
          };
          break;
          /// iTC[0]=3: \f$u(x,y)= x\f$
     case 3:
          u = [&](const VectorRd x) -> double {
               return x(0);
          };
          break;
          /// iTC[0]=4: \f$u(x,y)= y\f$
     case 4:
          u = [&](const VectorRd x) -> double {
               return x(1);
          };
          break;
          /// iTC[0]=5: \f$u(x,y)= x^2 + y^2\f$
     case 5:
          u = [&](const VectorRd x) -> double {
               return pow(x(0), 2) + pow(x(1), 2);
          };
          break;
          /// iTC[0]=6: \f$u(x,y)= e^x \sin(\pi x) \sin(\pi y)\f$
     case 6:
          u = [&](const VectorRd x) -> double {
               return exp(x(0)) * sin(pi * x(0)) * sin(pi * x(1));
          };
          break;
          /// iTC[0]=7: \f$u(x,y)= \frac{1}{4} \sin(2 \pi x) \sin(2 \pi y)\f$
     case 7:
          u = [&](const VectorRd x) -> double {
               return 0.25 * sin(2 * pi * x(0)) * sin(2 * pi * x(1));
          };
          break;
          /// iTC[0]=8: \f$u(x,y)= 5 x y (x-1) (y-1)\f$
     case 8:
          u = [&](const VectorRd x) -> double {
               return 5 * x(0) * x(1) * (x(0) - 1) * (x(1) - 1);
          };
          break;
          /// iTC[0]=9: \f$u(x,y)= \frac{1}{2} - \frac{1}{2}e^{\sin(\pi x) \sin(\pi y)}\f$
     case 9:
          u = [&](const VectorRd x) -> double {
               return 0.5 - 0.5 * exp(sin(pi * x(0)) * sin(pi * x(1)));
          };
          break;
          /// iTC[0]=10: \f$u(x,y)= 1 - e^{\frac{1}{4} \sin(2 \pi x) \sin(2 \pi y)}\f$
     case 10:
          u = [&](const VectorRd x) -> double {
               return 1 - exp(0.25 * sin(2 * pi * x(0)) * sin(2 * pi * x(1)));
          };
          break;
          /// iTC[0]=11: \f$u(x,y)= 5 - 5e^{x y (x-1) (y-1)}\f$
     case 11:
          u = [&](const VectorRd x) -> double {
               return 5 - 5 * exp(x(0) * x(1) * (x(0) - 1) * (x(1) - 1));
          };
          break;
     default:
          break;
     }

     return u;
}

// Gradient of the solution
CellFType<VectorRd> TestCase::grad_sol()
{
     CellFType<VectorRd> G = [&](const VectorRd x, const Cell *cell) -> VectorRd {
          return VectorRd::Zero();
     };

     switch (m_iTC[0])
     {
     case 1:
          G = [&](const VectorRd x, const Cell *cell) -> VectorRd {
               VectorRd g(cos(pi * x(0)) * sin(pi * x(1)), sin(pi * x(0)) * cos(pi * x(1)));
               return pi * g;
          };
          break;
     case 2:
          G = [&](const VectorRd x, const Cell *cell) -> VectorRd {
               VectorRd g(sin(pi * x(0)) * cos(pi * x(1)), cos(pi * x(0)) * sin(pi * x(1)));
               return -pi * g;
          };
          break;
     case 3:
          G = [&](const VectorRd x, const Cell *cell) -> VectorRd {
               VectorRd g(1, 0);
               return g;
          };
          break;
     case 4:
          G = [&](const VectorRd x, const Cell *cell) -> VectorRd {
               VectorRd g(0, 1);
               return g;
          };
          break;
     case 5:
          G = [&](const VectorRd x, const Cell *cell) -> VectorRd {
               VectorRd g(x(0), x(1));
               return 2 * g;
          };
          break;
     case 6:
          G = [&](const VectorRd x, const Cell *cell) -> VectorRd {
               VectorRd g(sin(pi * x(0)) * sin(pi * x(1)) + pi * cos(pi * x(0)) * sin(pi * x(1)), pi * sin(pi * x(0)) * cos(pi * x(1)));
               return exp(x(0)) * g;
          };
          break;
     case 7:
          G = [&](const VectorRd x, const Cell *cell) -> VectorRd {
               VectorRd g(cos(2 * pi * x(0)) * sin(2 * pi * x(1)), sin(2 * pi * x(0)) * cos(2 * pi * x(1)));
               return 0.5 * pi * g;
          };
          break;
     case 8:
          G = [&](const VectorRd x, const Cell *cell) -> VectorRd {
               VectorRd g(x(1) * (x(1) - 1) * (2 * x(0) - 1), x(0) * (x(0) - 1) * (2 * x(1) - 1));
               return 5 * g;
          };
          break;
     case 9:
          G = [&](const VectorRd x, const Cell *cell) -> VectorRd {
               VectorRd g(cos(pi * x(0)) * sin(pi * x(1)), sin(pi * x(0)) * cos(pi * x(1)));
               return -0.5 * pi * exp(sin(pi * x(0)) * sin(pi * x(1))) * g;
          };
          break;
     case 10:
          G = [&](const VectorRd x, const Cell *cell) -> VectorRd {
               VectorRd g(cos(2 * pi * x(0)) * sin(2 * pi * x(1)), sin(2 * pi * x(0)) * cos(2 * pi * x(1)));
               return -0.5 * pi * exp(0.25 * sin(2 * pi * x(0)) * sin(2 * pi * x(1))) * g;
          };
          break;
     case 11:
          G = [&](const VectorRd x, const Cell *cell) -> VectorRd {
               VectorRd g(x(1) * (x(1) - 1) * (2 * x(0) - 1), x(0) * (x(0) - 1) * (2 * x(1) - 1));
               return -5 * exp(x(0) * x(1) * (x(0) - 1) * (x(1) - 1)) * g;
          };
          break;
     default:
          break;
     }
     return G;
}

// Hessian of the solution
CellFType<MatrixRd> TestCase::hess_sol()
{
     CellFType<MatrixRd> H = [&](const VectorRd x, const Cell *cell) -> MatrixRd {
          return MatrixRd::Zero();
     };

     switch (m_iTC[0])
     {
     case 1:
          H = [&](const VectorRd x, const Cell *cell) -> MatrixRd {
               MatrixRd h;
               double a = -sin(pi * x(0)) * sin(pi * x(1));
               double b = cos(pi * x(0)) * cos(pi * x(1));
               h.row(0) << a, b;
               h.row(1) << b, a;
               return pi * pi * h;
          };
          break;
     case 2:
          H = [&](const VectorRd x, const Cell *cell) -> MatrixRd {
               MatrixRd h;
               double a = -cos(pi * x(0)) * cos(pi * x(1));
               double b = sin(pi * x(0)) * sin(pi * x(1));
               h.row(0) << a, b;
               h.row(1) << b, a;
               return pi * pi * h;
          };
          break;
     case 3:
          break;
     case 4:
          break;
     case 5:
          H = [&](const VectorRd x, const Cell *cell) -> MatrixRd {
               return 2 * MatrixRd::Identity();
          };
          break;
     case 6:
          H = [&](const VectorRd x, const Cell *cell) -> MatrixRd {
               MatrixRd h;
               double a = sin(pi * x(1)) * ((1 - pi * pi) * sin(pi * x(0)) + 2 * pi * cos(pi * x(0)));
               double b = pi * cos(pi * x(1)) * (sin(pi * x(0)) + pi * cos(pi * x(0)));
               double c = -pi * pi * sin(pi * x(0)) * sin(pi * x(1));
               h.row(0) << a, b;
               h.row(1) << b, c;
               return exp(x(0)) * h;
          };
          break;
     case 7:
          H = [&](const VectorRd x, const Cell *cell) -> MatrixRd {
               MatrixRd h;
               double a = -sin(2 * pi * x(0)) * sin(2 * pi * x(1));
               double b = cos(2 * pi * x(0)) * cos(2 * pi * x(1));
               h.row(0) << a, b;
               h.row(1) << b, a;
               return pi * pi * h;
          };
          break;
     case 8:
          H = [&](const VectorRd x, const Cell *cell) -> MatrixRd {
               MatrixRd h;
               double a = 2 * x(1) * (x(1) - 1);
               double b = (2 * x(0) - 1) * (2 * x(1) - 1);
               double c = 2 * x(0) * (x(0) - 1);
               h.row(0) << a, b;
               h.row(1) << b, c;
               return 5 * h;
          };
          break;
     case 9:
          H = [&](const VectorRd x, const Cell *cell) -> MatrixRd {
               MatrixRd h;
               double a = sin(pi * x(1)) * (-sin(pi * x(0)) + cos(pi * x(0)) * cos(pi * x(0)) * sin(pi * x(1)));
               double b = cos(pi * x(0)) * cos(pi * x(1)) * (1 + sin(pi * x(0)) * sin(pi * x(1)));
               double c = sin(pi * x(0)) * (-sin(pi * x(1)) + cos(pi * x(1)) * cos(pi * x(1)) * sin(pi * x(0)));
               h.row(0) << a, b;
               h.row(1) << b, c;
               return -0.5 * pi * pi * exp(sin(pi * x(0)) * sin(pi * x(1))) * h;
          };
          break;
     case 10:
          H = [&](const VectorRd x, const Cell *cell) -> MatrixRd {
               MatrixRd h;
               double a = sin(2 * pi * x(1)) * (-4 * sin(2 * pi * x(0)) + cos(2 * pi * x(0)) * cos(2 * pi * x(0)) * sin(2 * pi * x(1)));
               double b = cos(2 * pi * x(0)) * cos(2 * pi * x(1)) * (4 + sin(2 * pi * x(0)) * sin(2 * pi * x(1)));
               double c = sin(2 * pi * x(0)) * (-4 * sin(2 * pi * x(1)) + cos(2 * pi * x(1)) * cos(2 * pi * x(1)) * sin(2 * pi * x(0)));
               h.row(0) << a, b;
               h.row(1) << b, c;
               return -0.25 * pi * pi * exp(0.25 * sin(2 * pi * x(0)) * sin(2 * pi * x(1))) * h;
          };
          break;
     case 11:
          H = [&](const VectorRd x, const Cell *cell) -> MatrixRd {
               MatrixRd h;
               double a = x(1) * (x(1) - 1) * (2 + x(1) * (x(1) - 1) * (2 * x(0) - 1) * (2 * x(0) - 1));
               double b = (2 * x(0) - 1) * (2 * x(1) - 1) * (1 + x(0) * x(1) * (x(0) - 1) * (x(1) - 1));
               double c = x(0) * (x(0) - 1) * (2 + x(0) * (x(0) - 1) * (2 * x(1) - 1) * (2 * x(1) - 1));
               h.row(0) << a, b;
               h.row(1) << b, c;
               return -5 * exp(x(0) * x(1) * (x(0) - 1) * (x(1) - 1)) * h;
          };
          break;
     default:
          break;
     }
     return H;
}

//////////////////////////// DIFFUSION /////////////////////////////

// Diffusion matrix
CellFType<MatrixRd> TestCase::diff()
{
     CellFType<MatrixRd> K = [&](const VectorRd x, const Cell *cell) -> MatrixRd {
          return MatrixRd::Identity();
     };
     switch (m_iTC[1])
     {
          /// iTC[1]=1: Diff = Id
     case 1:
          break;
          /// iTC[1]=2: Diff = \f$\left[\begin{array}{cc}y^2+1 & -xy\\  -xy & x^2+1\end{array}\right]\f$
     case 2:
          K = [&](const VectorRd x, const Cell *cell) -> MatrixRd {
               MatrixRd k;
               k.row(0) << pow(x(1), 2) + 1, -x(0) * x(1);
               k.row(1) << -x(0) * x(1), pow(x(0), 2) + 1;
               return k;
          };
          break;
          /// iTC[1]=3: Diff=\f$\left[\begin{array}{cc}\lambda & 0\\ 0 & 1\end{array}\right]\f$ if \f$y<1/2\f$, Diff=Id if \f$y\ge 1/2\f$. Only valid with iTC[0]=2 (\f$\partial_x u\f$ must vanish along \f$y=1/2\f$).
     case 3:
          K = [&](const VectorRd x, const Cell *cell) -> MatrixRd {
               MatrixRd k = MatrixRd::Identity();
               if (cell->center_mass().y() < 0.5)
               {
                    k.row(0) << _lambda, 0;
               }
               return k;
          };
          break;
          /// iTC[1]=4: rotating diffusion. Diff=\f$\left[\begin{array}{cc}\epsilon \bar{x}^2 +  \bar{y}^2 & (\epsilon-1)\bar{x}\bar{y}\\ (\epsilon-1)\bar{x}\bar{y} & \bar{x}^2+\epsilon \bar{y}^2\end{array}\right]\f$, where \f$\bar{x}=x+0.1\f$ and \f$\bar{y}=y+0.1\f$.
     case 4:
     {
          K = [&](const VectorRd x, const Cell *cell) -> MatrixRd {
               MatrixRd k;
               double barx = x(0) + 0.1;
               double bary = x(1) + 0.1;
               k.row(0) << eps * barx * barx + bary * bary, (eps - 1) * barx * bary;
               k.row(1) << (eps - 1) * barx * bary, barx * barx + eps * bary * bary;
               return k;
          };
          break;
     }
          /// iTC[1]=5: Diff=\f$\left[\begin{array}{cc}\lambda & 0\\ 0 & 1\end{array}\right]\f$.
     case 5:
          K = [&](const VectorRd x, const Cell *cell) -> MatrixRd {
               MatrixRd k = MatrixRd::Identity();
               k.row(0) << _lambda, 0;
               return k;
          };
          break;
     case 6:
          K = [&](const VectorRd x, const Cell *cell) -> MatrixRd {
               MatrixRd k;
               k.row(0) << 2, 1;
               k.row(1) << 1, 1;
               return k;
          };
          break;
     default:
          break;
     }
     return K;
}

// Divergence by row of the diffusion matrix
CellFType<VectorRd> TestCase::div_diff()
{
     CellFType<VectorRd> divK = [&](const VectorRd x, const Cell *cell) -> VectorRd {
          return VectorRd::Zero();
     };

     switch (m_iTC[1])
     {
     case 1:
          break;
     case 2:
          divK = [&](const VectorRd x, const Cell *cell) -> VectorRd {
               VectorRd divk(-x(0), -x(1));
               return divk;
          };
          break;
     case 3:
          break;
     case 4:
     {
          divK = [&](const VectorRd x, const Cell *cell) -> VectorRd {
               double barx = x(0) + 0.1;
               double bary = x(1) + 0.1;
               VectorRd divk(2 * eps * barx + (eps - 1) * barx, (eps - 1) * bary + 2 * eps * bary);
               return divk;
          };
          break;
     }
     case 5:
          break;
     case 6:
          break;
     default:
          break;
     }
     return divK;
}

// Advection term
CellFType<VectorRd> TestCase::advec()
{
     // Default must be zero
     CellFType<VectorRd> A = [&](const VectorRd x, const Cell *cell) -> VectorRd {
          VectorRd a(0, 0);
          return a;
     };
     switch (m_iTC[2])
     {
     case 1:
          A = [&](const VectorRd x, const Cell *cell) -> VectorRd {
               VectorRd a(1, 1);
               return a;
          };
          break;
     case 2:
          A = [&](const VectorRd x, const Cell *cell) -> VectorRd {
               return x;
          };
          break;
     case 3:
          A = [&](const VectorRd x, const Cell *cell) -> VectorRd {
               VectorRd a(x(1), x(0));
               return a;
          };
          break;
     case 4:
          A = [&](const VectorRd x, const Cell *cell) -> VectorRd {
               VectorRd a(pow(x(0), 2) + pow(x(1), 2), pow(x(0), 2) + pow(x(1), 2));
               return a;
          };
          break;
     case 5:
          A = [&](const VectorRd x, const Cell *cell) -> VectorRd {
               VectorRd a(pow(eps, -1), 0);
               return a;
          };
          break;
     case 6:
          A = [&](const VectorRd x, const Cell *cell) -> VectorRd {
               VectorRd a(exp(x(0)), exp(x(1)));
               return a;
          };
          break;
     default:
          break;
     }
     return A;
}

// Divergence of advection term
CellFType<double> TestCase::div_advec()
{
     CellFType<double> divA = [&](const VectorRd x, const Cell *cell) -> double {
          return 0;
     };
     switch (m_iTC[2])
     {
     case 1:
          break;
     case 2:
          divA = [&](const VectorRd x, const Cell *cell) -> double {
               return 2;
          };
          div_advec_zero = false;
          break;
     case 3:
          break;
     case 4:
          divA = [&](const VectorRd x, const Cell *cell) -> double {
               return 2 * x(0) + 2 * x(1);
          };
          div_advec_zero = false;
          div_advec_const = false;
          break;
     case 5:
          break;
     case 6:
          divA = [&](const VectorRd x, const Cell *cell) -> double {
               return exp(x(0)) + exp(x(1));
          };
          div_advec_zero = false;
          div_advec_const = false;
          break;
     default:
          break;
     }
     return divA;
}

// Reaction term
CellFType<double> TestCase::reac()
{
     // Default must be zero
     CellFType<double> R = [&](const VectorRd x, const Cell *cell) -> double {
          return 0;
     };
     switch (m_iTC[3])
     {
     case 1:
          R = [&](const VectorRd x, const Cell *cell) -> double {
               return 1;
          };
          break;
     case 2:
          R = [&](const VectorRd x, const Cell *cell) -> double {
               if (cell->center_mass().y() < 0.5)
               {
                    return x(0) + eps;
               }
               else
               {
                    return x(1) + eps;
               }
          };
          reac_const = false;
          break;
     case 3:
          R = [&](const VectorRd x, const Cell *cell) -> double {
               return x(0) + x(1) + eps;
          };
          reac_const = false;
          break;
     case 4:
          R = [&](const VectorRd x, const Cell *cell) -> double {
               return exp(x(0) * x(1));
          };
          reac_const = false;
          break;
     case 5:
     {
          R = [&](const VectorRd x, const Cell *cell) -> double {
               return sin(pi * x(0)) * sin(pi * x(1)) + eps;
          };
          reac_const = false;
          break;
     }
     case 6:
          R = [&](const VectorRd x, const Cell *cell) -> double {
               double xcoord = cell->center_mass().x();
               double ycoord = cell->center_mass().y();
               if ((xcoord < 0.5 && ycoord < 0.5) || (xcoord > 0.5 && ycoord > 0.5))
               {
                    return (0.5 - xcoord) * (0.5 - ycoord) + eps;
               }
               else
               {
                    return -(0.5 - xcoord) * (0.5 - ycoord) + eps;
               }
          };
          reac_const = false;
          break;
     case 7:
          R = [&](const VectorRd x, const Cell *cell) -> double {
               return eps;
          };
          break;
     case 8:
          R = [&](const VectorRd x, const Cell *cell) -> double {
               return pow(eps, -1);
          };
          break;
     default:
          break;
     }
     return R;
}

///////////////////////////// SOURCE TERM ///////////////////////////

// Diffusion source term
CellFType<double> TestCase::diff_source()
{
     CellFType<double> f = [&](const VectorRd x, const Cell *cell) -> double {
          MatrixRd AHu = diff()(x, cell) * hess_sol()(x, cell);
          return -AHu.trace() - div_diff()(x, cell).dot(grad_sol()(x, cell));
     };

     return f;
}

// Diffusion-Advection-Reaction source term
CellFType<double> TestCase::diff_advec_reac_source()
{
     CellFType<double> f = [&](const VectorRd x, const Cell *cell) -> double {
          return diff_source()(x, cell) + div_advec()(x, cell) * sol()(x) + advec()(x, cell).dot(grad_sol()(x, cell)) + reac()(x, cell) * sol()(x);
     };

     return f;
}

///////////////////////////// VALIDATION ////////////////////////////

void TestCase::validate()
{

     if (m_iTC[0] > 11 || m_iTC[0] < 1 || m_iTC[1] > 6 || m_iTC[1] < 1 || m_iTC[2] > 6 || m_iTC[2] < 0 || m_iTC[3] > 8 || m_iTC[3] < 0 || (m_iTC[1] == 3 && m_iTC[0] != 2))
     {
          std::cout << "Incorrect choice of test cases: iTC = " << m_iTC[0] << ", " << m_iTC[1] << ", " << m_iTC[2] << ", " << m_iTC[3] << "\n";
          exit(EXIT_FAILURE);
     }
}
