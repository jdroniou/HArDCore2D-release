#ifndef _MATH_HPP
#define _MATH_HPP

#include <cmath>

namespace Math
{
     /// Free math functions and global variables ///

  inline constexpr double PI = M_PI;

     inline const size_t factorial(size_t n)
     {
          return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
     }

     inline const size_t nChoosek(size_t n, size_t k)
     {
          if (k > n)
          {
               return 0;
          }
          if (k * 2 > n)
          {
               k = n - k;
          }
          if (k == 0)
          {
               return 1;
          }

          size_t result = n;
          for (size_t i = 2; i <= k; ++i)
          {
               result *= (n - i + 1);
               result /= i;
          }
          return result;
     }

     template <typename T>
     const int sgn(T val)
     {
          return (T(0) < val) - (val < T(0));
     }

     // branch -1: -2\pi < \theta \le 0
     // branch 0:  -\pi < \theta \le \pi
     // branch 1:  0 \le \theta < 2\pi
     inline double atan2(double y, double x, int branch = 0)
     {
          double val = (((x != 0) || (y != 0)) ? std::atan2(y, x) : 0.0);
          if (branch == 1)
          {
               if (y < 0)
               {
                    val += 2.0 * PI;
               }
          }
          if (branch == -1)
          {
               if (y > 0)
               {
                    val -= 2.0 * PI;
               }
               else if ((y == 0) && (x < 0))
               {
                    val -= 2.0 * PI;
               }
          }
          return val;
     }
}

#endif
