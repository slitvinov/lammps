#ifndef LMP_MATH_SPECIAL_H
#define LMP_MATH_SPECIAL_H 
#include <cmath>
namespace LAMMPS_NS {
namespace MathSpecial {
  extern double factorial(const int n);
  extern double exp2_x86(double x);
  extern double fm_exp(double x);
  extern double erfcx_y100(const double y100);
  static inline double my_erfcx(const double x)
  {
    if (x >= 0.0)
      return erfcx_y100(400.0 / (4.0 + x));
    else
      return 2.0 * exp(x * x) - erfcx_y100(400.0 / (4.0 - x));
  }
  static inline double expmsq(double x)
  {
    x *= x;
    x *= 1.4426950408889634074;
#if defined(__BYTE_ORDER__) && __BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__
    return (x < 1023.0) ? exp2_x86(-x) : 0.0;
#else
    return (x < 1023.0) ? exp2(-x) : 0.0;
#endif
  }
  static inline double square(const double &x)
  {
    return x * x;
  }
  static inline double cube(const double &x)
  {
    return x * x * x;
  }
  static inline double powsign(const int n)
  {
    return (n & 1) ? -1.0 : 1.0;
  }
  static inline double powint(const double &x, const int n)
  {
    double yy, ww;
    if (n == 0) return 1.0;
    if (x == 0.0) return 0.0;
    int nn = (n > 0) ? n : -n;
    ww = x;
    for (yy = 1.0; nn != 0; nn >>= 1, ww *= ww)
      if (nn & 1) yy *= ww;
    return (n > 0) ? yy : 1.0 / yy;
  }
  static inline double powsinxx(const double &x, int n)
  {
    double yy, ww;
    if (x == 0.0) return 1.0;
    ww = sin(x) / x;
    for (yy = 1.0; n != 0; n >>= 1, ww *= ww)
      if (n & 1) yy *= ww;
    return yy;
  }
}
}
#endif
