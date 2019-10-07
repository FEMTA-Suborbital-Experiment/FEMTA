

#include "mathematics.h"


compiler_constant int gcd(int a, int b) {
  // greatest common divisor
  // algorithm adapted from Dr. Samuel Wagstaff
  
  int g = a;
  int t = b;
  
  while (t > 0) {
    int q = g / t;
    int w = g - q * t;
    g = t;
    t = w;
  }
  
  return g;
}
