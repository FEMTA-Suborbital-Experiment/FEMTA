
int gcd(int a, int b) {
  // returns greatest common divisor, as adapted from Dr.
  // Samuel Wagstaff's CS 355 course and discovered by Euclid.
  
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

int lcm(int a, int b) {          // the least common multiple
  return (a * b) / gcd(a, b);    // relationship between gcd and lcm
}
