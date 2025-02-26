#include "Sturm.hh"
#include <iostream>

int main()
{
  Sturm::Poly<double> p1(3); p1 << 1.0, -3.0, 2.0; // p1(x) = 1 - 3x + 2x^2
  std::cout << "p1(x) = " << p1 << std::endl;
  Sturm::Sequence<double> seq(p1);
  std::cout << "seq = " << seq << std::endl;
  return 0;
}
