#include "../src/householder.hpp"
#include <iostream>

int main() {
  mtrx::dense<double> A(3, 3, {1, 3, 2, 4, 2, 6, 3, 5, 1});
  A.print();
  hldr::hh_QR(A).first.print();
  std::cout << std::endl;
  hldr::hh_QR(A).second.print();
  std::cout << std::endl;
  std::cout << hldr::QR_solve(A, {2, 3, 1}) << std::endl;
}