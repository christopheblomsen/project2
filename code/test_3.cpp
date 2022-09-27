//g++ test_3.cpp src/utils.cpp -I include -larmadillo -o test_3

#include <iostream>
#include <armadillo>
#include "utils.hpp"

int main(){
  //define a matrix to test the function max_offdiag_symmetric()
  int n = 4;
  arma::mat A = arma::mat(n, n);
  arma::mat R = arma::mat(n, n);
   R.eye();
  A.eye();
  A(2, 1) = -0.7;
  A(1, 2) = -0.7;
  A(0, 3) = 0.5;
  A(3, 0) = 0.5;
  std::cout << A << std::endl;
  int k;
  int l;
  double m = max_offdiag_symmetric(A, k, l);

  //test jacobi rotation
  jacobi_rotate(A, R, k, l);
  std::cout << A<< std::endl;
  std::cout << R << std::endl;
  return 0;
}
