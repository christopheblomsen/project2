#include "utils.hpp"
#include <armadillo>
#include <iomanip>

int main(){

  //Problem 2
  arma::mat A_2 = create_tridiagonal(6, -1, 2, -1);
  
  std::cout << "Problem 2" << std::endl;
  std::cout << "Tridiagonal matrix A \n";
  std::cout << A_2 << std::endl;

  // Armadillo solution using eig_sym function
  arma::mat eigvector_2;
  arma::vec eigvalue_2;

  arma::eig_sym(eigvalue_2, eigvector_2, A_2);

  //Analytical solution 
  arma::vec eig_values_2 = analytical_eigenvalues(6, 2, -1);
  arma::mat eig_vectors_2 = analytical_eigenvectors(6);

  std::cout << "Analytical solution\n";
  std::cout << arma::normalise(eig_values_2) << std::endl;
  std::cout << arma::normalise(eig_vectors_2) << std::endl;

  std::cout << "Armadillo solution\n";
  std::cout << arma::normalise(eigvalue_2) <<std::endl;
  std::cout << eigvector_2 << std::endl;
  std::cout << "End Problem 2\n";

  //problem 4
  int n = 6;
  arma::mat A = create_tridiagonal(n, -1, 2, -1);
  arma::mat R = arma::mat(n, n);
  R.eye();

  //declare paramaters for simulations
  double eps = 1e-8;
  int maxiter = 10000;
  int iterations = 0;
  bool converged = false;
  arma::vec eigenvalues(n);
  arma::mat eigenvectors(n, n);
  jacobi_eigensolver(A, eps, eigenvalues, eigenvectors, maxiter, iterations, converged);

  //analytical solution
 arma::vec L = analytical_eigenvalues(n, 2, -1);
 arma::mat u = analytical_eigenvectors(n);

 //save results to compare solutions in python-script
 eigenvectors.save("eig_vec_num_4.txt", arma::raw_ascii);
 u.save("eig_vec_ana_4.txt", arma::raw_ascii);


  //problem 5
 int n_max = 10; //actual max is multiplied by 10
 int width = 24;
 //printing to screen
 std::cout << std::setw(width) << "N" << std::setw(width) << "Iterations" << std::setw(width) <<
   "log10(Iter)/log10(N)" << std::setw(width) << "Converged" << std::endl;

 //looping to test different n
 for (int i  = 1; i<=n_max; i++ ){
   n = 10*i;
   //resetting simulation parameters
  reset_simulation(n, A, R,iterations, converged, eigenvalues, eigenvectors);
  jacobi_eigensolver(A, eps, eigenvalues, eigenvectors, maxiter, iterations, converged);
  double scaling = std::log10(iterations)/std::log10(n);
  //printing results to screen
  std::cout  << std::setw(width) << n << std::setw(width) <<  iterations << std::setw(width) <<
    scaling << std::setw(width) << converged << std::endl;
  }

//problem 6
//simulating for 10 steps
   n = 9;
   reset_simulation(n, A, R,iterations, converged, eigenvalues, eigenvectors);
   jacobi_eigensolver(A, eps, eigenvalues, eigenvectors, maxiter, iterations, converged);

  //adding boundary points to eigenvectors
  arma::mat  first_modes = add_boundary_vals(n, eigenvectors);
  arma::vec x = arma::linspace(0, 1, n+2);

  first_modes.save("eig_vec_num_6a.txt", arma::raw_ascii);
  x.save("x_axis_6a.txt", arma::raw_ascii);

  //analytical solution
  arma::mat u_10 = analytical_eigenvectors(n);
  arma::mat u_10_bound = add_boundary_vals(n, u_10);
  u_10_bound = arma::normalise(u_10_bound);
  u_10_bound(arma::span::all, 1) *= -1;

  u_10_bound.save("eig_vec_ana_6a.txt", arma::raw_ascii);

//simulating with 100 steps
   n = 99;
   reset_simulation(n, A, R,iterations, converged, eigenvalues, eigenvectors);
   jacobi_eigensolver(A, eps, eigenvalues, eigenvectors, maxiter, iterations, converged);

  //adding boundary points to eigenvectors
  first_modes = add_boundary_vals(n, eigenvectors);
  x = arma::linspace(0, 1, n+2);

  first_modes.save("eig_vec_num_6b.txt", arma::raw_ascii);
  x.save("x_axis_6b.txt", arma::raw_ascii);

  //analytical solution
  arma::mat u_100 = analytical_eigenvectors(n);
  arma::mat u_100_bound = add_boundary_vals(n, u_100);
  u_100_bound = arma::normalise(u_100_bound);
  u_100_bound.save("eig_vec_ana_6b.txt", arma::raw_ascii);

  return 0;
}
