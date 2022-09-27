#include <armadillo>

//calculate the analytic eigenvals and return as a vector
arma::vec analytical_eigenvalues(int n, double d, double a);

//calculate analytic eig.vecs and return as a matrix, columns are the eigvecs
arma::mat analytical_eigenvectors(int n );

// Sort the column of a matrix according to the absolute value of an array 
// in descending order
void sorted_matrix(const arma::mat &A, const arma::vec &lam);

//find the largest off diagonal element in an arma matrix
//return the value and write the indecies to memory
double max_offdiag_symmetric(const arma::mat &A, int &k , int &l );


// Create a general tridiagonal matrix of size NxN
arma::mat create_tridiagonal(const int n, const double a, const double d, const double e);

//perform one jacobi rotation
void jacobi_rotate(arma::mat &A, arma::mat &R, int &k, int &l);

//solve the eigenvalue problem
void jacobi_eigensolver(const arma::mat& A, double eps, arma::vec& eigenvalues, arma::mat& eigenvectors,
                        const int maxiter, int& iterations, bool& converged);

//reset simulation variables
void reset_simulation(int& n, arma::mat& A, arma::mat& R, int& iterations, bool &converged,
                      arma::vec& eigenvalues, arma::mat& eigenvectors);

//add boundary vals and return first three modes
arma::mat add_boundary_vals(int &n, arma::mat &eigenvectors);

//return a string with scientific format
std::string sci_format(const double d, const int width, const int prec);

// output
void output(std::string filename, const arma::vec& x, const arma::vec& u,
           const int n, int width, int prec);
