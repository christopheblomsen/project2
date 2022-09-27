#include "utils.hpp"
#include <iostream>
#include <armadillo>
#include <cmath>
#include <iomanip>

//calculate analytical eigenvalues, n is number of points, d is diagonal, and a is super/sub-diagonal
arma::vec analytical_eigenvalues(int n, double d, double a){
  double h = 1./n;
  d = d/std::pow(h, 2);
  a = a/std::pow(h, 2);
  arma::vec i = arma::regspace(1, n);
  arma::vec L = d + 2.*a*arma::cos((i*M_PI)/(n+1));
  return L;
}



//calculate analytical eigenvectors
arma::mat analytical_eigenvectors(int n){
  arma::mat u(n, n);
  arma::vec i = arma::regspace(1, n);
  for (int j = 0; j<n; j++){
    u(arma::span::all, j) = arma::sin((i*(j+1)*M_PI)/(n+1));
  }
  return u;
}

//Sorting  algorithm modifies matrix A and corresponding eigenvalues lam to follow ascending order
void sorted_matrix(arma::mat &A, arma::vec &lam){
  // A is the matrix we sent to sort the column considering the elements of lam,
  // We get the index in descending order of the vector lam if we want in ascending 
  // just take off "descend"
  arma::uvec indices = arma::sort_index(lam, "ascend"); // Code for task 5 and 6 works with ascending order
  
  double M = lam.size();

  // New sorted matrix B and vector b
  arma::mat B = arma::mat(M, M);
  arma::vec b = arma::vec(M);

  for (int i = 0; i < M; i++){
    B(arma::span::all, i) = A(arma::span::all, indices(i));
    b(i) = lam(indices(i));
  }
  
  //modifying matrix and vector
  A = B;
  lam = b;
}


//find the largest (in absolute value) element of an arma matrix
//return the value and write the indecies to memory
//
double max_offdiag_symmetric(const arma::mat &A, int& k, int& l){
  int n = A.n_rows; //assume nxn symmetric matrix
 double b = 0; //biggest absolute value
  // iterate over upper triangle
  for (int i = 0; i<n; i++){
    for (int j = i+1; j<n; j++){
        if (std::abs(A(i, j)) > std::abs(b)){
          b = A(i, j);
          k = i;
          l = j;
        }
    }
  }
  double m = A(k, l);
  return m;
}



//give size of matrix n, and coefficients a, d, e which signatures the diagonal signature
//since h varies with matrix size we dont pass a/h^2 as an argument
arma::mat create_tridiagonal(const int n, const double a, const double d, const double e)
{   
    //assume domain from 0 to 1 so that h=1/n
    double h = 1./n;
    //std::cout <<"stepsize:" <<  h << std::endl;
    // Initialize the matrix and its diagonal
    arma::mat A = arma::mat(n, n, arma::fill::eye) * d/std::pow(h, 2);

    // Iteration to fill a general tridiagonal
    for (int i=0; i < n-1; i++){
        A(i, i+1) = e/std::pow(h, 2);
        A(i+1, i) = a/std::pow(h, 2);
    }
    return A;
}




// perform one jacobi rotation
void jacobi_rotate(arma::mat &A, arma::mat &R, int &k, int &l){
  // performs steps 3.1 to 3.4 as presented in lectures
  //updates A and R to A(m+1) and R(m+1)
  //
  //calculate tau and t, c, s
  int n = A.n_rows;
  double tau = (A(l, l) - A(k, k))/(2*A(k, l));
  double t;
  if (tau >= 0){
    t = 1/(tau + std::sqrt(1+ tau*tau));
    }
    else{
      t = -1/(-tau + std::sqrt(1 + tau*tau));
    }
  //std::cout << tau << std::endl;
  double c = 1/std::sqrt(1 + t*t);
  double s = c*t;

  // update values at kk, kl, lk, ll
  double a_kk = A(k, k);
  A(k, k) = A(k, k)*std::pow(c, 2) - 2*A(k, l)*c*s + A(l, l)*std::pow(s, 2);
  A(l, l) = A(l, l)*std::pow(c, 2) + 2*A(k, l)*c*s + a_kk*std::pow(s, 2);
  A(k, l) = 0;
   A(l, k) = 0;

   //update rows and cols k, l
   double a_ik = 0;
   for (int i = 0; i < n; i++){
     if (i != k and i != l){
       a_ik = A(i, k); //save for later
       A(i, k) = A(i, k)*c - A(i, l)*s;
        A(k, i) = A(i, k);
        A(i, l) = A(i, l)*c + a_ik*s;
        A(l, i) = A(i, l);
     }
   }

   //update R
   double r_ik = 0;
   for (int i = 0; i < n; i++){
     r_ik = R(i, k);
     R(i, k) = R(i, k)*c - R(i, l)*s;
      R(i, l) = R(i, l)*c + r_ik*s;
   }
}



// Returns the eigenvalues and eigenvectors of A, in ascending order of eigenvalues

//solve the eigenvalue problem
void jacobi_eigensolver(const arma::mat& A, double eps, arma::vec& eigenvalues, arma::mat& eigenvectors,
                        const int maxiter, int& iterations, bool& converged){
 int k = 0;
  int l = 0;
  arma::mat D = A; // this is not efficient but might be necessary if we want to keep A const
  arma::mat R(arma::size(A), arma::fill::zeros);
  R.eye();

  //find the indecies of largest offdiag element
  double a_max = max_offdiag_symmetric(A, k, l);

  //loop to perform the rotations
  while (std::abs(a_max) > eps and iterations <= maxiter){
    jacobi_rotate(D, R, k, l); // changed D to be A instead
    a_max = max_offdiag_symmetric(D, k, l); // --//--
    iterations +=1;
  }

  // Store solution
  eigenvalues = arma::diagvec(D);
  eigenvectors = R;
  
  // Sorting the eigenvalues and eigenvectors with ascending order
  sorted_matrix(eigenvectors, eigenvalues);
  
  //update convergence
  if (iterations < maxiter){
    converged = true;
  }

}


//resetting variables that are modifide during simulations
void reset_simulation(int& n, arma::mat& A, arma::mat& R, int& iterations, bool &converged,
                      arma::vec& eigenvalues, arma::mat& eigenvectors){
  A  = create_tridiagonal(n, -1, 2, -1);
  R = arma::mat(n, n);
  R.eye();
  iterations = 0;
  converged = false;
  eigenvalues.set_size(n);
   eigenvectors.set_size(n, n);
  eigenvalues.fill(0);
   eigenvectors.fill(0);
}

//add boundary values and return the first three modes
arma::mat add_boundary_vals(int &n, arma::mat &eigenvectors){
  arma::rowvec zeros = arma::rowvec(n, arma::fill::zeros);
  eigenvectors.insert_rows(0, zeros);
  eigenvectors.insert_rows(n+1, zeros);
  arma::mat first_modes = eigenvectors.cols(0, 2);
  return first_modes;
}


// Return a string with double in scientific notation
/*
** Parameters:
** -----------
** d
**     the number to be printed
** width
**     the width of the whitespace
** prec
**     the precision of the number
 */
std::string sci_format(const double d, const int width, const int prec){
    std::stringstream ss;
    ss << std::setw(width) << std::setprecision(prec) << std::scientific << d;
    return ss.str();
}

// Makes an output file
/*
** Parameters:
** -----------
** filename
**         Name of the file
** v
**         v vector
** x
**         corresponding x vector to v
** n
**         length of the vectors
** width
**         whitespace width in output
*p* prec
**         precision in output
 */
void output(std::string filename, const arma::vec& x, const arma::vec& u,
            const int n, int width, int prec){
    std::ofstream outfile;
    outfile.open(filename);
    // Header for the file
    outfile << "  x" << std::setw(width) << "u" << std::setw(width) << std::endl;
    // Couldn't get a pretty output without the loop, might revisit
    for (int i=0; i<n; i++){
        outfile << sci_format(x[i], width, prec) << sci_format(u[i], width, prec) << std::endl;
    }
    outfile.close();
}
