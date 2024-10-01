#include <RcppArmadillo.h>
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]
using namespace arma;

mat crossprod(mat X, mat Y){
    return(X.t() * Y);
}

// [[Rcpp::export]]
mat projection(mat A){
  return(A * solve(crossprod(A, A), A.t()));
}


double plogis(double x){
  return(std::exp(x)/(1+std::exp(x)));
}

// [[Rcpp::export]]
mat modelMatrix(mat X, mat Y){
  int p = X.n_cols;
  int q = Y.n_cols;
  int n = X.n_rows;
  mat Xcopy = repmat(X, 1, q);
  mat Ycopy(n, p*q, fill::zeros);
  for(int k=0; k<q; k++){
    int h = (k/p);
    Ycopy.col(k) = Y.col(h);
  }
  mat S = Xcopy % Ycopy; 
  return(S);
}



// [[Rcpp::export]]
List extractMatrix(int j, Mat<int> grid, mat X)
{
  uvec ind_jj = find(grid(span::all, 0) == j && grid(span::all, 1) == j);
  uvec ind_jjprime1 = find(grid(span::all, 0) == j && grid(span::all, 1) !=j);
  uvec ind_jjprime2 = find(grid(span::all, 0) != j && grid(span::all, 1) == j);
  uvec ind_jjprime = join_vert(ind_jjprime1, ind_jjprime2);
  
  Mat<int> gridvec1 = grid.rows(ind_jjprime1);
  Col<int> gridij1 = gridvec1(span::all, 1);
  uvec grid1 = conv_to< uvec >::from(gridij1);
  mat Xij1 = X.cols(grid1);
  
  Mat<int> gridvec2 = grid.rows(ind_jjprime2);
  Col<int> gridij2 = gridvec2(span::all, 0);
  uvec grid2 = conv_to< uvec >::from(gridij2);
  mat Xij2 = X.cols(grid2);
  
  
  mat Xij = join_horiz(Xij1, Xij2);
  uvec ind_return = join_vert(grid1, grid2);
  
  return(Rcpp::List::create(
                            Rcpp::Named("ind_jj") = ind_jj,
                            Rcpp::Named("ind_jjprime") = ind_jjprime,
                            Rcpp::Named("Xminusj") = Xij
                            ));
}

// [[Rcpp::export]]
double pseudologLLH(mat X, mat fy, mat tau0, mat Gamma, mat C, Mat<int> grid){
    int n = fy.n_rows;
    int q = X.n_cols;
    mat vary_y = fy * C * Gamma.t();
    mat fC = fy * C;
    double llh = 0.0;
    for(int j=0; j<q; j++){
        List Xindex = extractMatrix(j, grid, X);
        for(int i=0; i < n; ++i){
            mat vary_y_i = vary_y.row(i);
            mat Xj = Xindex[2];
            uvec ind_jj = Xindex[0];
            uvec ind_jjprime = Xindex[1];
        
            mat Xij = Xj.row(i);
            
            double eps_ij = as_scalar(tau0.rows(ind_jj)) + 
              as_scalar(vary_y_i.cols(ind_jj)) + 
              as_scalar(crossprod(tau0.rows(ind_jjprime), Xij.t())) + 
              as_scalar(crossprod(vary_y_i.cols(ind_jjprime).t(), Xij.t()));
            
            //eps_ij = as_scalar(eps_ij);
            
            llh += X(i,j) * eps_ij - std::log(1 + std::exp(eps_ij));
        }
    }
        
     return(llh);
}

// [[Rcpp::export]]
mat globalupdateGamma(mat X, mat fy, mat currenttau0, mat currentGamma, mat currentC, Mat<int> grid){
    int d = currentGamma.n_cols;
    int q = X.n_cols;
    int kk = currentGamma.n_rows;
    int n = fy.n_rows;
    mat vary_y = fy * currentC * currentGamma.t();
    mat J(kk*d, kk*d, fill::zeros);
    mat S_Gamma(kk, 1, fill::zeros);
    mat fC = fy * currentC;
    for(int j=0; j<q; j++){
        List Xindex = extractMatrix(j, grid, X);
        for(int i=0; i < n; ++i){
            mat vary_y_i = vary_y.row(i);
            mat Xj = Xindex[2];
            uvec ind_jj = Xindex[0];
            uvec ind_jjprime = Xindex[1];
        
            mat Xij = Xj.row(i);
            
            double eps_ij = as_scalar(currenttau0.rows(ind_jj)) +
              as_scalar(vary_y_i.cols(ind_jj)) +
              as_scalar(crossprod(currenttau0.rows(ind_jjprime), Xij.t())) +
              as_scalar(crossprod(vary_y_i.cols(ind_jjprime).t(), Xij.t()));
           
            double probs_ij = plogis(eps_ij);
            double vars_ij = probs_ij*(1-probs_ij);
          
            S_Gamma.rows(ind_jj) += (X(i, j) - probs_ij) * fC.row(i);
            S_Gamma.rows(ind_jjprime) += (X(i, j) - probs_ij) * as_scalar(fC.row(i)) * Xij.t();
          
          mat onesmatrix = ones(1, 1);
          mat fiC = fC.row(i);
          mat M = modelMatrix(fiC, Xij);
          mat XmatrixInvolved = join_horiz(fiC, modelMatrix(fiC, Xij));  
          
          
          uvec tmp_index = join_vert(ind_jj, ind_jjprime);
          //tmp_index.print();
          J(tmp_index, tmp_index) -=  vars_ij * crossprod(XmatrixInvolved, XmatrixInvolved);
        }
      }
      mat newGamma = currentGamma - solve(J, S_Gamma);
      return(newGamma);
}


// [[Rcpp::export]]
mat globalupdateTau0(mat X, mat fy, mat currenttau0, mat currentGamma, mat currentC, Mat<int> grid){
  int d = currentGamma.n_cols;
  int q = X.n_cols;
  int kk = currentGamma.n_rows;
  int n = fy.n_rows;
  mat vary_y = fy * currentC * currentGamma.t();
  mat J(kk, kk, fill::zeros);
  mat S_tau0(kk, 1, fill::zeros);
  mat fC = fy * currentC;
  for(int j=0; j<q; j++){
    List Xindex = extractMatrix(j, grid, X);
    for(int i=0; i < n; ++i){
      mat vary_y_i = vary_y.row(i);
      mat Xj = Xindex[2];
      uvec ind_jj = Xindex[0];
      uvec ind_jjprime = Xindex[1];
      
      mat Xij = Xj.row(i);
      
      double eps_ij = as_scalar(currenttau0.rows(ind_jj)) +
        as_scalar(vary_y_i.cols(ind_jj)) +
        as_scalar(crossprod(currenttau0.rows(ind_jjprime), Xij.t())) +
        as_scalar(crossprod(vary_y_i.cols(ind_jjprime).t(), Xij.t()));
      
      double probs_ij = plogis(eps_ij);
      double vars_ij = probs_ij*(1-probs_ij);
      
      S_tau0.rows(ind_jj) += (X(i, j) - probs_ij);
      S_tau0.rows(ind_jjprime) += (X(i, j) - probs_ij)  *  Xij.t();
      
      mat onesmatrix = ones(1, 1);
      mat XmatrixInvolved = join_horiz(onesmatrix, Xij);
      uvec tmp_index = join_vert(ind_jj, ind_jjprime);
      J(tmp_index, tmp_index) -=  vars_ij * crossprod(XmatrixInvolved, XmatrixInvolved);
    }
  }
  mat newtau0 = currenttau0 - solve(J, S_tau0);
  return(newtau0);
}


// [[Rcpp::export]]
mat globalupdateC(mat X, mat fy, mat currenttau0, mat currentGamma, mat currentC, Mat<int> grid){
  int q = X.n_cols;
  int kk = currentGamma.n_rows;
  int n = fy.n_rows;
  mat vary_y = fy * currentC * currentGamma.t();
  mat J(currentC.n_rows*currentC.n_cols, currentC.n_rows*currentC.n_cols, fill::zeros);
  mat S_C(currentC.n_rows, currentC.n_cols, fill::zeros);
  mat fC = fy * currentC;
  for(int j=0; j<q; j++){
    List Xindex = extractMatrix(j, grid, X);
    for(int i=0; i < n; ++i){
      mat vary_y_i = vary_y.row(i);
      mat Xj = Xindex[2];
      uvec ind_jj = Xindex[0];
      uvec ind_jjprime = Xindex[1];
      
      mat Xij = Xj.row(i);
      
      double eps_ij = as_scalar(currenttau0.rows(ind_jj)) +
        as_scalar(vary_y_i.cols(ind_jj)) +
        as_scalar(crossprod(currenttau0.rows(ind_jjprime), Xij.t())) +
        as_scalar(crossprod(vary_y_i.cols(ind_jjprime).t(), Xij.t()));
      
      double probs_ij = plogis(eps_ij);
      double vars_ij = probs_ij*(1-probs_ij);
      
      mat depsdC = fy.row(i).t() * currentGamma.rows(ind_jj) +
        kron(currentGamma.rows(ind_jjprime), fy.row(i)).t() * Xij.t();
        
      //cout << "number of rows: " << depsdC.n_rows << endl;  
      //cout << "number of columns: " << depsdC.n_cols << endl;
      
      S_C += (X(i, j) - probs_ij) * depsdC;
      
      J -=  vars_ij * kron(depsdC, depsdC.t());
    }
  }
  mat Scvech = reshape(S_C, currentC.n_rows*currentC.n_cols, 1);
  mat newCvech = currentC - solve(J, Scvech);
  mat newC = reshape(newCvech, currentC.n_rows, currentC.n_cols);
  return(newC);
}

