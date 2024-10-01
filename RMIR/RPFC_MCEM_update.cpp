#include <RcppArmadillo.h>
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]
using namespace arma;

static double const log2pi = std::log(2.0 * M_PI);

void inplace_tri_mat_mult(arma::rowvec &x, arma::mat const &trimat){
  arma::uword const n = trimat.n_cols;
  
  for(unsigned j = n; j-- > 0;){
    double tmp(0.);
    for(unsigned i = 0; i <= j; ++i)
      tmp += trimat.at(i, j) * x[i];
    x[j] = tmp;
  }
}


arma::vec dmvnrm(arma::mat const &x,  
                           arma::rowvec const &mean,  
                           arma::mat const &sigma, 
                           bool const logd = false) { 
    using arma::uword;
    uword const n = x.n_rows, 
             xdim = x.n_cols;
    arma::vec out(n);
    arma::mat const rooti = arma::inv(trimatu(arma::chol(sigma)));
    double const rootisum = arma::sum(log(rooti.diag())), 
                constants = -(double)xdim/2.0 * log2pi, 
              other_terms = rootisum + constants;
    
    arma::rowvec z;
    for (uword i = 0; i < n; i++) {
        z = (x.row(i) - mean);
        inplace_tri_mat_mult(z, rooti);
        out(i) = other_terms - 0.5 * arma::dot(z, z);     
    }  
      
    if (logd)
      return out;
    return exp(out);
}
// [[Rcpp::export]]
double log_dmatnorm(mat Zi, mat M, mat U, mat V){
  int n = Zi.n_rows;
  int p = Zi.n_cols;
  
  double logdet_U = log_det_sympd(U);
  double logdet_V = log_det_sympd(V);

  arma::mat diff = Zi - M;
  arma::mat U_inv = arma::inv(U);
  arma::mat V_inv = arma::inv(V);
  double trace_term = arma::trace(U_inv * diff * V_inv * diff.t());
  
  double log_dens_value = (-n * p / 2.0) * std::log(2 * M_PI) 
    - (p / 2.0) * logdet_U
  - (n / 2.0) * logdet_V
  - (1 / 2.0) * trace_term;
  return(log_dens_value);
} 
// [[Rcpp::export]]
Rcpp::List convert_to_Z(mat Xi, mat fyi){ 
  // Step 1: Get number of rows and columns from Xi
  int m = Xi.n_rows;
  int p = Xi.n_cols;
  
  // Step 2: Calculate the column means (mui)
  arma::rowvec mui = arma::mean(Xi, 0);  // equivalent to colMeans in R
  
  // Step 3: Subtract mui from each column (sweep operation)
  arma::mat Zi_temp = Xi.each_row() - mui;  // Xi - mui across rows
  
  // Step 4: Remove the last row to avoid multicollinearity
  int ni = m - 1;  // ni is m - 1
  arma::mat Zi = Zi_temp.rows(0, m - 2);
  
  // Step 5: Create matrix Ai (correlation matrix among rows of Z)
  arma::mat Ai = arma::eye(ni, ni) - 1.0 / (ni + 1.0);  // Identity matrix of size ni, with adjustment
  
  // Step 6: Center fyi and remove the last row, equivalent to scale(fyi, scale=FALSE, center=TRUE)
  arma::rowvec mean_fyi = arma::mean(fyi, 0);  // Calculate the mean of fyi
  arma::mat fyi_centered = fyi.each_row() - mean_fyi;  
  
  // Step 7: Remove the last row and transpose
  arma::mat fyic = fyi_centered.rows(0, m - 2);
  
  List L = List::create(
    Named("Zi") = Zi , _["fyic"] = fyic, _["Ai"] = Ai);
  return(L);
}  

                   
// [[Rcpp::export]]
double llhOnecluster(mat Xi, mat fyi, mat Gammai, mat C, 
                             mat Delta){
  
  Rcpp::List tmp = convert_to_Z(Xi, fyi);
  mat Ai = tmp["Ai"];
  mat fyic = tmp["fyic"];
  mat Zi = tmp["Zi"];
  
  mat Ctrans = C.t();
  mat Gammatrans = Gammai.t();
  mat M = fyic * Ctrans * Gammatrans;
  
  // cout << "Number of rows for M:" << M.n_rows << endl;
  // cout << "Number of columns for M:" << M.n_cols << endl;
  // cout << "Number of rows for Zi:" << Zi.n_rows << endl;
  // cout << "Number of columns for Zi:" << Zi.n_cols << endl;
  // 
  
  
  double llh = log_dmatnorm(Zi, M, Ai, Delta);
    
  return(llh);
 }

// [[Rcpp::export]]
vec computeMarginalLLH(mat X, DataFrame Y, mat fy, 
                       mat Gamma0, mat C, mat Delta, List GeneratedGamma, int B){
  
  IntegerVector clusterdf = Y["cluster"]; 
  ivec clusterIndex = as<ivec>(wrap(clusterdf));
  // clusterIndex.print();
  int n = max(clusterdf);
  vec marginalLLh(n, fill::zeros); 
  for(int i=0; i < n; i++){
    uvec clusteriindex = find(clusterIndex == i+1);
    mat Xi = X.rows(clusteriindex);
    mat fyi = fy.rows(clusteriindex);
    double term1 = 0;
    for(int b = 0; b < B; b++){
      mat Gammaib =  GeneratedGamma[b];
      term1 += exp(llhOnecluster(Xi, fyi, Gammaib, C, Delta));
    }
    marginalLLh(i) = term1/B;
  }
  return(log(marginalLLh));    
}

// [[Rcpp::export]]
mat computeWeightsEstep(mat X, DataFrame Y, mat fy, 
                       mat Gamma0, mat C, mat Delta, List GeneratedGamma, int B){
  
  IntegerVector clusterdf = Y["cluster"]; 
  ivec clusterIndex = as<ivec>(wrap(clusterdf));
  // clusterIndex.print();
  int n = max(clusterdf);
  mat wb(n, B, fill::zeros);
  vec rowSum(n, fill::zeros);
  for(int i=0; i < n; i++){
    uvec clusteriindex = find(clusterIndex == i+1);
    mat Xi = X.rows(clusteriindex);
    mat fyi = fy.rows(clusteriindex);
    for(int b = 0; b < B; b++){
      mat Gammaib =  GeneratedGamma[b];
      wb(i, b) = exp(llhOnecluster(Xi, fyi, Gammaib, C, Delta));
    }
    rowSum(i) = sum(wb.row(i));
    if (rowSum(i) !=0){
      wb.row(i) = wb.row(i)/rowSum(i);
    } else{
      wb(i, 0) = 1;
    }
  }
  return(wb);    
}
// [[Rcpp::export]]
mat updateDelta(mat X, DataFrame y, mat fy, mat C, mat wb, 
               List GeneratedGamma){
  int N = X.n_rows;
  int B = wb.n_cols;
  IntegerVector clusterdf = y["cluster"]; 
  ivec clusterIndex = as<ivec>(wrap(clusterdf));
  int n = max(clusterdf);
  int p = X.n_cols;
  mat clusterSum(p, p, fill::zeros);
  
  for(int i=0; i< n; i++){
    uvec clusteriindex = find(clusterIndex == i+1);
    mat Xi = X.rows(clusteriindex); // n_i \times p
    int mi = Xi.n_rows-1;
    mat fyi = fy.rows(clusteriindex); // n_i \times d
    Rcpp::List tmp = convert_to_Z(Xi, fyi);
    mat fyic = tmp["fyic"];
    mat Zi = tmp["Zi"];
    mat Ai_inv = arma::eye(mi, mi) + arma::ones(mi, mi);
    
    for(int b=0; b<B; b++){
      mat Gammab = GeneratedGamma[b]; // p \times d
      mat tmp = Zi - fyic * C.t() * Gammab.t(); // ni \times p;
      clusterSum +=  wb(i, b)* (tmp.t() * Ai_inv *tmp);
    }
  }
  mat Delta = clusterSum/N;
  return(Delta);
}

// [[Rcpp::export]]
mat updateC(mat X, DataFrame y, mat fy, mat Delta,
            mat wb, List GeneratedGamma){

  int N = X.n_rows;
  int B = wb.n_cols;
  IntegerVector clusterdf = y["cluster"]; 
  ivec clusterIndex = as<ivec>(wrap(clusterdf));
  int n = max(clusterdf);
  int p = X.n_cols;
  int r = fy.n_cols;
  //cout << "Number of columns:" << r << endl;
  mat Gammab = GeneratedGamma[1];
  int d = Gammab.n_cols;
  
  mat lhs(r*d, r*d, fill::zeros);
  mat rhs(r*d, 1, fill::zeros);
  
  for(int i=0; i< n; i++){
    uvec clusteriindex = find(clusterIndex == i+1);
    mat Xi = X.rows(clusteriindex);
    int mi = Xi.n_rows -1;
    mat fyi = fy.rows(clusteriindex); // n_i \times d
    Rcpp::List tmp = convert_to_Z(Xi, fyi);
    mat fyic = tmp["fyic"];
    mat Zi = tmp["Zi"];
    mat Ai_inv = arma::eye(mi, mi) + arma::ones(mi, mi);
    // cout << "Number of columns:" << fyii.n_cols << endl;
    //cout << "Number of rows:" << fyii.n_rows << endl;
    mat Fyii  = fyic.t() * Ai_inv * fyic; // r \times r
    mat Gii = Zi.t() * Ai_inv * fyic; // p \times r
      
    for(int b=0; b<B; b++){
      mat Gammab = GeneratedGamma[b]; // p \times d
      mat tmp = Gammab.t() * solve(Delta, Gammab); // d \times d
      mat tmp2 = kron(Fyii, tmp); // rd \times rd
      lhs +=  wb(i, b)*tmp2;
      mat tmp3 = Gammab.t() * solve(Delta, Gii);
      mat tmp4 = reshape(tmp3, r*d, 1);
      rhs += wb(i, b)*tmp4;
    }
  }
  mat vecC = solve(lhs, rhs);
  mat hatC = reshape(vecC, d, r);
  return(hatC);
}


  