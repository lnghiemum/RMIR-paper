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
double llhOnecluster(mat Xi, mat fyi, mat Wi, rowvec muX, mat Gammai, mat C, 
                             mat beta, rowvec muW, mat Delta){
                             
  int m = Xi.n_rows;  
  int p = Xi.n_cols;
  
  mat muXrep = repmat(muX, m, 1);
  mat term2_part1 = Xi - muXrep;
  mat term2_part2 = fyi * C.t() * Gammai.t();

  mat muWrep = repmat(muW, m, 1); 
  mat tmp2 = Wi - muWrep; 
  mat term3 = tmp2 * beta.t(); 
  mat meanTerm = term2_part1 - term2_part2 - term3;
  vec d = dmvnrm(meanTerm, rowvec(p, fill::zeros), Delta, true);
  double term1 = sum(d);
  return(term1);
 }

// [[Rcpp::export]]
vec computeMarginalLLH(mat X, DataFrame Y, mat fy, mat W, rowvec muX, 
                       mat Gamma0, mat C, rowvec muW, mat beta, mat Delta, List GeneratedGamma, int B){
  
  IntegerVector clusterdf = Y["cluster"]; 
  ivec clusterIndex = as<ivec>(wrap(clusterdf));
  // clusterIndex.print();
  int n = max(clusterdf);
  vec marginalLLh(n, fill::zeros); 
  for(int i=0; i < n; i++){
    uvec clusteriindex = find(clusterIndex == i+1);
    mat Xi = X.rows(clusteriindex);
    mat fyi = fy.rows(clusteriindex);
    mat Wi = W.rows(clusteriindex);
    double term1 = 0;
    for(int b = 0; b < B; b++){
      mat Gammaib =  GeneratedGamma[b];
      term1 += exp(llhOnecluster(Xi, fyi, Wi, muX, Gammaib, C, beta, muW, Delta));
    }
    marginalLLh(i) = term1/B;
  }
  return(log(marginalLLh));    
}

// [[Rcpp::export]]
mat computeWeightsEstep(mat X, DataFrame Y, mat fy, mat W, mat muX, 
                       mat Gamma0, mat C, mat muW, mat beta, mat Delta, List GeneratedGamma, int B){
  
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
    mat Wi = W.rows(clusteriindex);
    for(int b = 0; b < B; b++){
      mat Gammaib =  GeneratedGamma[b];
      wb(i, b) = exp(llhOnecluster(Xi, fyi, Wi, muX, Gammaib, C, beta, muW, Delta));
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
rowvec updatemuX(mat X, DataFrame y, mat fy, mat W, 
                 mat beta, mat C, rowvec muW, mat wb, List GeneratedGamma){
  int N = X.n_rows;
  IntegerVector clusterdf = y["cluster"]; 
  ivec clusterIndex = as<ivec>(wrap(clusterdf));
  int n = max(clusterdf);
  int p = X.n_cols;
  int B = wb.n_cols;
  mat muWrep = repmat(muW, N, 1);
  mat term1 = X - (W-muWrep) * beta.t(); 
  mat Cf = fy * C.t(); 
  rowvec clusterSum(1, p, fill::zeros);
  
  for(int i=0; i< n; i++){
    uvec clusteriindex = find(clusterIndex == i+1);
    mat term1ii = term1.rows(clusteriindex);
    mat Cfii = Cf.rows(clusteriindex);
    for(int b=0; b<B; b++){
      mat Gammab = GeneratedGamma[b];
      mat tmp = term1ii - Cfii * Gammab.t(); // n_i \times p
      // cout << "Values:" << sum(tmp, 0) << endl;
      clusterSum +=  wb(i, b) * sum(tmp, 0);
    }
    // cout << "Values:" << clusterSum << endl;
  }
  rowvec muX = clusterSum/N;
  //cout << "Values:" << muX << endl;
  return(muX);
}
// [[Rcpp::export]]
rowvec updatemuW(mat X, DataFrame y, mat fy, mat W, rowvec muX,
                 mat beta, mat C, mat wb, List GeneratedGamma){
  int N = X.n_rows;
  IntegerVector clusterdf = y["cluster"]; 
  ivec clusterIndex = as<ivec>(wrap(clusterdf));
  int n = max(clusterdf);
  int q = W.n_cols;
  int B = wb.n_cols;
  int p = X.n_cols;
  mat clusterSum(n, p, fill::zeros);
  mat muXrep = repmat(muX, N, 1);
  mat btb = beta.t() * beta;
  mat term1 = X - muXrep - W * beta.t(); // N \times p
  
  mat Cf = fy * C.t(); // # N \times d
  for(int i=0; i< n; i++){
    uvec clusteriindex = find(clusterIndex == i+1);
    mat term1ii = term1.rows(clusteriindex);
    mat Cfii = Cf.rows(clusteriindex);
    for(int b=0; b<B; b++){
      mat Gammab = GeneratedGamma[b];
      mat tmp = term1ii - Cfii * Gammab.t();
      clusterSum.row(i) +=  wb(i, b) * sum(tmp, 0);
    }
  }
  
  vec muW = -solve(btb, beta.t()) * sum(clusterSum, 0).t()/N;  
  rowvec hatmuW = muW.t();
  return(hatmuW);
}
// [[Rcpp::export]]
mat updatebeta(mat X, DataFrame y, mat fy, mat W, 
               rowvec muX, mat C, rowvec muW, mat wb, List GeneratedGamma){
  int N = X.n_rows;
  IntegerVector clusterdf = y["cluster"]; 
  ivec clusterIndex = as<ivec>(wrap(clusterdf));
  int n = max(clusterdf);
  int q = W.n_cols;
  int p = X.n_cols;
  int B = wb.n_cols;
  
  mat muWmat = repmat(muW, N, 1);
  mat muXmat = repmat(muX, N, 1);
  mat Wc = W - muWmat;  
  mat WctWc = Wc.t() * Wc; // q \times q
  mat term1 = X - muXmat;  
  mat term2 = fy * C.t();// N \times d
  
  mat clusterSum(p, q, fill::zeros);
  
  for(int i=0; i< n; i++){
    uvec clusteriindex = find(clusterIndex == i+1);
    mat term1ii = term1.rows(clusteriindex); // n_i \times p
    mat Cfii = term2.rows(clusteriindex); // n_i \times d
    mat Wcii = Wc.rows(clusteriindex); // n_i \times q
    for(int b=0; b<B; b++){
      mat Gammab = GeneratedGamma[b]; // p \times d
      mat tmp = term1ii - Cfii * Gammab.t();
      mat tmp2 = tmp.t() * Wcii; // p \times q  
      clusterSum +=  wb(i, b) * tmp2;
    }
  }
  mat tmp = clusterSum.t(); 
  mat beta = solve(WctWc, tmp).t();
  return(beta);
}
// [[Rcpp::export]]
mat updateDelta(mat X, DataFrame y, mat fy, mat W, 
               rowvec muX, mat C, rowvec muW, mat beta, mat wb, 
               List GeneratedGamma){
  int N = X.n_rows;
  int B = wb.n_cols;
  IntegerVector clusterdf = y["cluster"]; 
  ivec clusterIndex = as<ivec>(wrap(clusterdf));
  int n = max(clusterdf);
  int p = X.n_cols;
  
  mat muWmat = repmat(muW, N, 1);
  mat muXmat = repmat(muX, N, 1);
  mat Wc = W - muWmat;  
  mat term1 = X - muXmat - Wc * beta.t();  
  mat term2 = fy * C.t();// N \times d
  mat clusterSum(p, p, fill::zeros);
  
  for(int i=0; i< n; i++){
    uvec clusteriindex = find(clusterIndex == i+1);
    mat term1ii = term1.rows(clusteriindex); // n_i \times p
    mat term2ii = term2.rows(clusteriindex); // n_i \times d
    mat Wcii = Wc.rows(clusteriindex); // n_i \times q
    for(int b=0; b<B; b++){
      mat Gammab = GeneratedGamma[b]; // p \times d
      mat tmp = term1ii - term2ii * Gammab.t(); // ni \times p;
      clusterSum +=  wb(i, b)* (tmp.t()*tmp);
    }
  }
  mat Delta = clusterSum/N;
  return(Delta);
}
// [[Rcpp::export]]
mat updateC(mat X, DataFrame y, mat fy, mat W, rowvec muX, 
            mat beta, rowvec muW, mat Delta,
            mat wb, List GeneratedGamma){

  int N = X.n_rows;
  int B = wb.n_cols;
  IntegerVector clusterdf = y["cluster"]; 
  ivec clusterIndex = as<ivec>(wrap(clusterdf));
  int n = max(clusterdf);
  int q = W.n_cols;
  int p = X.n_cols;
  int r = fy.n_cols;
  //cout << "Number of columns:" << r << endl;
  mat Gammab = GeneratedGamma[1];
  int d = Gammab.n_cols;
  mat muWmat = repmat(muW, N, 1);
  mat muXmat = repmat(muX, N, 1);
  mat Wc = W - muWmat;  
  mat term1 = X - muXmat - Wc * beta.t();  
  
  mat lhs(r*d, r*d, fill::zeros);
  mat rhs(r*d, 1, fill::zeros);
  
  for(int i=0; i< n; i++){
    uvec clusteriindex = find(clusterIndex == i+1);
    mat term1ii = term1.rows(clusteriindex); // n_i \times p
    mat fyii = fy.rows(clusteriindex); // n_i \times r
    mat Fyii = fyii.t() * fyii; // r \times r
    mat Gii = fyii.t() * term1ii; // r \times p
    // cout << "Number of columns:" << fyii.n_cols << endl;
    //cout << "Number of rows:" << fyii.n_rows << endl;
    
    for(int b=0; b<B; b++){
      mat Gammab = GeneratedGamma[b]; // p \times d
      mat tmp = Gammab.t() * solve(Delta, Gammab); // d \times d
      mat tmp2 = kron(Fyii, tmp); // rd \times rd
      lhs +=  wb(i, b)*tmp2;
      mat tmp3 = Gammab.t() * solve(Delta, Gii.t());
      mat tmp4 = reshape(tmp3, r*d, 1);
      rhs += wb(i, b)*tmp4;
    }
  }
  mat vecC = solve(lhs, rhs);
  mat hatC = reshape(vecC, d, r);
  return(hatC);
}


  