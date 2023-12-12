#include <Rcpp.h>
#include <iostream>
#include <string>
#include "spams.h"

// [[Rcpp::export]]
void proximalGraph(
    Rcpp::NumericVector& U,
    int& p,
    std::string& regul,
    Rcpp::IntegerMatrix& grp,
    Rcpp::IntegerMatrix& grpV,
    Rcpp::NumericVector& etaG,
    double lam) {
  
  // read in U and convert to spams::Matrix<double> alpha0
  Matrix<double> alpha0(p, 1);
  for (int r = 0; r < p; r++) {
    alpha0(r, 0) = U(r);
  }
  
  // grp dimensions
  int gr = grp.nrow();
  int gc = grp.ncol();
  // read in grp and convert to spams::Matrix<bool> grp_dense
  // then to spams::SpMatrix<bool> groups
  Matrix<bool> grp_dense(gr, gc);
  for (int r = 0; r < gr; r++) {
    for (int c = 0; c < gc; c++) {
      grp_dense(r, c) = (grp(r, c) != 0);
    }
  }
  SpMatrix<bool> groups;
  grp_dense.toSparse(groups);
  
  // grpV dimensions
  int gvr = grpV.nrow();
  int gvc = grpV.ncol();
  // read in grpV and convert to spams::Matrix<bool> grpV_dense
  // then to spams::SpMatrix<bool> groups_var
  Matrix<bool> grpV_dense(gvr, gvc);
  for (int r = 0; r < gvr; r++) {
    for (int c = 0; c < gvc; c++) {
      grpV_dense(r, c) = (grpV(r, c) != 0);
    }
  }
  SpMatrix<bool> groups_var;
  grpV_dense.toSparse(groups_var);
  
  // read in etaG and convert to spams::Vector<double> eta_g
  int n_etaG = etaG.length();
  Vector<double> eta_g(&etaG[0], n_etaG);
  
  // Initialize alpha - proximal operator
  Matrix<double> alpha(p, 1);
  alpha.setZeros();
  
  // call _proximalGraph
  _proximalGraph(&alpha0, &alpha,
                 &eta_g, &groups, &groups_var,
                 1, lam, 0.0, 0.0,
                 false, false, &regul[0],
                                     false, false, true, true,
                                     1, false);
  
  // put updated alpha back into U
  for (int r = 0; r < p; r++) {
    U(r) = alpha(r, 0);
  }
}

// [[Rcpp::export]]
void proximalTree(
    Rcpp::NumericVector& U,
    int& p,
    std::string& regul,
    Rcpp::IntegerMatrix& grp,
    Rcpp::IntegerVector& own_var,
    Rcpp::IntegerVector& N_own_var,
    Rcpp::NumericVector& etaG,
    double lam) {
    
  // read in U and convert to spams::Matrix<double> alpha0
  Matrix<double> alpha0(p, 1);
  for (int r = 0; r < p; r++) {
    alpha0(r, 0) = U(r);
  }
  
  // grp dimensions
  int gr = grp.nrow();
  int gc = grp.ncol();
  // read in grp and convert to spams::Matrix<bool> grp_dense
  // then to spams::SpMatrix<bool> groups
  Matrix<bool> grp_dense(gr, gc);
  for (int r = 0; r < gr; r++) {
    for (int c = 0; c < gc; c++) {
      grp_dense(r, c) = (grp(r, c) != 0);
    }
  }
  SpMatrix<bool> groups;
  grp_dense.toSparse(groups);

  // read in own_var and convert to spams::Vector<int> ov
  int n_ov = own_var.length();
  Vector<int> ov(&own_var[0], n_ov);

  // read in N_own_var and convert to spams::Vector<int> nov
  int n_nov = N_own_var.length();
  Vector<int> nov(&N_own_var[0], n_nov);
  
  // read in etaG and convert to spams::Vector<int> eta_g
  int n_etaG = etaG.length();
  Vector<double> eta_g(&etaG[0], n_etaG);
  
  // Initialize alpha - proximal operator
  Matrix<double> alpha(p, 1);
  alpha.setZeros();
  
  // call _proximalTree
  _proximalTree(&alpha0, &alpha,
                &eta_g, &groups,
                &ov, &nov,
                1, lam, 0.0, 0.0,
                false, false, &regul[0],
                false, false, true, true,
                1, false);
  
  // put updated alpha back into U
  for (int r = 0; r < p; r++) {
    U(r) = alpha(r, 0);
  }
}

// [[Rcpp::export]]
Rcpp::List sox_cpp(Rcpp::NumericMatrix& x,
                   Rcpp::NumericVector& start,
                   Rcpp::NumericVector& stop,
                   Rcpp::NumericVector& event,
                   int& n_unique,
                   std::string& regul,
                   Rcpp::NumericVector& lam,
                   Rcpp::IntegerMatrix& grp,
                   Rcpp::IntegerMatrix& grpV,
                   Rcpp::IntegerVector& own_var,
                   Rcpp::IntegerVector& N_own_var,
                   Rcpp::NumericVector& etaG,
                   Rcpp::NumericVector& init,
                   Rcpp::Function& l_ld,
                   double& init_stepsize,
                   double& ls_shrink,
                   double& partol,
                   int& maxit,
                   bool& verbose) {
  
  
  int p;
  p = x.ncol();
  int i;
  int nlam;
  nlam = lam.length();
  double lam_tmp;
  int l;
  
  Rcpp::NumericVector beta = clone(init);
  Rcpp::NumericVector lkhd_grad (p + 1);
  double lkhd;
  Rcpp::NumericVector grad (p);
  
  Rcpp::NumericVector beta_tmp (p);
  Rcpp::NumericVector lkhd_grad_tmp (p + 1);
  double lkhd_tmp;
  Rcpp::NumericVector grad_tmp (p);
  
  double stepsize;
  stepsize = init_stepsize;
  bool ls_fail;
  Rcpp::NumericVector par_diff (p);
  double par_diff_norm;
  
  Rcpp::NumericMatrix estimates (p, nlam);
  Rcpp::NumericVector iterations (nlam);
  
  for (l = 0; l < nlam; l++) {
    lam_tmp = lam[l];
    i = 0;
    
    while (i < maxit) {
      lkhd_grad = l_ld(beta, x, n_unique, start, stop, event);
      lkhd = lkhd_grad[0];
      grad = lkhd_grad[Rcpp::Range(1, p)];
      
      beta_tmp = beta - grad * stepsize;
      if (regul == "graph") {
        proximalGraph(beta_tmp, p, regul, grp, grpV, etaG, stepsize * lam_tmp);
      } else if (regul == "tree-l2") {
        proximalTree(beta_tmp, p, regul, grp, own_var, N_own_var, etaG, stepsize * lam_tmp);
      }
      
      par_diff = beta_tmp - beta;
      
      lkhd_grad_tmp = l_ld(beta_tmp, x, n_unique, start, stop, event);
      lkhd_tmp = lkhd_grad_tmp[0];
      // grad_tmp = lkhd_grad_tmp[Rcpp::Range(1, p)];
      
      par_diff_norm = sum(pow(par_diff, 2));
      ls_fail = ( lkhd_tmp > (lkhd + sum(grad * par_diff) + par_diff_norm/2/stepsize) );
      
      while (ls_fail) {
        stepsize = ls_shrink * stepsize;
        
        beta_tmp = beta - grad * stepsize;
        if (regul == "graph") {
          proximalGraph(beta_tmp, p, regul, grp, grpV, etaG, stepsize * lam_tmp);
        } else if (regul == "tree-l2") {
          proximalTree(beta_tmp, p, regul, grp, own_var, N_own_var, etaG, stepsize * lam_tmp);
        }
        
        par_diff = beta_tmp - beta;
        
        lkhd_grad_tmp = l_ld(beta_tmp, x, n_unique, start, stop, event);
        lkhd_tmp = lkhd_grad_tmp[0];
        // grad_tmp = lkhd_grad_tmp[Rcpp::Range(1, p)];
        
        par_diff_norm = sum(pow(par_diff, 2));
        ls_fail = ( lkhd_tmp > (lkhd + sum(grad * par_diff) + par_diff_norm/2/stepsize) );
      }
      
      beta = clone(beta_tmp);
      
      i++;

      if (verbose) {
        Rcpp::Rcout << "Lambda = " << lam_tmp << "\n";
        Rcpp::Rcout << "Iteration " << i;
        Rcpp::Rcout << "; Stepsize = " << stepsize << "\n";
        Rcpp::Rcout << "||beta - beta_prev||/p = " << par_diff_norm/p << "\n\n";
      }
      
      if ( par_diff_norm/p < partol ) break;
    }
    estimates.column(l) = beta;
    iterations[l] = i;
  }
  
  Rcpp::List result = Rcpp::List::create(Rcpp::Named("Estimates") = estimates,
                                         Rcpp::Named("Iterations") = iterations);
  return result;
}


