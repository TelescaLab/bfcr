#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#include <cmath>
#include "Utility.h"
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]


// [[Rcpp::export]]
void updateProjBeta(arma::cube& Lambda, arma::mat& Theta, arma::mat& Eta, arma::vec& Delta,
                double Prec, arma::mat& X, arma::mat& Y, arma::mat B, arma::mat& Proj, double beta){
  arma::uword p = B.n_cols;
  for(arma::uword i = 0; i < Y.n_rows; i++){
    arma::vec b = Prec * B.t() * Y.row(i).t() + arma::diagmat(Delta) * Theta * X.row(i).t();
    for(arma::uword k = 0; k < Lambda.n_slices; k++){
      b = b + arma::diagmat(Delta) * Lambda.slice(k) * X.row(i).t() * Eta(i, k);
    }
    arma::mat C = (Prec * B.t() * B + arma::diagmat(Delta));
    arma::mat mychol = arma::chol(C, "lower");
    arma::vec w = solve(arma::trimatl(mychol), b);
    arma::vec mu = solve(arma::trimatu(mychol.t()), w);
    arma::vec z = 1 / sqrt(beta) * arma::randn<arma::vec>(p);
    arma::vec v = arma::solve(arma::trimatu(mychol.t()), z);
    Proj.row(i) = (mu + v).t();
  }
}

// [[Rcpp::export]]
double updatePrecPBeta(arma::mat& Proj, arma::mat& Y, arma::mat& B, double beta){
  double a = 0;
  double b = 0;
  //double my_sum = arma::sum(arma::square(Y.t() - B * Proj.t()));
  double my_sum = arma::accu(arma::square(Y.t() - B * Proj.t()));
  return(R::rgamma(beta * (a + Y.n_elem/2) + 1 - beta, 1.0 / (beta * (b + 1.0/2.0 * my_sum))));
}


// [[Rcpp::export]]
void updateTauBeta(arma::mat& Theta, arma::cube& Lambda, arma::mat& Tau, double beta){
  double t_alpha = 1;
  double t_beta =  .000000005;
  //double t_beta = 1;
  arma::mat P = getPenalty2(Lambda.n_rows, 2);
  arma::uword R = Lambda.n_slices;
  arma::uword D = Lambda.n_cols;
  for(arma::uword d = 0; d < D; d++){
    Tau(0, d) = R::rgamma(beta * (t_alpha + (Lambda.n_rows - 1.0) / 2.0) + 1 - beta, 1.0 / (beta * (t_beta + 1.0 / 2.0 * arma::as_scalar(Theta.col(d).t() * P * Theta.col(d)))));
  }
  /*
  for(arma::uword r = 0; r < R; r++){
    for(arma::uword d = 0; d < D; d++){
      Tau(r + 1, d) = R::rgamma(beta * (t_alpha + (Lambda.n_rows - 1.0) / 2.0) + 1 - beta, 1.0 / (beta * (t_beta + 1.0 / 2.0 * arma::as_scalar(Lambda.slice(r).col(d).t() * P * Lambda.slice(r).col(d)))));
    }
  }
  */
  double temp_beta;
  for(arma::uword d = 0; d < D; d++){
    temp_beta = 0;
    for(arma::uword r = 1; r < R; r++){
      temp_beta = temp_beta + 1.0 / 2.0 * arma::as_scalar(Lambda.slice(r).col(d).t() * P * Lambda.slice(r).col(d));
    }
    Tau(1, d) = R::rgamma(beta * (t_alpha + R * (Lambda.n_rows - 1.0) / 2.0) + 1 - beta, 1.0 / (beta * (t_beta + temp_beta)));
    for(arma::uword r = 0; r < R; r++){
      Tau(r + 1, d) = Tau(1, d);
    }
  }
}

/*
// [[Rcpp::export]]
double updateProjT(arma::cube& Lambda, arma::mat& Theta, arma::mat& Eta, arma::vec& Delta,
                 double Prec, arma::mat& X, arma::mat& Y, arma::mat B, arma::mat& Proj, double beta){
  arma::uword p = B.n_cols;
  double u = R::runif(0.0, 1.0);
  double f = 0;
  double g = 0;
  arma::mat C = Prec * B.t() * B + arma::diagmat(Delta);
  arma::mat C_beta = beta * C;
  arma::mat C_sample;
  double val, val_beta, sign;
  arma::log_det(val_beta, sign, C_beta);
  arma::log_det(val, sign, C);
  arma::vec v;
  for(arma::uword i = 0; i < Y.n_rows; i++){
    arma::vec b = Prec * B.t() * Y.row(i).t() + arma::diagmat(Delta) * Theta * X.row(i).t();
    for(arma::uword k = 0; k < Lambda.n_slices; k++){
      b = b + arma::diagmat(Delta) * Lambda.slice(k) * X.row(i).t() * Eta(i, k);
    }
    
    
    arma::mat mychol = arma::chol(C, "lower");
    arma::vec w = solve(arma::trimatl(mychol), b);
    arma::vec mu = solve(arma::trimatu(mychol.t()), w);
    arma::vec z = arma::randn<arma::vec>(p);
    if(u < 0.5){
      v = 1/sqrt(beta) * arma::solve(arma::trimatu(mychol.t()), z);
    } else{
      v = arma::solve(arma::trimatu(mychol.t()), z);
    }
    
    Proj.row(i) = (mu + v).t();
    if(u < 0.5){
      f = (1.0/2.0 * val) - 1.0/2.0 *
                  arma::as_scalar((Proj.row(i).t() - mu).t() * C * (Proj.row(i).t() - mu)) + f;
      g = (1.0/2.0 * val_beta) -  1.0/2.0 *
                  arma::as_scalar((Proj.row(i).t() - mu).t() * C_beta * (Proj.row(i).t() - mu)) + g;
    } 
  }
  return(.5 + .5*exp(g-f));
}
*/

// [[Rcpp::export]]
arma::uword choose_coordinate(arma::vec log_weights){
  arma::uvec sequence = arma::linspace<arma::uvec>(1, log_weights.n_elem, log_weights.n_elem);
  double max = log_weights.max();
  arma::vec normalized = arma::exp(log_weights - max);
  normalized = normalized / arma::sum(normalized);
  arma::uvec out = Rcpp::RcppArmadillo::sample(sequence, 1, false, normalized);
  return(out(0));
}

// [[Rcpp::export]]
double updateThetaLambdaPT(arma::cube& Lambda, arma::mat& Theta, arma::mat& Eta,
                         arma::vec& Delta, arma::mat& Proj, arma::mat& Tau, arma::mat& X,
                         double beta){
  arma::uword d = X.n_cols;
  arma::uword p = Theta.n_rows;
  arma::uword K = Lambda.n_slices;
  arma::uword n = Proj.n_rows;
  arma::mat P = getPenalty2(p, 1);
  arma::mat X_eta(n, (1 + K) * d);
  double val, val_beta, sign;
  arma::mat C, C_beta, mychol, mychol_beta;
  arma::vec b, b_beta, w, mu, z, v, w_beta, mu_beta, z_beta, v_beta, result;
  double pr = 0;
  double u = R::runif(0.0, 1.0);
  X_eta.cols(0, d - 1) = X;
  for(arma::uword k = 0; k < K; k++){
    X_eta.cols(d*(k+1), d*(k+2)-1) = arma::diagmat(Eta.col(k)) * X;
  }
  b = arma::vectorise(arma::diagmat(Delta) * Proj.t() * X_eta);
  C = arma::kron( X_eta.t() * X_eta, arma::diagmat(Delta)) + arma::kron(arma::diagmat(arma::vectorise(Tau.t())), P);
  mychol = arma::chol(C, "lower");
  w = solve(arma::trimatl(mychol), b);
  mu = solve(arma::trimatu(mychol.t()), w);
  z = arma::randn<arma::vec>(p*d*(K+1));
  v = arma::solve(arma::trimatu(mychol.t()), z);
  result = mu + v;
  if(u < 0.5){
    
    b_beta = beta * arma::vectorise(arma::diagmat(Delta) * Proj.t() * X_eta);
    C_beta = beta * arma::kron( X_eta.t() * X_eta, arma::diagmat(Delta)) + beta * arma::kron(arma::diagmat(arma::vectorise(Tau.t())), P);
    arma::log_det(val_beta, sign, C_beta);
    arma::log_det(val, sign, C);
    mychol_beta = arma::chol(C_beta, "lower");
    w_beta = solve(arma::trimatl(mychol_beta), b_beta);
    mu_beta = solve(arma::trimatu(mychol_beta.t()), w_beta);
    z_beta = arma::randn<arma::vec>(p*d*(K+1));
    v_beta = arma::solve(arma::trimatu(mychol_beta.t()), z_beta);
    result = mu_beta + v_beta;
    double f = 1.0 / 2.0 * val - 1.0 / 2.0 * arma::as_scalar((result - mu).t() * C * (result - mu));
    double g = 1.0 / 2.0 * val_beta - 1.0 / 2.0 * arma::as_scalar((result - mu_beta).t() * C_beta * (result - mu_beta));
    pr = g - f + log(1 + .5 / exp(g-f));

  }

  Theta = arma::reshape(result.rows(0, p*d-1), p, d);
  for(arma::uword k = 0; k < K; k++){
    Lambda.slice(k) = arma::reshape(result.rows(p*d*(k+1), p*d*(k+2) - 1), p, d);
  }
  return(pr);
}

// [[Rcpp::export]]
double updateEtaPT(arma::cube& Lambda, arma::mat& Theta, arma::mat& Eta,
                 arma::vec& Delta, arma::mat& Proj, arma::mat& X, double beta){
  arma::uword n = Proj.n_rows;
  arma::mat Xtilde(Proj.n_cols, Lambda.n_slices);
  double pr = 0;
  double u = R::runif(0.0, 1.0);
  double val, val_beta, sign;
  arma::mat mychol_beta, C_beta;
  arma::vec b_beta;
  arma::vec w_beta;
  arma::vec mu_beta;
  arma::vec z_beta;
  arma::vec v_beta;
  for(arma::uword i = 0; i < n; i++){
    for(arma::uword k = 0; k < Lambda.n_slices; k++){
      Xtilde.col(k) = Lambda.slice(k) * X.row(i).t();
    }
    
    arma::vec b = Xtilde.t() * arma::diagmat(Delta) * (Proj.row(i).t() - Theta * X.row(i).t());
    arma::mat C = Xtilde.t() * arma::diagmat(Delta) * Xtilde + arma::eye(Lambda.n_slices, Lambda.n_slices);

    arma::mat mychol = arma::chol(C, "lower");
    arma::vec w = solve(arma::trimatl(mychol), b);
    arma::vec mu = solve(arma::trimatu(mychol.t()), w);
    arma::vec z = arma::randn<arma::vec>(Lambda.n_slices);
    arma::vec v = arma::solve(arma::trimatu(mychol.t()), z);
    Eta.row(i) = (mu + v).t();
    
    if(u < .5){
      arma::vec b_beta = Xtilde.t() * arma::diagmat(Delta) * (Proj.row(i).t() - Theta * X.row(i).t());
      arma::mat C_beta = beta * Xtilde.t() * arma::diagmat(Delta) * Xtilde + beta * arma::eye(Lambda.n_slices, Lambda.n_slices);
      
      arma::mat mychol_beta = arma::chol(C_beta, "lower");
      arma::vec w_beta = solve(arma::trimatl(mychol_beta), b_beta);
      arma::vec mu_beta = solve(arma::trimatu(mychol_beta.t()), w_beta);
      arma::vec z_beta = arma::randn<arma::vec>(Lambda.n_slices);
      arma::vec v_beta = arma::solve(arma::trimatu(mychol_beta.t()), z_beta);
      Eta.row(i) = (mu_beta + v_beta).t();
    
      arma::log_det(val, sign, C);
      arma::log_det(val_beta, sign, C_beta);
      double f = 1.0 / 2.0 * val
      - 1.0 / 2.0 * arma::as_scalar((Eta.row(i).t() - mu).t() * C * (Eta.row(i).t() - mu));
      double g = 1.0 / 2.0 * val_beta
        - 1.0 / 2.0 * arma::as_scalar((Eta.row(i).t() - mu_beta).t() * C_beta * (Eta.row(i).t() - mu_beta));
      pr = pr + (g - f) + log(1 + .5 / exp(g-f));
    } 
  }
  return(pr);
  
}

// [[Rcpp::export]]
void updateProj(arma::cube& Lambda, arma::mat& Theta, arma::mat& Eta, arma::vec& Delta, double Prec, arma::mat& X, arma::mat& Y, arma::mat B, arma::mat& Proj){
  arma::uword p = B.n_cols;
  for(arma::uword i = 0; i < Y.n_rows; i++){
    arma::vec b = Prec * B.t() * Y.row(i).t() + arma::diagmat(Delta) * Theta * X.row(i).t();
    for(arma::uword k = 0; k < Lambda.n_slices; k++){
      b = b + arma::diagmat(Delta) * Lambda.slice(k) * X.row(i).t() * Eta(i, k);
    }
    arma::mat C = Prec * B.t() * B + arma::diagmat(Delta);
    arma::mat mychol = arma::chol(C, "lower");
    arma::vec w = solve(arma::trimatl(mychol), b);
    arma::vec mu = solve(arma::trimatu(mychol.t()), w);
    arma::vec z = arma::randn<arma::vec>(p);
    arma::vec v = arma::solve(arma::trimatu(mychol.t()), z);
    Proj.row(i) = (mu + v).t();
  }
}
// [[Rcpp::export]]
void updateThetaLambdaP(arma::cube& Lambda, arma::mat& Theta, arma::mat& Eta, arma::vec& Delta, arma::mat& Proj, arma::mat& Tau, arma::mat& X){
  arma::uword d = X.n_cols;
  arma::uword p = Theta.n_rows;
  arma::uword K = Lambda.n_slices;
  arma::uword n = Proj.n_rows;
  arma::mat P = getPenalty2(p, 2);
  arma::mat X_eta(n, (1 + K) * d);
  X_eta.cols(0, d - 1) = X;
  for(arma::uword k = 0; k < K; k++){
    X_eta.cols(d*(k+1), d*(k+2)-1) = arma::diagmat(Eta.col(k)) * X;
  }
  arma::vec b = arma::vectorise(arma::diagmat(Delta) * Proj.t() * X_eta);
  arma::mat C = arma::kron( X_eta.t() * X_eta, arma::diagmat(Delta)) + arma::kron(arma::diagmat(arma::vectorise(Tau.t())), P);
  arma::mat mychol = arma::chol(C, "lower");
  arma::vec w = solve(arma::trimatl(mychol), b);
  arma::vec mu = solve(arma::trimatu(mychol.t()), w);
  arma::vec z = arma::randn<arma::vec>(p*d*(K+1));
  arma::vec v = arma::solve(arma::trimatu(mychol.t()), z);

  arma::vec result = mu + v;

  Theta = arma::reshape(result.rows(0, p*d-1), p, d);
  for(arma::uword k = 0; k < K; k++){
    Lambda.slice(k) = arma::reshape(result.rows(p*d*(k+1), p*d*(k+2) - 1), p, d);
  }
}

// [[Rcpp::export]]
void updateEtaP(arma::cube& Lambda, arma::mat& Theta, arma::mat& Eta, arma::vec& Delta, arma::mat& Proj, arma::mat& X){
  arma::uword n = Proj.n_rows;
  arma::mat Xtilde(Proj.n_cols, Lambda.n_slices);
  for(arma::uword i = 0; i < n; i++){
    for(arma::uword k = 0; k < Lambda.n_slices; k++){
      Xtilde.col(k) = Lambda.slice(k) * X.row(i).t();
    }
    arma::vec b = Xtilde.t() * arma::diagmat(Delta) * (Proj.row(i).t() - Theta * X.row(i).t());
    arma::mat C = Xtilde.t() * arma::diagmat(Delta) * Xtilde + arma::eye(Lambda.n_slices, Lambda.n_slices);
    arma::mat mychol = arma::chol(C, "lower");
    arma::vec w = solve(arma::trimatl(mychol), b);
    arma::vec mu = solve(arma::trimatu(mychol.t()), w);
    arma::vec z = arma::randn<arma::vec>(Lambda.n_slices);
    arma::vec v = arma::solve(arma::trimatu(mychol.t()), z);
    Eta.row(i) = (mu + v).t();
  }
}

// [[Rcpp::export]]
void updateDelta(arma::mat& Proj, arma::mat& Theta, arma::cube& Lambda, arma::mat& Eta, arma::vec& Delta, arma::mat& X){
  arma::uword p = Proj.n_cols;
  arma::uword n = Proj.n_rows;
  arma::uword K = Lambda.n_slices;
  //double a = .1;
  //double b = .1;
  double a = 1;
  double b = .1;
  arma::vec my_sum = arma::zeros<arma::vec>(p);
  arma::vec ptm;
  for(arma::uword i = 0; i < n; i++){
    ptm = Theta * X.row(i).t();
    for(arma::uword k = 0; k < K; k++){
      ptm = ptm + Lambda.slice(k) * X.row(i).t() * Eta(i, k); 
    }
    my_sum = my_sum + arma::square(Proj.row(i).t() - ptm);
  }
  for(arma::uword lp = 0; lp < p; lp++){
    Delta(lp) = R::rgamma(a + n/2, 1.0 / (b + 1.0/2.0 * my_sum(lp)));
  }
}

// [[Rcpp::export]]
double updatePrecP(arma::mat& Proj, arma::mat& Y, arma::mat& B){
  double a = 0;
  double b = 0;
  //double my_sum = arma::sum(arma::square(Y.t() - B * Proj.t()));
  double my_sum = arma::accu(arma::square(Y.t() - B * Proj.t()));
  return(R::rgamma(a + Y.n_elem/2, 1.0 / (b + 1.0/2.0 * my_sum)));
}
// [[Rcpp::export]]
void updateLambda(arma::mat& Y, arma::cube& Lambda, arma::vec& r, arma::mat& Gamma, arma::mat& X, arma::mat& B, double prec){
  arma::mat G;
  arma::mat condMean = arma::zeros<arma::mat>(Y.n_rows, Y.n_cols);
  arma::uword Q = Lambda.n_slices;
  arma::uword counter = 0;
  arma::mat mu;
  arma::mat Z;
  arma::mat D;
  arma::mat Pb;
  
  //arma::mat Dchol;
  //arma::mat Pbchol;
  for(arma::uword outer = 0; outer < Q; outer++){
    for(arma::uword q = 0; q < Q; q++){
      if(q != counter){
        condMean = condMean + arma::diagmat(Gamma.col(q)) * X * arma::trans(Lambda.slice(q)) * B.t();
      }
    }
    G = arma::diagmat(Gamma.col(outer));
    mu = arma::trans(prec * arma::solve(X.t() * arma::square(arma::diagmat(G)) * X, X.t()) * 
      arma::diagmat(G) * (Y - condMean) * B * arma::inv(prec * B.t() * B + r(outer) * getPenalty2(B.n_cols, 2)));
    Z = arma::randn<arma::mat>(B.n_cols, X.n_cols);
    D = (arma::chol(X.t() * arma::square(arma::diagmat(G)) * X));
    Pb = (arma::chol(prec * B.t() * B + r(outer) * getPenalty2(B.n_cols, 2), "lower"));
    //Dchol = arma::chol(D);
    //Pbchol = arma::chol(Pb, "lower");
    Lambda.slice(outer) = mu + arma::solve(arma::trimatl(Pb), arma::solve(arma::trimatl(D.t()), Z.t()).t());
    condMean.zeros();
    counter++;
  }
}

// [[Rcpp::export]]
void updateLambda2(arma::mat& Y, arma::cube& Lambda, arma::mat& Tau, arma::mat& Gamma, arma::mat& X, arma::mat& B, double prec, arma::mat& Theta){
  //arma::mat P = getPenalty(Lambda.n_rows);
  arma::mat P = getPenalty2(Lambda.n_rows, 2);
  arma::mat PriorPrecision;
  arma::mat Precision;
  arma::mat mychol;
  arma::vec mymean;
  arma::vec Ltemp;
  arma::uword Q = Lambda.n_slices;
  arma::uword counter = 0;
  arma::mat condMean = arma::zeros<arma::mat>(Y.n_rows, Y.n_cols);
  for(arma::uword i = 0; i < Q; i++){
    for(arma::uword q = 0; q < Q; q++){
      if(q != counter){
        condMean = condMean + arma::diagmat(Gamma.col(q)) * X * arma::trans(Lambda.slice(q)) * B.t();
      }
    }
    //PriorPrecision = arma::kron(P, arma::diagmat(Tau.row(i+1))) + arma::kron(arma::eye<arma::mat>(Lambda.n_rows, Lambda.n_rows), arma::diagmat(c.row(i)));
    PriorPrecision = arma::kron(P, arma::diagmat(Tau.row(i+1)));
    Precision = prec * arma::kron(B.t() * B, X.t() * arma::diagmat(arma::square(Gamma.col(i))) * X) + PriorPrecision;
    mychol = arma::chol(Precision, "lower");
    mymean = arma::solve(arma::trimatu(mychol.t()), arma::solve(arma::trimatl(mychol), prec * arma::vectorise(X.t() * arma::diagmat(Gamma.col(i)) * (Y - X * Theta.t() * B.t() - condMean) * B)));
    
    Ltemp = mymean + arma::solve(arma::trimatl(mychol), arma::randn<arma::vec>(Lambda.n_rows * Lambda.n_cols));
    Lambda.slice(i) = arma::trans(arma::reshape(Ltemp,Lambda.n_cols,Lambda.n_rows));
    //Lambda.slice(i) = arma::orth(Lambda.slice(i));
    //Ltemp = mymean;
    //Lambda.slice(i) = arma::trans(arma::reshape(Ltemp,Lambda.n_cols,Lambda.n_rows));
    counter++;
    condMean.zeros();
   }

}

// Doesn't work so well
void updateLambda3(arma::mat& Y, arma::cube& Lambda, arma::mat& Tau, arma::mat& Eta, arma::mat& X, arma::mat& B, double prec, arma::mat& Theta){
  arma::mat P = getPenalty2(Theta.n_rows, 2);
  arma::uword n = Y.n_rows;
  arma::uword d = Theta.n_cols;
  arma::mat V = arma::zeros<arma::mat>(Theta.n_rows * d, Theta.n_rows * d);
  arma::vec Mu = arma::zeros<arma::vec>(Theta.n_rows * d);
  for(arma::uword k_curr = 0; k_curr < Lambda.n_slices; k_curr++){
    for(arma::uword i = 0; i < n; i++){
      arma::mat Vx = arma::zeros<arma::mat>(B.n_rows, B.n_rows);
      for(arma::uword k = 0; k < Lambda.n_slices; k++){
        if(k != k_curr){
          Vx = Vx + B * Lambda.slice(k) * X.row(i).t() * X.row(i) * Lambda.slice(k).t() * B.t();
        }
      }
      Vx = Vx + 1/prec * arma::eye(B.n_rows, B.n_rows);
      arma::mat Vx_inv = arma::inv_sympd(Vx);
      V = V + arma::kron(X.row(i).t() * Eta(i, k_curr), B.t()) * Vx_inv * arma::kron(X.row(i) * Eta(i, k_curr), B);
      Mu = Mu + arma::kron(X.row(i).t() * Eta(i, k_curr), B.t()) * Vx_inv * (Y.row(i).t() - (B * Theta * X.row(i).t()));
    }
    //V = V + arma::kron(P, arma::diagmat(Tau.row(0)));
    //V = V + arma::kron(P, arma::diagmat(Tau.row(k_curr + 1)));
    arma::mat mymean = arma::solve(V, Mu);
    //arma::vec mymean = arma::solve(arma::trimatu(mychol.t()), arma::solve(arma::trimatl(mychol), Mu));
    Lambda.slice(k_curr) = arma::reshape(arma::mvnrnd(mymean, arma::inv_sympd(V)), Theta.n_rows, Theta.n_cols);
    V.zeros();
    Mu.zeros();
  }
}

// [[Rcpp::export]]
void updateTheta(arma::mat& Y, arma::cube& Lambda, arma::mat& Tau, arma::mat& Gamma, arma::mat& X, arma::mat& B, double prec, arma::mat& Theta){
  arma::mat P = getPenalty2(Theta.n_rows, 2);
  arma::mat PriorPrecision = arma::kron(P, arma::diagmat(Tau.row(0)));
  arma::mat Precision = prec * arma::kron(B.t() * B, X.t() * X) + PriorPrecision;
  arma::mat mychol = arma::chol(Precision, "lower");
  arma::mat partialmean = X.t() * Y * B;
  for(arma::uword i = 0; i < Lambda.n_slices; i++){
    partialmean = partialmean - X.t() * arma::diagmat(Gamma.col(i)) * X * Lambda.slice(i).t() * B.t() * B;
  }
  arma::vec mymean = arma::solve(arma::trimatu(mychol.t()), arma::solve(arma::trimatl(mychol), prec * arma::vectorise(partialmean)));
  Theta = arma::trans(arma::reshape(mymean + arma::solve(arma::trimatl(mychol), arma::randn<arma::vec>(Theta.n_rows * Theta.n_cols)), Theta.n_cols, Theta.n_rows));
}

// Marginalizing over distribution of random effects...
// [[Rcpp::export]]
void updateTheta2(arma::mat& Y, arma::cube& Lambda, arma::mat& Tau, arma::mat& X, arma::mat& B, double prec, arma::mat& Theta){
  arma::mat P = getPenalty2(Theta.n_rows, 2);
  arma::uword n = Y.n_rows;
  arma::uword d = Theta.n_cols;
  arma::mat V = arma::zeros<arma::mat>(Theta.n_rows * d, Theta.n_rows * d);
  arma::vec Mu = arma::zeros<arma::vec>(Theta.n_rows * d);
  for(arma::uword i = 0; i < n; i++){
    arma::mat Vx = arma::zeros<arma::mat>(B.n_rows, B.n_rows);
    for(arma::uword k = 0; k < Lambda.n_slices; k++){
      Vx = Vx + B * Lambda.slice(k) * X.row(i).t() * X.row(i) * Lambda.slice(k).t() * B.t();
    }
    Vx = Vx + 1/prec * arma::eye(B.n_rows, B.n_rows);
    arma::mat Vx_inv = arma::inv_sympd(Vx);
    V = V + arma::kron(X.row(i).t(), B.t()) * Vx_inv * arma::kron(X.row(i), B);
    Mu = Mu + arma::kron(X.row(i).t(), B.t()) * Vx_inv * Y.row(i).t();
  }
  //V = V + arma::kron(P, arma::diagmat(Tau.row(0)));
  V = V + arma::kron(arma::diagmat(Tau.row(0)), P);
  arma::mat mychol = arma::chol(V, "lower");
  arma::vec mymean = arma::solve(arma::trimatu(mychol.t()), arma::solve(arma::trimatl(mychol), Mu));
  Theta = arma::reshape(arma::mvnrnd(mymean, arma::inv_sympd(V)), Theta.n_rows, Theta.n_cols);
  //Rcpp::Rcout << arma::reshape(mymean, Theta.n_rows, Theta.n_cols);
  //Theta = arma::reshape(mymean + arma::solve(arma::trimatl(mychol), arma::randn(Theta.n_cols * Theta.n_rows)), Theta.n_rows, Theta.n_cols);
}

// [[Rcpp::export]]
void updateThetaLambda(arma::mat &Y, arma::cube& Lambda, arma::mat& Eta, arma::mat& Tau, arma::mat& X, arma::mat& B, double prec, arma::mat& Theta){
  arma::uword n = Y.n_rows;
  arma::uword p = Theta.n_rows;
  arma::uword d = Theta.n_cols;
  arma::uword K = Lambda.n_slices;
  arma::mat P = getPenalty2(p, 1);
  arma::mat X_eta(n, (1 + K) * d);
  arma::vec result((K+1)*d*p);
  X_eta.cols(0, d - 1) = X;
  for(arma::uword k = 0; k < K; k++){
    X_eta.cols(d*(k+1), d*(k+2)-1) = arma::diagmat(Eta.col(k)) * X;
  }
  arma::mat Prior = arma::kron(arma::diagmat(arma::vectorise(Tau.t())), P);
  arma::mat Precision = prec * arma::kron(X_eta.t() * X_eta, B.t() * B) + Prior;
  arma::mat mychol = arma::chol(Precision, "lower");
  arma::vec w = solve(arma::trimatl(mychol), prec * arma::vectorise(B.t() * Y.t() * X_eta));
  arma::vec mu = solve(arma::trimatu(mychol.t()), w);
  arma::vec z = arma::randn<arma::vec>(p*d*(K+1));
  arma::vec v = arma::solve(arma::trimatu(mychol.t()), z);
  result = mu + v;
  //arma::vec mymean = arma::solve(arma::trimatu(mychol.t()), arma::solve(arma::trimatl(mychol), prec * arma::vectorise(B.t() * Y.t() * X_eta)));
  //result = mymean + arma::solve(arma::trimatl(mychol), arma::randn((1+K)*p*d));
  //result = arma::mvnrnd(arma::inv_sympd(Precision) * prec * arma::vectorise(B.t() * Y.t() * X_eta), arma::inv_sympd(Precision));
  Theta = arma::reshape(result.rows(0, p*d-1), p, d);
  for(arma::uword k = 0; k < K; k++){
    Lambda.slice(k) = arma::reshape(result.rows(p*d*(k+1), p*d*(k+2) - 1), p, d);
  }
}

// [[Rcpp::export]]
void updateThetaLambdaMH(arma::mat& Y, arma::mat& Theta, arma::cube& Lambda, arma::mat& Tau, double prec, arma::mat& X, arma::mat& B, double noise, arma::uword n){
  arma::uword D = X.n_cols;
  arma::mat Theta_Proposal = Theta;
  arma::cube Lambda_Proposal = Lambda;
  double J, P, A;
  for(arma::uword d = 0; d < D; d++){
    arma::mat Lambda_temp = Lambda.col(d);
    double prior_prop_lambda = 0;
    double prior_curr_lambda = 0;
    //Rcpp::List Proposals = Proposal(Theta.col(d), arma::mat(Lambda.col(d)), noise, n);
    //arma::vec new_theta = Proposals["Theta_Proposal"];
    //arma::mat new_lambda = Proposals["Lambda_Proposal"];
    //arma::vec new_theta = Theta.col(d) + noise * arma::randn<arma::vec>(Lambda.n_rows);
    arma::mat new_lambda = Lambda_temp + noise * arma::randn<arma::mat>(Lambda.n_rows, Lambda.n_slices);
    
    //Theta_Proposal.col(d) = new_theta;
    Lambda_Proposal.col(d) = new_lambda;
    //J = cpploglik_bayes(Theta_Proposal, Lambda, prec, X, B, Y, 12);
    
    if(d == 0){
    //  P = cpploglik_bayes(Theta, Lambda, prec, X, B, Y, 12);
    }
    /*
    A = J - P -1/2 * Tau(0, d) * arma::as_scalar(Theta.col(d).t() * getPenalty2(Lambda.n_rows, 2) * Theta.col(d) + 1/2 *
      Theta_Proposal.col(d).t() * getPenalty2(Lambda.n_rows, 2) * Theta_Proposal.col(d));
    if(R::runif(0.0, 1.0) < exp(A)){
       Theta.col(d) = Theta_Proposal.col(d);
       P = J;
    } else{
      Theta_Proposal.col(d) = Theta.col(d);
    }
    */
    //J = cpploglik_bayes(Theta, Lambda_Proposal, prec, X, B, Y, 12);
    for(arma::uword k = 0; k < Lambda.n_slices; k++){
      prior_prop_lambda = prior_prop_lambda + Tau(k + 1, d) * arma::as_scalar(Lambda_Proposal.slice(k).col(d).t() *
        getPenalty2(Lambda.n_rows, 2) * Lambda_Proposal.slice(k).col(d));
      prior_curr_lambda = prior_curr_lambda + Tau(k + 1, d) * arma::as_scalar(Lambda.slice(k).col(d).t() *
        getPenalty2(Lambda.n_rows, 2) * Lambda.slice(k).col(d));
    }
    A = J - P -1/2 * prior_curr_lambda + 1/2 * prior_prop_lambda;
    if(R::runif(0.0, 1.0) < exp(A)){
      Lambda.col(d) = Lambda_Proposal.col(d);
      P = J;
    }
    else{
      Lambda_Proposal.col(d) = Lambda.col(d);
    }
  }
  //arma::mat Lambda_current.col
  //Rcpp::List Proposals = Proposal(Theta.col(0), Lambda.slice(0).col(0));
  //arma::vec Theta_Proposal = d["Theta_Proposal"];
}
// [[Rcpp::export]]
void updateEta(arma::mat& Y, arma::cube& Lambda, arma::vec& Sigma, arma::mat& Eta, arma::mat& X, arma::mat& B, double prec, arma::mat& Theta){
  arma::uword nsubj = Y.n_rows;
  arma::mat Xtilde(Y.n_cols, Lambda.n_slices);
  arma::mat PriorPrecision = arma::eye<arma::mat>(Lambda.n_slices, Lambda.n_slices);
  arma::mat Precision(Lambda.n_slices, Lambda.n_slices);
  arma::mat mychol;
  for(arma::uword i = 0; i < nsubj; i++){
    for(arma::uword k = 0; k < Lambda.n_slices; k++){
      Xtilde.col(k) = B * Lambda.slice(k) * X.row(i).t();
    }
    Precision = prec * Xtilde.t() * Xtilde + PriorPrecision;/*
    arma::mat mychol = arma::chol(Precision, "lower");
    arma::vec w = arma::solve(arma::trimatl(mychol), prec * Xtilde.t() * (Y.row(i).t() - B * Theta * X.row(i).t()));
    arma::vec mu = arma::solve(trimatu(mychol.t()), w);
    arma::vec z = arma::randn<arma::vec>(Lambda.n_slices);
    arma::vec v = arma::solve(trimatu(mychol.t()), z);
    Eta.row(i) = arma::trans(mu + v);*/
    Eta.row(i) = arma::mvnrnd(prec * arma::inv_sympd(Precision) * Xtilde.t() * (Y.row(i).t() - B * Theta * X.row(i).t()), arma::inv_sympd(Precision)).t();
    //mychol = arma::chol(Precision, "lower");
    //Eta.row(i) = arma::trans(arma::solve(arma::trimatu(mychol.t()), arma::solve(arma::trimatl(mychol), prec * Xtilde.t() * (Y.row(i).t() - B * Theta * X.row(i).t()))) +
     //+ arma::solve(mychol, arma::randn<arma::vec>(Lambda.n_slices)));
  }
  //for(arma::uword i = 0; i < Eta.n_cols; i++){
  //  Eta.col(i) = Eta.col(i) / arma::stddev(Eta.col(i));
  //}
}

// [[Rcpp::export]]
void updateEta2(arma::mat& Y, arma::cube& Lambda, arma::vec& Sigma, arma::mat& Eta, arma::mat& X, arma::mat& B, double prec, arma::mat& Theta){
  arma::uword n = Y.n_rows;
  arma::uword K = Lambda.n_slices;
  arma::mat Ltilde(Lambda.n_rows, K);
  arma::vec EtaM(K);
  arma::mat EtaV(K, K);
  for(arma::uword i = 0; i < n; i++){
    Ltilde.zeros();
    for(arma::uword k = 0; k < K; k++){
      Ltilde.col(k) = Lambda.slice(k) * arma::trans(X.row(i));
    }
    EtaV = arma::inv_sympd(arma::eye<arma::mat>(K, K) + prec * arma::trans(Ltilde) * arma::trans(B) * B * Ltilde);
    EtaM = prec * EtaV * arma::trans(Ltilde) * arma::trans(B) * (arma::trans(Y.row(i)) - B * Theta * arma::trans(X.row(i)));
    Eta.row(i) = arma::mvnrnd(EtaM, EtaV).t();
  }
}

// [[Rcpp::export]]
void updateEta3(arma::mat& Y, arma::cube& Lambda, arma::mat& Eta, arma::mat& X, arma::mat& B, double prec, arma::mat& Theta){
  arma::uword n = Y.n_rows;
  arma::uword K = Lambda.n_slices;
  arma::vec fullmean;
  double fullprec;
  arma::vec condmean;
  double condprec;
  for(arma::uword i = 0; i < n; i++){
    fullmean = B * Theta * X.row(i).t();
    fullprec = 1;
    for(arma::uword kt = 0; kt < K; kt++){
      fullmean = fullmean + B * Lambda.slice(kt) * X.row(i).t() * Eta(i, kt);
      //fullprec = fullprec + prec * arma::as_scalar(X.row(i) * Lambda.slice(kt).t() * B.t() * B * Lambda.slice(kt) * X.row(i).t());
    }
    for(arma::uword k = 0; k < K; k++){
      condmean = fullmean - B * Lambda.slice(k) * X.row(i).t() * Eta(i, k);
      //condprec = fullprec - prec * arma::as_scalar(X.row(i) * Lambda.slice(k).t() * B.t() * B * Lambda.slice(k) * X.row(i).t());
      condprec = 1 + prec * arma::as_scalar(X.row(i) * Lambda.slice(k).t() * B.t() * B * Lambda.slice(k) * X.row(i).t());
      Eta(i, k) = arma::as_scalar(prec/condprec * X.row(i) * Lambda.slice(k).t() * B.t() * (Y.row(i).t()-condmean)) + 1/std::sqrt(double(condprec)) * R::rnorm(0, 1);
      fullmean = condmean + B * Lambda.slice(k) * X.row(i).t() * Eta(i, k);
    }
  }
}
// [[Rcpp::export]]
double updatePrec(arma::mat& Y, arma::cube& Lambda, arma::mat Gamma, arma::mat& X, arma::mat& B, arma::mat& Theta){
  arma::uword nsubj = Y.n_rows;
  arma::uword k = Lambda.n_slices;
  arma::vec mean;
  //double a = .0001;
  //double b = .0001;
  double a = .00001;
  double b = .00001;
  double mycumsum = 0;
  double totalN = Y.n_elem;
  for(arma::uword i = 0; i < nsubj; i++){
    mean = B * Theta * X.row(i).t();
    for(arma::uword j = 0; j < k; j++){
      mean = mean + Gamma(i, j) * B * Lambda.slice(j) * X.row(i).t();
    }
    mycumsum = mycumsum + arma::sum((Y.row(i).t() - mean) % (Y.row(i).t() - mean));
  }
  return(R::rgamma(a + totalN/2, 1.0 / (b + 1.0 / 2.0 * mycumsum)));
}



// [[Rcpp::export]]
void updateTau(arma::mat& Theta, arma::cube& Lambda, arma::mat& Tau){
  double t_alpha = 1;
  double t_beta = .000000005;
  //double t_beta = 1;
  arma::mat P = getPenalty2(Lambda.n_rows, 2);
  arma::uword R = Lambda.n_slices;
  arma::uword D = Lambda.n_cols;
  for(arma::uword d = 0; d < D; d++){
    Tau(0, d) = R::rgamma(t_alpha + (Lambda.n_rows - 1.0) / 2.0, 1.0 / (t_beta + 1.0 / 2.0 * arma::as_scalar(Theta.col(d).t() * P * Theta.col(d))));
  }
  /*
  for(arma::uword r = 0; r < R; r++){
    for(arma::uword d = 0; d < D; d++){
      Tau(r + 1, d) = R::rgamma(t_alpha + (Lambda.n_rows - 1.0) / 2.0, 1.0 / (t_beta + 1.0 / 2.0 * arma::as_scalar(Lambda.slice(r).col(d).t() * P * Lambda.slice(r).col(d))));
    }
  }
  */
  double temp_beta;
  for(arma::uword d = 0; d < D; d++){
    temp_beta = 0;
    for(arma::uword r = 1; r < R; r++){
      temp_beta = temp_beta + 1.0 / 2.0 * arma::as_scalar(Lambda.slice(r).col(d).t() * P * Lambda.slice(r).col(d));
    }
    Tau(1, d) = R::rgamma(t_alpha + R * (Lambda.n_rows - 1.0) / 2.0, 1.0 / (t_beta + temp_beta));
    for(arma::uword r = 0; r < R; r++){
      Tau(r + 1, d) = Tau(1, d);
    }
  }
}

// [[Rcpp::export]]
void updateSigma(arma::mat& Eta, arma::vec& Sigma){
  arma::uword K = Sigma.n_elem;
  double n = Eta.n_rows;
  double a = 1;
  double b = 1;
  for(arma::uword k = 0; k < K; k++){
    //Sigma(k) = gam_trunc_left(a + n / 2.0, 1.0 / (b + 1.0 / 2.0 * arma::sum(arma::square(Eta.col(k)))), Sigma(k - 1));
    Sigma(k) = R::rgamma(a + n / 2.0, 1.0 / (b + 1.0 / 2.0 * arma::sum(arma::square(Eta.col(k)))));
  }
}

void updatec(arma::cube& Lambda, arma::mat& c){
  double a = .1;
  double b = .1;
  arma::uword K = Lambda.n_slices;
  arma::uword D = Lambda.n_cols;
  for(arma::uword k = 0; k < K; k++){
    for(arma::uword d = 0; d < D; d++){
      c(k, d) = R::rgamma(a + Lambda.n_rows / 2.0, 1.0 / (b + 1.0 / 2.0 * arma::sum(arma::square(Lambda.slice(k).col(d)))));
      /*
      if(k == 0){
        c(k, d) = R::rgamma(a + Lambda.n_rows / 2.0, 1.0 / (b + 1.0 / 2.0 * arma::sum(arma::square(Lambda.slice(k).col(d)))));
      } else {
      c(k, d) = gam_trunc_left(a + Lambda.n_rows / 2.0, 1.0 / (b + 1.0 / 2.0 * arma::sum(arma::square(Lambda.slice(k).col(d)))), c(k - 1, d));
      }
      */
    }
  }
}

double gam_trunc_left(double a, double b,  double cut){ 
  double u, pcut, y; 
  pcut = R::pgamma(cut,a,b,1,0);
  if(pcut>0.99){
    return(cut);
  } 
  u = R::runif(0.0,1.0); 
  u = pcut + (1-pcut)*u; 
  y = R::qgamma(u, a, b, 1, 0);
  return y; 
} 

// [[Rcpp::export]]
void updateSigBeta(arma::vec& sigma, arma::vec& SigBeta, double Phi, arma::mat& X){
  double llsum = 0;
  double llsum_prop = 0;
  arma::vec SigBetaProp = SigBeta;
  double A;
  for(arma::uword i = 0; i < sigma.n_elem; i++){
    llsum = llsum + R::dgamma(sigma(i),
                              std::pow(std::exp(arma::dot(X.row(i).t(), SigBeta)), 2.0) / Phi,
                              Phi / std::exp(arma::dot(X.row(i).t(), SigBeta)),
                              true);
  }
  for(arma::uword d = 0; d < X.n_cols; d++){
    SigBetaProp(d) = SigBetaProp(d) + R::rnorm(0, .1);
    Rcpp::Rcout << SigBetaProp(d) << std::endl;
    for(arma::uword i = 0; i < sigma.n_elem; i++){
      llsum_prop = llsum_prop + R::dgamma(sigma(i),
                                          std::pow(std::exp(arma::dot(X.row(i).t(), SigBetaProp)), 2.0) / Phi,
                                          Phi / std::exp(arma::dot(X.row(i).t(), SigBetaProp)),
                                          true);
    }
    A = (llsum_prop - llsum) +
      R::dnorm(SigBeta(d), 0, .001, true) -
      R::dnorm(SigBetaProp(d), 0, .001, true);
    if(R::runif(0.0, 1.0) < exp(A)){
      Rcpp::Rcout << "Changed at " << d << std::endl;
      SigBeta(d) = SigBetaProp(d);
      llsum = llsum_prop;
    }
    
  }
}