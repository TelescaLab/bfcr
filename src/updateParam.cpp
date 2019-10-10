#include <RcppArmadillo.h>
#include <cmath>
#include "Utility.h"
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]

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
  arma::mat P = getPenalty2(p, 2);
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
      P = cpploglik_bayes(Theta, Lambda, prec, X, B, Y, 12);
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
    J = cpploglik_bayes(Theta, Lambda_Proposal, prec, X, B, Y, 12);
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
  double t_beta = .000005;
  //double t_beta = 1;
  arma::mat P = getPenalty2(Lambda.n_rows, 2);
  arma::uword R = Lambda.n_slices;
  arma::uword D = Lambda.n_cols;
  for(arma::uword d = 0; d < D; d++){
    Tau(0, d) = R::rgamma(t_alpha + (Lambda.n_rows - 1.0) / 2.0, 1.0 / (t_beta + 1.0 / 2.0 * arma::as_scalar(Theta.col(d).t() * P * Theta.col(d))));
  }
  for(arma::uword r = 0; r < R; r++){
    for(arma::uword d = 0; d < D; d++){
      Tau(r + 1, d) = R::rgamma(t_alpha + (Lambda.n_rows - 1.0) / 2.0, 1.0 / (t_beta + 1.0 / 2.0 * arma::as_scalar(Lambda.slice(r).col(d).t() * P * Lambda.slice(r).col(d))));
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