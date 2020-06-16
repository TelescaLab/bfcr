#include <RcppArmadillo.h>
#ifdef _OPENMP
 #include <omp.h>
#endif
// [[Rcpp::plugins(openmp)]]
using namespace Rcpp;

// [[Rcpp::export]]
arma::uvec armadillo_modulus2(arma::uvec indicies, arma::uword n){
  return(indicies - n * arma::floor(indicies / n));
}

// [[Rcpp::export]]
List get_omnibus_fit(List mod){
  arma::field<arma::cube> ThetaF = mod["Theta"];
  arma::field<arma::mat> PrecF = mod["Prec"];
  arma::mat Ymat = mod["Y"];
  arma::mat B = mod["B"];
  arma::uvec missing = arma::find_nonfinite(Ymat);
  arma::uvec missing_sub = armadillo_modulus2(missing, Ymat.n_rows);
  arma::uvec missing_time = arma::floor(missing / Ymat.n_rows);
  arma::uword nchains = arma::size(ThetaF)(0);
  arma::uword iter = ThetaF(0).n_slices;
  arma::uword subjects = Ymat.n_rows;
  arma::uword time_points = Ymat.n_cols;
  arma::vec statistic_rep(nchains * iter);
  arma::vec statistic_obs(nchains * iter);
  statistic_rep.zeros();
  statistic_obs.zeros();
  for(arma::uword u = 0; u < nchains; u++){
    for(arma::uword i = 0; i < iter; i++){
      arma::mat fit = ThetaF(u).slice(i) * B.t();
      for(arma::uword s = 0; s < subjects; s++){
        for(arma::uword m = 0; m < missing_sub.n_elem; m++){
          Ymat(missing_sub(m), missing_time(m)) = fit(missing_sub(m), missing_time(m)) + R::rnorm(0, std::pow(PrecF(u)(missing_sub(m), i), -.5));
        }
        arma::vec rep_temp = arma::square(arma::randn<arma::vec>(time_points));
        arma::rowvec obs_temp = arma::square(Ymat.row(s) -
          ThetaF(u).slice(i).row(s) * B.t());
        statistic_rep(u * iter + i) = statistic_rep(u * iter + i) + std::pow(PrecF(u)(s, i), -1.0 / 2.0) * arma::sum(rep_temp);
        statistic_obs(u * iter + i) = statistic_obs(u * iter + i) +
          std::pow(PrecF(u)(s, i), 1.0 / 2.0) * arma::sum(obs_temp);
        
      }
    }
  }
  return(List::create(Named("statistic_obs", statistic_obs),
                      Named("statistic_rep", statistic_rep)));
}

// [[Rcpp::export]]
List get_omnibus_fit2(List mod){
  arma::field<arma::cube> LambdaF = mod["Lambda"];
  arma::field<arma::cube> BetaF = mod["Beta"];
  arma::field<arma::mat> PrecF = mod["Prec"];
  arma::field<arma::cube> EtaF = mod["Eta"];
  
  arma::mat Ymat = mod["Y"];
  arma::mat B = mod["B"];
  arma::mat X = mod["X"];
  arma::mat Z = mod["Z"];
  arma::uvec missing = arma::find_nonfinite(Ymat);
  arma::uvec missing_sub = armadillo_modulus2(missing, Ymat.n_rows);
  arma::uvec missing_time = arma::floor(missing / Ymat.n_rows);
  arma::uword nchains = arma::size(BetaF)(0);
  arma::uword iter = BetaF(0).n_slices;
  arma::uword subjects = Ymat.n_rows;
  arma::uword time_points = Ymat.n_cols;
  arma::vec statistic_rep(nchains * iter);
  arma::vec statistic_obs(nchains * iter);
  statistic_rep.zeros();
  statistic_obs.zeros();
  arma::mat fit;
  for(arma::uword u = 0; u < nchains; u++){
    for(arma::uword i = 0; i < iter; i++){
      fit = X * BetaF(u, 0, 0).slice(i).t() * B.t();
      for(arma::uword k = 0; k < LambdaF(0, 0, 0).n_slices; k++){
        fit = fit + arma::diagmat(EtaF(u, 0, 0).slice(i).col(k)) * Z * LambdaF(i + u * iter, 0, 0).slice(k).t() * B.t();
      }
      for(arma::uword s = 0; s < subjects; s++){
        for(arma::uword m = 0; m < missing_sub.n_elem; m++){
          Ymat(missing_sub(m), missing_time(m)) = fit(missing_sub(m), missing_time(m)) + R::rnorm(0, std::pow(PrecF(u)(missing_sub(m), i), -.5));
        }
        arma::vec rep_temp = arma::square(arma::randn<arma::vec>(time_points));
        arma::rowvec obs_temp = arma::square(Ymat.row(s) -
          fit.row(s));
        statistic_rep(u * iter + i) = statistic_rep(u * iter + i) + std::pow(PrecF(u)(s, i), -1.0 / 2.0) * arma::sum(rep_temp);
        statistic_obs(u * iter + i) = statistic_obs(u * iter + i) +
          std::pow(PrecF(u)(s, i), 1.0 / 2.0) * arma::sum(obs_temp);
        
      }
    }
  }
  return(List::create(Named("statistic_obs", statistic_obs),
                      Named("statistic_rep", statistic_rep)));
}
// [[Rcpp::export]]
arma::mat DiffOp(arma::uword n){
  arma::mat D = arma::eye(n, n);
  for(arma::uword i = 1; i < n; i++){
    D(i, i - 1) = -1;
  }
  return(D);
}

// [[Rcpp::export]]
arma::mat getPenalty2(arma::uword n, arma::uword D){
  arma::mat DiffMat = DiffOp(n);
  DiffMat.rows(0, D-1).zeros();
  for(arma::uword d = 1; d < D; d++){
    DiffMat = DiffMat * DiffMat;
  }
  return(arma::trans(DiffMat) * DiffMat);
}

// [[Rcpp::export]]
arma::mat getPenalty(arma::uword n){
  arma::mat Penalty = arma::zeros<arma::mat>(n, n);
  for(arma::uword i = 0; i < n - 1; i++){
    Penalty(i + 1, i) = -1.0;
  }
  for(arma::uword i = 0; i < n - 1; i++){
    Penalty(i, i + 1) = -1.0;
  }
  for(arma::uword i = 1; i < n - 1; i++){
    Penalty(i, i) = 2.0;
  }
  Penalty(0, 0) = 1.0;
  Penalty(n - 1, n - 1) = 1.0;
  
  return(Penalty);
}

// [[Rcpp::export]]
arma::rowvec initializeY(arma::vec y, arma::vec observedTimes, arma::vec fullTimes){
  arma::uvec result(observedTimes.n_elem);
  arma::uword J = 0;
  for(arma::uword i = 0; i < observedTimes.n_elem; i++){
    for(arma::uword j = J; j < fullTimes.n_elem; j++){
      if(observedTimes(i) == fullTimes(j)){
        result(i) = j;
        J = j;
      }
    }
  }
  arma::vec yfull = arma::zeros<arma::vec>(fullTimes.n_elem);
  yfull(result) = y;
  return(yfull.t());
}

// [[Rcpp::export]]
arma::uvec getObservedOrder(arma::vec observedTimes, arma::vec fullTimes){
  arma::uvec result(observedTimes.n_elem);
  arma::uword J = 0;
  for(arma::uword i = 0; i < observedTimes.n_elem; i++){
    for(arma::uword j = J; j < fullTimes.n_elem; j++){
      if(observedTimes(i) == fullTimes(j)){
        result(i) = j;
        J = j;
      }
    }
  }
  return(result);
}

// [[Rcpp::export]]
void PredictY(arma::mat& ImputedY, arma::mat X, arma::mat B, arma::mat Theta, arma::mat Eta, arma::cube Lambda, double Prec){
  arma::mat mymean = X * Theta.t() * B.t();
  Rcpp::Rcout << arma::size(X) << std::endl;
  
  for(arma::uword k = 0; k < Lambda.n_slices; k++){
    mymean = mymean + arma::diagmat(Eta.col(k)) * X * Lambda.slice(k).t() * B.t();
  }
  mymean = mymean + (1.0 / sqrt(Prec)) * arma::randn(mymean.n_rows, mymean.n_rows);
  ImputedY = mymean;
}

// [[Rcpp::export]]
void PredictY2(arma::mat& ImputedY, arma::field<arma::uvec> observedOrder, arma::mat X, arma::mat B, arma::mat Theta, arma::cube Lambda, arma::mat Eta, double Prec){
  
  arma::mat temp = X * Theta.t() * B.t();
  for(arma::uword k = 0; k < Lambda.n_slices; k++){
    temp = temp + arma::diagmat(Eta.col(k)) * X * Lambda.slice(k).t() * B.t();
  }
  temp = temp + (1.0 / sqrt(Prec)) * arma::randn(X.n_rows, B.n_rows);
  for(arma::uword i = 0; i < X.n_rows; i++){
    for(arma::uword j = 0; j < observedOrder(i).n_elem; j++){
      temp(i, observedOrder(i)(j)) = ImputedY(i, observedOrder(i)(j));
    }
  }
  ImputedY = temp;
}
arma::uvec test(arma::field<arma::vec> observedTimes, arma::vec fullTimes){
  return(getObservedOrder(observedTimes(0), fullTimes));
}

// Propose a new column of Theta and new columns of Lambda 
// [[Rcpp::export]]
Rcpp::List Proposal(arma::vec Theta, arma::mat Lambda, double noise = .1, arma::uword samples = 200){
  arma::mat Y = arma::mvnrnd(Theta, Lambda * Lambda.t() + noise * arma::eye(Lambda.n_rows, Lambda.n_rows), samples);
  arma::mat coeff;
  arma::mat score;
  arma::vec latent;
  arma::princomp(coeff, score, latent, Y.t());
  arma::mat Lambda_Proposal = coeff.cols(0, Lambda.n_cols - 1) * arma::diagmat(arma::sqrt(latent.subvec(0, Lambda.n_cols - 1)));
  return(Rcpp::List::create(Rcpp::Named("Theta_Proposal", arma::mean(Y, 1)), Rcpp::Named("Lambda_Proposal", Lambda_Proposal)));
}

// [[Rcpp::export]]
double cpploglik_bayes(arma::mat &Theta, arma::cube &Lambda, double precision, arma::vec& Phi,
                 arma::mat &X, arma::mat &B, arma::mat &Y, int cores = 1){
  arma::uword K = Lambda.n_slices;
  arma::uword n = Y.n_rows;
  arma::uword tmax = Y.n_cols;
  arma::uword p = Lambda.n_rows;
  
  double loglik = 0;
  double constants = -double(tmax) / 2 * log(2 * PI);
#ifdef _OPENMP
  omp_set_num_threads(cores);
  #pragma omp parallel for reduction(+:loglik)
#endif
  for(arma::uword i = 0; i < n; i++){
    arma::vec mean(tmax);
    arma::mat cov(tmax, tmax);
    arma::mat LambdaCov = arma::zeros<arma::mat>(p, p);
    arma::mat precisionmat = 1/precision * arma::eye<arma::mat>(tmax, tmax) + B * arma::diagmat(1.0 / Phi) * B.t();
    arma::mat rooti;
    double rootisum;
    arma::mat rootsum;
    arma::vec z;
    for(arma::uword k = 0; k < K; k++){
      LambdaCov = Lambda.slice(k) * X.row(i).t() * X.row(i) * Lambda.slice(k).t() + LambdaCov;
    }
    mean = B * Theta * arma::trans(X.row(i));
    cov = B * LambdaCov * arma::trans(B) + precisionmat;
    LambdaCov.zeros();
    rooti = arma::trans(arma::inv(trimatu(arma::chol(cov))));
    rootisum = arma::sum(log(rooti.diag()));
    z = rooti * arma::trans(Y.row(i) - arma::trans(mean));
    loglik = constants - 0.5 * arma::sum(z % z) + rootisum + loglik;
  }
  
  return(loglik);
}

// [[Rcpp::export]]
void find_stepsize(arma::mat& Y, arma::mat& Theta, arma::cube& Lambda, double prec, arma::mat& X, arma::mat& B, double noise){
  arma::uword D = X.n_cols;
  arma::mat Theta_Proposal = Theta;
  arma::cube Lambda_Proposal = Lambda;
  double J, P;
  for(arma::uword d = 0; d < D; d++){
    arma::mat Lambda_temp = Lambda.col(d);
    arma::mat new_lambda = Lambda_temp + noise * arma::randn<arma::mat>(Lambda.n_rows, Lambda.n_slices);
    
    //Theta_Proposal.col(d) = new_theta;
    Lambda_Proposal.col(d) = new_lambda;
    //J = cpploglik_bayes(Theta_Proposal, Lambda, prec, X, B, Y, 12);
    
    if(d == 0){
      //P = cpploglik_bayes(Theta, Lambda, prec, X, B, Y, 12);
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
    Rcpp::Rcout << "Current likelihood: " << P << std::endl << "Proposal likelihood: " << J << std::endl << 
      "Probability of acceptance: " << exp(J-P)*100 << "%" << std::endl;

  }
}