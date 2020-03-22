#include <RcppArmadillo.h>
#include <omp.h>
#include "updateParam.h"
#include "Utility.h"
// [[Rcpp::plugins(openmp)]]
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]

//[[Rcpp::export]]
arma::mat posterior_predictive_bands(List mod, arma::vec quantiles){
  arma::field<arma::cube> BetaF = mod["Beta"];
  arma::field<arma::cube> ThetaF = mod["Theta"];
  arma::field<arma::mat> PrecF = mod["Prec"];
  arma::mat Ymat = mod["Y"];
  arma::mat B = mod["B"];
  arma::vec Yvec = arma::vectorise(arma::trans(Ymat));
  arma::uword nchains = arma::size(ThetaF)(0);
  arma::uword iter = ThetaF(0).n_slices;
  arma::uword subjects = ThetaF(0).n_rows;
  arma::uword time_points = B.n_rows;
  arma::vec ID(subjects * time_points);
  arma::vec time(subjects * time_points);
  arma::mat predictions(subjects * time_points, iter * nchains);
  arma::mat means(subjects * time_points, iter * nchains);
  arma::mat predictions_quant(subjects * time_points, arma::size(quantiles)(0));
  arma::mat means_quant(subjects * time_points, arma::size(quantiles)(0));
  for(arma::uword u = 0; u < nchains; u++){
    for(arma::uword i = 0; i < iter; i++){
      for(arma::uword s = 0; s < subjects; s++){
        means.col(u * iter + i).subvec(s * time_points, (s + 1) * time_points - 1) = B * ThetaF(u).slice(i).row(s).t();
        predictions.col(u * iter + i).subvec(s * time_points, (s + 1) * time_points - 1) = 
          means.col(u * iter + i).subvec(s * time_points, (s + 1) * time_points - 1) + arma::randn<arma::vec>(time_points) * std::pow(PrecF(u)(s, i), -1.0 / 2.0);
      }
    }
  }
  for(arma::uword s = 0; s < subjects; s++){
    for(arma::uword t = 0; t < time_points; t++){
      means_quant.row(s * time_points + t) = arma::quantile(means.row(s * time_points + t), quantiles);
      predictions_quant.row(s * time_points + t) = arma::quantile(predictions.row(s * time_points + t), quantiles);
      time(s * time_points + t) = t;
      ID(s * time_points + t) = s;
    }
  }
  
  return(arma::join_rows(arma::join_rows(ID, time), Yvec, predictions_quant, means_quant));
}