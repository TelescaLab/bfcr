#include <RcppArmadillo.h>
#include <omp.h>
#include "updateParam.h"
#include "Utility.h"
// [[Rcpp::plugins(openmp)]]
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]

//[[Rcpp::export]]
arma::mat get_posterior_predictive_bands(List mod, arma::vec quantiles){
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

// [[Rcpp::export]]
arma::mat get_posterior_predictive_bands2(List mod, arma::vec quantiles){
  arma::field<arma::cube> BetaF = mod["Beta"];
  arma::field<arma::mat> PrecF = mod["Prec"];
  arma::field<arma::cube> LambdaF = mod["Lambda"];
  arma::field<arma::cube> EtaF = mod["Eta"];
  arma::mat Ymat = mod["Y"];
  arma::mat B = mod["B"];
  arma::mat X = mod["X"];
  arma::mat Z = mod["Z"];
  arma::vec Yvec = arma::vectorise(arma::trans(Ymat));
  arma::uword nchains = arma::size(BetaF)(0);
  arma::uword iter = BetaF(0).n_slices;
  arma::uword subjects = Ymat.n_rows;
  arma::uword time_points = B.n_rows;
  arma::vec ID(subjects * time_points);
  arma::vec time(subjects * time_points);
  arma::mat predictions(subjects * time_points, iter * nchains);
  arma::mat means(subjects * time_points, iter * nchains);
  arma::mat predictions_quant(subjects * time_points, arma::size(quantiles)(0));
  arma::mat means_quant(subjects * time_points, arma::size(quantiles)(0));
  arma::mat fit;
  
  for(arma::uword u = 0; u < nchains; u++){
    for(arma::uword i = 0; i < iter; i++){
      fit = X * BetaF(u, 0, 0).slice(i).t() * B.t();
      for(arma::uword k = 0; k < LambdaF(0, 0, 0).n_slices; k++){
        fit = fit + arma::diagmat(EtaF(u, 0, 0).slice(i).col(k)) * Z * LambdaF(u * iter + i, 0, 0).slice(k).t() * B.t();
      }
      means.col(u * iter + i) = arma::vectorise(fit.t());
      for(arma::uword s = 0; s < subjects; s++){
        
        predictions.col(u * iter + i).subvec(s * time_points, (s + 1) * time_points - 1) = 
          means.col(u * iter + i).subvec(s * time_points, (s + 1) * time_points - 1) + arma::randn<arma::vec>(time_points) * std::pow(PrecF(u,0,0)(s, i), -1.0 / 2.0);
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

// [[Rcpp::export]]
arma::mat get_posterior_means(List mod, arma::vec xi, double alpha){
  arma::field<arma::cube> BetaF = mod["Beta"];
  arma::mat B = mod["B"];
  arma::uword nchains = arma::size(BetaF)(0);
  arma::uword iter = BetaF(0).n_slices;
  arma::vec Malpha(iter * nchains);
  arma::running_stat_vec<arma::vec> stats;
  arma::vec Mbeta(BetaF(0).n_rows);
  arma::mat quants(B.n_rows, 3);
  for(arma::uword u = 0; u < nchains; u++){
    for(arma::uword i = 0; i < iter; i++){
      //Mbeta = Mbeta + BetaF(u).slice(i) * xi;
      stats(BetaF(u).slice(i) * xi);
    }
  }
  for(arma::uword u = 0; u < nchains; u++){
    for(arma::uword i = 0; i < iter; i++){
      Malpha(u * iter + i) = arma::max(arma::abs(BetaF(u).slice(i) * xi - stats.mean()) / stats.stddev());
    }
  }
  arma::vec alpha_vec = {1 - alpha};
  double q_alpha = arma::as_scalar(arma::quantile(Malpha, alpha_vec));
  
  quants.col(0) = B * stats.mean() - q_alpha * B * stats.stddev();
  quants.col(1) = B * stats.mean();
  quants.col(2) = B * stats.mean() + q_alpha * B * stats.stddev();
  //Rcout << "mean = " << std::endl << stats.mean() << std::endl;
  //Rcout << "var = " << std::endl << stats.stddev() << std::endl;
  return(quants);
}

// [[Rcpp::export]]
List extract_eigenfn(arma::cube& Lambda, const arma::vec& Delta,
                     arma::mat& Psi, arma::mat& Psi_sqrt,
                     arma::mat& Psi_sqrt_inv, arma::mat& B,
                     arma::uword eigenvals, arma::vec z){
  arma::uword dim_latent = Lambda.n_rows;
  arma::uword dim_spline = B.n_rows;
  arma::mat eigenfn_latent(dim_latent, dim_latent);
  arma::vec eigenval_latent(dim_latent);
  arma::mat eigenfn_latent_ordered(dim_latent, eigenvals);
  arma::mat eigenfn_spline(dim_spline, eigenvals);
  arma::vec eigenval_spline(eigenvals);
  arma::vec eigenval_pve(eigenvals);
  arma::mat cov_latent = arma::zeros<arma::mat>(dim_latent, dim_latent);
  for(arma::uword k = 0; k < Lambda.n_slices; k++){
    cov_latent = cov_latent + Lambda.slice(k) * z * z.t() * Lambda.slice(k).t();
  }
  cov_latent = cov_latent + arma::diagmat(1.0 / Delta);
  arma::mat cov_latent_transformed = Psi_sqrt * cov_latent * Psi_sqrt;
  arma::eig_sym(eigenval_latent, eigenfn_latent, cov_latent_transformed);
  eigenfn_latent =  Psi_sqrt_inv * eigenfn_latent;

  for(arma::uword v = 0; v < eigenvals; v++){
    eigenfn_latent_ordered.col(v) = eigenfn_latent.col(dim_latent - v - 1);

    //eigenval_spline(v) = arma::as_scalar(eigenfn_latent.col(dim_latent - v - 1).t() *
    //  Psi * cov_latent * Psi * eigenfn_latent.col(dim_latent - v - 1));
    eigenval_spline(v) = eigenval_latent(dim_latent - v - 1);
    eigenval_pve(v) = eigenval_spline(v) / arma::sum(eigenval_latent);
  }
  double magnitude = arma::sum(eigenval_latent);
  return(List::create(Named("eigenfn_latent", eigenfn_latent_ordered),
                      Named("eigenval", eigenval_spline),
                      Named("eigenval_pve", eigenval_pve),
                      Named("cov_latent", cov_latent),
                      Named("magnitude", magnitude)));
}

// [[Rcpp::export]]
List get_posterior_eigen(List mod, arma::uword eigenvals, arma::vec zi, double alpha){
  arma::mat B = mod["B"];
  arma::field<arma::cube> LambdaF = mod["Lambda"];
  arma::field<arma::mat> DeltaF = mod["Delta"];
  arma::vec Time = mod["Time"];
  arma::uword nchains = arma::size(DeltaF)(0);
  arma::uword iter = DeltaF(0).n_cols;
  arma::uword dim = B.n_cols;
  
  arma::mat Malpha(iter * nchains, eigenvals);
  arma::running_stat_vec<arma::vec> stats_vec;
  arma::running_stat_vec<arma::vec> stats_val;
  arma::running_stat_vec<arma::vec> stats_cov;
  arma::mat Psi(dim, dim);
  arma::mat Psi_sqrt;
  arma::mat Psi_sqrt_inv;
  for(arma::uword j = 0; j < dim; j++){
    for(arma::uword i = 0; i < dim; i++){
      Psi(i, j) = arma::as_scalar(arma::trapz(Time, B.col(i) % B.col(j)));
    }
  }
  Psi_sqrt = arma::sqrtmat_sympd(Psi);
  Psi_sqrt_inv = arma::inv_sympd(Psi_sqrt);
  
  arma::mat temp_evec;
  arma::mat eval_mat(nchains * iter, eigenvals);
  arma::mat eval_pve_mat(nchains * iter, eigenvals);
  List eigen_list;
  arma::vec magnitude(nchains * iter);
  // Create running mean and standard deviations of eigenvectors

  arma::uword idx1, idx2;
  for(arma::uword u = 0; u < nchains; u++){
    for(arma::uword i = 0; i < iter; i++){
      eigen_list = extract_eigenfn(LambdaF(u * iter + i, 0, 0),
                      DeltaF(u, 0, 0).col(i),
                      Psi,
                      Psi_sqrt,
                      Psi_sqrt_inv,
                      B,
                      eigenvals,
                      zi);
      temp_evec = Rcpp::as<arma::mat>(eigen_list["eigenfn_latent"]);
      eval_mat.row(u * iter + i) = Rcpp::as<arma::rowvec>(eigen_list["eigenval"]);
      eval_pve_mat.row(u * iter + i) = Rcpp::as<arma::rowvec>(eigen_list["eigenval_pve"]);
      stats_cov(arma::vectorise(as<arma::vec>(eigen_list["cov_latent"])));
      magnitude(u * iter + i) = eigen_list["magnitude"];
      if ((i == 0) && (u == 0)) {
        stats_vec(arma::vectorise(temp_evec));
      } else {
        // align eigenvectors
        
        for (arma::uword k = 0; k < eigenvals; k++) {
          idx1 = k * dim;
          idx2 = (k + 1) * dim - 1;
          if (arma::sum(arma::square(temp_evec.col(k) + stats_vec.mean().subvec(idx1, idx2))) <
            arma::sum(arma::square(temp_evec.col(k) - stats_vec.mean().subvec(idx1, idx2)))) {
            temp_evec.col(k) = -temp_evec.col(k);
          }
        }

        stats_vec(arma::vectorise(temp_evec));
      }
    }
  }
  for (arma::uword u = 0; u < nchains; u++) {
    for (arma::uword i = 0; i < iter; i++) {
      eigen_list = extract_eigenfn(LambdaF(u * iter + i, 0, 0),
                                   DeltaF(u, 0, 0).col(i),
                                   Psi,
                                   Psi_sqrt,
                                   Psi_sqrt_inv,
                                   B,
                                   eigenvals,
                                   zi);
      temp_evec = Rcpp::as<arma::mat>(eigen_list["eigenfn_latent"]);

      for (arma::uword k = 0; k < eigenvals; k++) {
        idx1 = k * dim;
        idx2 = (k + 1) * dim - 1;
        if (arma::sum(arma::square(temp_evec.col(k) + stats_vec.mean().subvec(idx1, idx2))) <
          arma::sum(arma::square(temp_evec.col(k) - stats_vec.mean().subvec(idx1, idx2)))) {
          temp_evec.col(k) = -temp_evec.col(k);
        }
        Malpha(u * iter + i, k) = arma::max((temp_evec.col(k) -
          stats_vec.mean().subvec(idx1, idx2)) / stats_vec.stddev().subvec(idx1, idx2));
      }
    }
  }
  arma::vec alpha_vec = {1 - alpha};
  arma::vec alpha_vec_eval = {alpha / 2.0, 0.5, 1 - alpha / 2.0};
  double q_alpha;
  arma::mat lower(Time.n_elem, eigenvals);
  arma::mat mean(Time.n_elem, eigenvals);
  arma::mat upper(Time.n_elem, eigenvals);
  arma::mat eigenbands(Time.n_elem, eigenvals * 3);
  arma::mat eigenval_intervals(3, eigenvals);
  arma::mat eigenval_pve_intervals(3, eigenvals);
  arma::vec magnitude_interval(3);
  for (arma::uword k = 0; k < eigenvals; k++) {
    idx1 = k * dim;
    idx2 = (k + 1) * dim - 1;
    q_alpha = arma::as_scalar(arma::quantile(Malpha.col(k), alpha_vec));
    lower.col(k) = B * (stats_vec.mean().subvec(idx1, idx2) -
      q_alpha * stats_vec.stddev().subvec(idx1, idx2));
    mean.col(k) = B * (stats_vec.mean().subvec(idx1, idx2));
    upper.col(k) = B * (stats_vec.mean().subvec(idx1, idx2) + 
      q_alpha * stats_vec.stddev().subvec(idx1, idx2));
    eigenval_intervals.col(k) = arma::quantile(eval_mat.col(k), alpha_vec_eval);
    eigenval_pve_intervals.col(k) = arma::quantile(eval_pve_mat.col(k), alpha_vec_eval);
  }
  magnitude_interval = arma::quantile(magnitude, alpha_vec_eval);
  arma::mat mean_cov = B * arma::reshape(stats_cov.mean(), dim, dim) * B.t();
  return(List::create(Named("lower", lower),
                      Named("mean", mean),
                      Named("upper", upper),
                      Named("eigenval_intervals", eigenval_intervals),
                      Named("eigenval_pve_intervals", eigenval_pve_intervals),
                      Named("surface", mean_cov),
                      Named("magnitude", magnitude_interval),
                      Named("raw_magnitude", magnitude)));
}


// [[Rcpp::export]]
List extract_eigenfn2(arma::cube& Lambda,
                     arma::mat& Psi, arma::mat& Psi_sqrt,
                     arma::mat& Psi_sqrt_inv, arma::mat& B,
                     arma::uword eigenvals, arma::vec z){
  arma::uword dim_latent = Lambda.n_rows;
  arma::uword dim_spline = B.n_rows;
  arma::mat eigenfn_latent(dim_latent, dim_latent);
  arma::vec eigenval_latent(dim_latent);
  arma::mat eigenfn_latent_ordered(dim_latent, eigenvals);
  arma::mat eigenfn_spline(dim_spline, eigenvals);
  arma::vec eigenval_spline(eigenvals);
  arma::vec eigenval_pve(eigenvals);
  arma::mat cov_latent = arma::zeros<arma::mat>(dim_latent, dim_latent);
  for(arma::uword k = 0; k < Lambda.n_slices; k++){
    cov_latent = cov_latent + Lambda.slice(k) * z * z.t() * Lambda.slice(k).t();
  }
  arma::mat cov_latent_transformed = Psi_sqrt * cov_latent * Psi_sqrt;
  arma::eig_sym(eigenval_latent, eigenfn_latent, cov_latent_transformed);
  eigenfn_latent =  Psi_sqrt_inv * eigenfn_latent;
  
  for(arma::uword v = 0; v < eigenvals; v++){
    eigenfn_latent_ordered.col(v) = eigenfn_latent.col(dim_latent - v - 1);
    
    //eigenval_spline(v) = arma::as_scalar(eigenfn_latent.col(dim_latent - v - 1).t() *
    //  Psi * cov_latent * Psi * eigenfn_latent.col(dim_latent - v - 1));
    eigenval_spline(v) = eigenval_latent(dim_latent - v - 1);
    eigenval_pve(v) = eigenval_spline(v) / arma::sum(eigenval_latent);
  }
  double magnitude = arma::sum(eigenval_latent);
  return(List::create(Named("eigenfn_latent", eigenfn_latent_ordered),
                      Named("eigenval", eigenval_spline),
                      Named("eigenval_pve", eigenval_pve),
                      Named("cov_latent", cov_latent),
                      Named("magnitude", magnitude)));
}

// [[Rcpp::export]]
List get_posterior_eigen2(List mod, arma::uword eigenvals, arma::vec zi, double alpha){
  arma::mat B = mod["B"];
  arma::field<arma::cube> LambdaF = mod["Lambda"];
  arma::field<arma::cube> BetaF = mod["Beta"];
  arma::vec Time = mod["Time"];
  arma::uword nchains = arma::size(BetaF)(0);
  arma::uword iter = BetaF(0, 0, 0).n_slices;
  arma::uword dim = B.n_cols;
  
  arma::mat Malpha(iter * nchains, eigenvals);
  arma::running_stat_vec<arma::vec> stats_vec;
  arma::running_stat_vec<arma::vec> stats_val;
  arma::running_stat_vec<arma::vec> stats_cov;
  arma::mat Psi(dim, dim);
  arma::mat Psi_sqrt;
  arma::mat Psi_sqrt_inv;
  for(arma::uword j = 0; j < dim; j++){
    for(arma::uword i = 0; i < dim; i++){
      Psi(i, j) = arma::as_scalar(arma::trapz(Time, B.col(i) % B.col(j)));
    }
  }
  Psi_sqrt = arma::sqrtmat_sympd(Psi);
  Psi_sqrt_inv = arma::inv_sympd(Psi_sqrt);
  
  arma::mat temp_evec;
  arma::mat eval_mat(nchains * iter, eigenvals);
  arma::mat eval_pve_mat(nchains * iter, eigenvals);
  List eigen_list;
  arma::vec magnitude(nchains * iter);
  // Create running mean and standard deviations of eigenvectors
  
  arma::uword idx1, idx2;
  for(arma::uword u = 0; u < nchains; u++){
    for(arma::uword i = 0; i < iter; i++){
      eigen_list = extract_eigenfn2(LambdaF(u * iter + i, 0, 0),
                                   Psi,
                                   Psi_sqrt,
                                   Psi_sqrt_inv,
                                   B,
                                   eigenvals,
                                   zi);
      temp_evec = Rcpp::as<arma::mat>(eigen_list["eigenfn_latent"]);
      eval_mat.row(u * iter + i) = Rcpp::as<arma::rowvec>(eigen_list["eigenval"]);
      eval_pve_mat.row(u * iter + i) = Rcpp::as<arma::rowvec>(eigen_list["eigenval_pve"]);
      stats_cov(arma::vectorise(as<arma::vec>(eigen_list["cov_latent"])));
      magnitude(u * iter + i) = eigen_list["magnitude"];
      if ((i == 0) && (u == 0)) {
        stats_vec(arma::vectorise(temp_evec));
      } else {
        // align eigenvectors
        
        for (arma::uword k = 0; k < eigenvals; k++) {
          idx1 = k * dim;
          idx2 = (k + 1) * dim - 1;
          if (arma::sum(arma::square(temp_evec.col(k) + stats_vec.mean().subvec(idx1, idx2))) <
            arma::sum(arma::square(temp_evec.col(k) - stats_vec.mean().subvec(idx1, idx2)))) {
            temp_evec.col(k) = -temp_evec.col(k);
          }
        }
        
        stats_vec(arma::vectorise(temp_evec));
      }
    }
  }
  for (arma::uword u = 0; u < nchains; u++) {
    for (arma::uword i = 0; i < iter; i++) {
      eigen_list = extract_eigenfn2(LambdaF(u * iter + i, 0, 0),
                                   Psi,
                                   Psi_sqrt,
                                   Psi_sqrt_inv,
                                   B,
                                   eigenvals,
                                   zi);
      temp_evec = Rcpp::as<arma::mat>(eigen_list["eigenfn_latent"]);
      
      for (arma::uword k = 0; k < eigenvals; k++) {
        idx1 = k * dim;
        idx2 = (k + 1) * dim - 1;
        if (arma::sum(arma::square(temp_evec.col(k) + stats_vec.mean().subvec(idx1, idx2))) <
          arma::sum(arma::square(temp_evec.col(k) - stats_vec.mean().subvec(idx1, idx2)))) {
          temp_evec.col(k) = -temp_evec.col(k);
        }
        Malpha(u * iter + i, k) = arma::max((temp_evec.col(k) -
          stats_vec.mean().subvec(idx1, idx2)) / stats_vec.stddev().subvec(idx1, idx2));
      }
    }
  }
  arma::vec alpha_vec = {1 - alpha};
  arma::vec alpha_vec_eval = {alpha / 2.0, 0.5, 1 - alpha / 2.0};
  double q_alpha;
  arma::mat lower(Time.n_elem, eigenvals);
  arma::mat mean(Time.n_elem, eigenvals);
  arma::mat upper(Time.n_elem, eigenvals);
  arma::mat eigenbands(Time.n_elem, eigenvals * 3);
  arma::mat eigenval_intervals(3, eigenvals);
  arma::mat eigenval_pve_intervals(3, eigenvals);
  arma::vec magnitude_interval(3);
  for (arma::uword k = 0; k < eigenvals; k++) {
    idx1 = k * dim;
    idx2 = (k + 1) * dim - 1;
    q_alpha = arma::as_scalar(arma::quantile(Malpha.col(k), alpha_vec));
    lower.col(k) = B * (stats_vec.mean().subvec(idx1, idx2) -
      q_alpha * stats_vec.stddev().subvec(idx1, idx2));
    mean.col(k) = B * (stats_vec.mean().subvec(idx1, idx2));
    upper.col(k) = B * (stats_vec.mean().subvec(idx1, idx2) + 
      q_alpha * stats_vec.stddev().subvec(idx1, idx2));
    eigenval_intervals.col(k) = arma::quantile(eval_mat.col(k), alpha_vec_eval);
    eigenval_pve_intervals.col(k) = arma::quantile(eval_pve_mat.col(k), alpha_vec_eval);
  }
  magnitude_interval = arma::quantile(magnitude, alpha_vec_eval);
  arma::mat mean_cov = B * arma::reshape(stats_cov.mean(), dim, dim) * B.t();
  return(List::create(Named("lower", lower),
                      Named("mean", mean),
                      Named("upper", upper),
                      Named("eigenval_intervals", eigenval_intervals),
                      Named("eigenval_pve_intervals", eigenval_pve_intervals),
                      Named("surface", mean_cov),
                      Named("magnitude", magnitude_interval),
                      Named("raw_magnitude", magnitude)));
}
