#include <RcppArmadillo.h>
#include "Utils.h"

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
Rcpp::List get_posterior_subject_bands_cpp(List mcmc_output,
                                           double alpha){
  Rcpp::List samples = mcmc_output["samples"];
  Rcpp::List data = mcmc_output["data"];
  Rcpp::List control = mcmc_output["control"];
  arma::cube beta = samples["beta"];
  arma::mat varphi = samples["varphi"];
  arma::field<arma::cube> lambda = samples["lambda"];
  arma::cube eta = samples["eta"];
  arma::uword kdim = data["latent_dimension"];
  arma::mat basis = data["basis"];
  arma::mat response = data["response"];
  arma::mat design_mean = data["design_mean"];
  arma::mat design_var = data["design_var"];
  arma::uword iterations = control["iterations"];
  arma::uword burnin = control["burnin"];
  arma::uword num_subjects = response.n_rows;
  arma::running_stat_vec<arma::vec> stats;
  arma::mat m_alpha(num_subjects, iterations - burnin);
  arma::mat current_lower_dim_fit;
  arma::mat current_fit;
  
  arma::vec response_vectorized = arma::vectorise(arma::trans(response));
  for (arma::uword iter = burnin; iter < iterations; iter++) {
    current_lower_dim_fit = beta.slice(iter) * design_mean.t();
    for (arma::uword k = 0; k < kdim; k++) {
      current_lower_dim_fit = current_lower_dim_fit + 
        lambda(iter).slice(k) * design_var.t() * 
        arma::diagmat(eta.slice(iter).col(k));
    }
    current_fit = basis * current_lower_dim_fit;
    stats(arma::vectorise(current_fit));
  }
  
  arma::mat subj_mean = arma::reshape(stats.mean(), basis.n_rows, num_subjects);
  arma::mat subj_sd = arma::reshape(stats.stddev(), basis.n_rows, num_subjects);
  arma::uword counter = 0;
  for (arma::uword iter = burnin; iter < iterations; iter++) {
    current_lower_dim_fit = beta.slice(iter) * design_mean.t();
    for (arma::uword k = 0; k < kdim; k++) {
      current_lower_dim_fit = current_lower_dim_fit + 
        lambda(iter).slice(k) * design_var.t() * 
        arma::diagmat(eta.slice(iter).col(k));
    }
    current_fit = basis * current_lower_dim_fit;
    m_alpha.col(counter) =
      arma::max(arma::abs(current_fit - subj_mean) /
        subj_sd, 0).t();
    counter++;
  }
  arma::vec alpha_vec = {1-alpha};
  arma::vec m_star = arma::quantile(m_alpha, alpha_vec, 1);
  arma::mat mean_overall = subj_mean;
  arma::mat lower = mean_overall - subj_sd * arma::diagmat(m_star);
  arma::mat upper = mean_overall + subj_sd * arma::diagmat(m_star);
  Rcpp::List posterior_bands = 
    Rcpp::List::create(Rcpp::Named("subject_means", mean_overall),
                       Rcpp::Named("subject_lower", lower),
                       Rcpp::Named("subject_upper", upper),
                       Rcpp::Named("current_lower_dim_fit", current_lower_dim_fit));
  return(posterior_bands);
}

// [[Rcpp::export]]
arma::mat get_posterior_means_cpp_correct(List mcmc_results, arma::vec xi, double alpha){
  Rcpp::List data = mcmc_results["data"];
  Rcpp::List samples = mcmc_results["samples"];
  Rcpp::List control = mcmc_results["control"];
  arma::cube beta = samples["beta"];
  arma::mat basis = data["basis"];
  arma::uword iterations = control["iterations"];
  arma::uword burnin = control["burnin"];
  arma::uword basis_dim = basis.n_cols;
  arma::vec m_alpha(iterations);
  arma::running_stat_vec<arma::vec> stats;
  arma::vec m_beta(basis_dim);
  arma::uword num_time_pts = basis.n_rows;
  arma::mat quants(num_time_pts, 3);
  arma::vec current_mean;
  for(arma::uword i = 0; i < iterations; i++){
    current_mean = basis * beta.slice(i) * xi;
    stats(current_mean);
  }
  arma::uword counter = 0;
  for(arma::uword i = burnin; i < iterations; i++){
    current_mean = basis * beta.slice(i) * xi;
    m_alpha(counter) = 
      arma::max(arma::abs(current_mean - stats.mean()) / stats.stddev());
    counter++;
  }
  
  arma::vec alpha_vec = {1 - alpha};
  double m_star = arma::as_scalar(arma::quantile(m_alpha, alpha_vec));
  
  quants.col(0) = stats.mean() - m_star * stats.stddev();
  quants.col(1) = stats.mean();
  quants.col(2) = stats.mean() + m_star * stats.stddev();
  return(quants);
}

// [[Rcpp::export]]
List get_posterior_coefs(List mod, double alpha) {
  arma::field<arma::cube> BetaF = mod["Beta"];
  arma::mat B = mod["B"];
  arma::uword D1 = arma::size(BetaF(0,0,0))(1);
  arma::uword nchains = arma::size(BetaF)(0);
  arma::uword iter = BetaF(0).n_slices;
  arma::vec alpha_vec = {1 - alpha};
  arma::vec Malpha(iter * nchains);
  arma::mat lower(B.n_rows, D1);
  arma::mat mean(B.n_rows, D1);
  arma::mat upper(B.n_rows, D1);
  for (arma::uword d = 0; d < D1; d++) {
    arma::running_stat_vec<arma::vec> stats;
    for (arma::uword u = 0; u < nchains; u++) {
      for(arma::uword i = 0; i < iter; i++) {
        stats(BetaF(u, 0, 0).slice(i).col(d));
      }
    }
    
    for (arma::uword u = 0; u < nchains; u++) {
      for (arma::uword i = 0; i < iter; i++) {
        Malpha(u * iter + i) = arma::max(arma::abs(BetaF(u, 0, 0).slice(i).col(d) - stats.mean()) / stats.stddev());
      }
    }
    
    double q_alpha = arma::as_scalar(arma::quantile(Malpha, alpha_vec));
    lower.col(d) = B * stats.mean() - q_alpha * B * stats.stddev();
    mean.col(d) = B * stats.mean();
    upper.col(d) = B * stats.mean() + q_alpha * B * stats.stddev();
    
  }
  return(List::create(Named("lower", lower),
                      Named("mean", mean),
                      Named("upper", upper)));
}


List extract_eigenfn(arma::cube& Lambda,
                     arma::mat& Psi, arma::mat& Psi_sqrt,
                     arma::mat& Psi_sqrt_inv, arma::mat& B,
                     arma::uword eigenvals, arma::vec z,
                     arma::vec time){
  arma::uword dim_latent = Lambda.n_rows;
  arma::uword dim_spline = B.n_rows;
  arma::mat eigenfn_latent(dim_latent, dim_latent);
  arma::mat eigenfn_spline(dim_spline, dim_latent);
  arma::mat eigenfn_spline_ordered(dim_spline, eigenvals);
  arma::vec eigenval_latent(dim_latent);
  arma::mat eigenfn_latent_ordered(dim_latent, eigenvals);
  arma::vec eigenval_spline(eigenvals);
  arma::vec eigenval_pve(eigenvals);
  arma::mat cov_latent = arma::zeros<arma::mat>(dim_latent, dim_latent);
  arma::mat cov_spline;
  for(arma::uword k = 0; k < Lambda.n_slices; k++){
    cov_latent = cov_latent + Lambda.slice(k) * z * z.t() * Lambda.slice(k).t();
  }
  arma::mat cov_latent_transformed = Psi_sqrt * cov_latent * Psi_sqrt;
  arma::eig_sym(eigenval_latent, eigenfn_latent, cov_latent_transformed);
  eigenfn_spline =  B * Psi_sqrt_inv * eigenfn_latent;
  cov_spline = B * cov_latent * B.t();
  for(arma::uword v = 0; v < eigenvals; v++){
    eigenfn_spline_ordered.col(v) = eigenfn_spline.col(dim_latent - v - 1);
    eigenval_spline(v) = eigenval_latent(dim_latent - v - 1);
    eigenval_pve(v) = eigenval_spline(v) / arma::sum(eigenval_latent);
    
  }
  double magnitude = arma::sum(eigenval_latent);
  return(List::create(Named("eigenfn_spline", eigenfn_spline_ordered),
                      Named("eigenval", eigenval_spline),
                      Named("eigenval_pve", eigenval_pve),
                      Named("cov_spline", cov_spline),
                      Named("magnitude", magnitude)));
}


// [[Rcpp::export]]
arma::mat arma_cov2cor(arma::mat V){
  arma::mat cor(V.n_rows, V.n_rows);
  cor.diag() = arma::ones<arma::vec>(V.n_rows);
  for(arma::uword i = 0; i < V.n_rows - 1; i++){
    for(arma::uword j = i + 1; j < V.n_rows; j++){
      cor(j, i) = V(j, i) * std::pow(V(i, i) * V(j,j), -.5);
    }
  }
  return(arma::symmatl(cor));
}

/*
// [[Rcpp::export]]
List get_posterior_eigen_cpp(Rcpp::List mcmc_results,
                          arma::uword eigenvals,
                          arma::vec zi, double alpha = 0.05){
  Rcpp::List control = mcmc_results["control"];
  Rcpp::List data = mcmc_results["data"];
  Rcpp::List samples = mcmc_results["samples"];
  arma::mat basis = data["basis"];
  arma::field<arma::cube> lambda = samples["lambda"];
  arma::cube beta = samples["beta"];
  arma::vec time = data["time"];
  arma::uword iterations = control["iterations"];
  arma::uword burnin = control["burnin"];
  arma::uword num_post_iter = iterations - burnin;
  arma::uword basis_dim = basis.n_cols;
  arma::mat m_alpha(num_post_iter, eigenvals);
  arma::running_stat_vec<arma::vec> stats_vec;
  arma::running_stat_vec<arma::vec> stats_val;
  arma::running_stat_vec<arma::vec> stats_cov;
  arma::mat psi(basis_dim, basis_dim);
  arma::mat psi_sqrt;
  arma::mat psi_sqrt_inv;
  for(arma::uword j = 0; j < basis_dim; j++){
    for(arma::uword i = 0; i < basis_dim; i++){
      psi(i, j) = 
        arma::as_scalar(arma::trapz(time, basis.col(i) % basis.col(j)));
    }
  }
  psi_sqrt = arma::sqrtmat_sympd(psi);
  psi_sqrt_inv = arma::inv_sympd(psi_sqrt);
  
  arma::mat temp_evec;
  arma::mat eval_mat(num_post_iter, eigenvals);
  arma::mat eval_pve_mat(num_post_iter, eigenvals);
  List eigen_list;
  arma::vec magnitude(num_post_iter);
  arma::mat tempcov;
  arma::uword idx1, idx2;
  arma::uword counter = 0;
  for(arma::uword i = burnin; i < iterations; i++){
    eigen_list = extract_eigenfn(lambda(i),
                                  psi,
                                  psi_sqrt,
                                  psi_sqrt_inv,
                                  basis,
                                  eigenvals,
                                  zi,
                                  time);
    temp_evec = Rcpp::as<arma::mat>(eigen_list["eigenfn_latent"]);
    eval_mat.row(counter) = Rcpp::as<arma::rowvec>(eigen_list["eigenval"]);
    eval_pve_mat.row(counter) = Rcpp::as<arma::rowvec>(eigen_list["eigenval_pve"]);
    stats_cov(arma::vectorise(as<arma::vec>(eigen_list["cov_latent"])));
    magnitude(counter) = eigen_list["magnitude"];
    if (counter == 0) stats_vec(arma::vectorise(temp_evec));
    else {
      // align eigenvectors
      
      for (arma::uword k = 0; k < eigenvals; k++) {
        idx1 = k * basis_dim;
        idx2 = (k + 1) * basis_dim - 1;
        if (arma::sum(arma::square(temp_evec.col(k) + stats_vec.mean().subvec(idx1, idx2))) <
          arma::sum(arma::square(temp_evec.col(k) - stats_vec.mean().subvec(idx1, idx2)))) {
          temp_evec.col(k) = -temp_evec.col(k);
        }
      }
      Rcpp::Rcout << "size: " << arma::size(temp_evec) << "\n";
      stats_vec(arma::vectorise(temp_evec));
    }
    counter++;
  }
  counter = 0;
  for (arma::uword i = burnin; i < iterations; i++) {
    eigen_list = extract_eigenfn(lambda(i),
                                  psi,
                                  psi_sqrt,
                                  psi_sqrt_inv,
                                  basis,
                                  eigenvals,
                                  zi,
                                  time);
    temp_evec = Rcpp::as<arma::mat>(eigen_list["eigenfn_latent"]);
    
    for (arma::uword k = 0; k < eigenvals; k++) {
      idx1 = k * basis_dim;
      idx2 = (k + 1) * basis_dim - 1;
      if (arma::sum(arma::square(temp_evec.col(k) + stats_vec.mean().subvec(idx1, idx2))) <
        arma::sum(arma::square(temp_evec.col(k) - stats_vec.mean().subvec(idx1, idx2)))) {
        temp_evec.col(k) = -temp_evec.col(k);
      }
      m_alpha(counter, k) = arma::max((temp_evec.col(k) -
        stats_vec.mean().subvec(idx1, idx2)) / stats_vec.stddev().subvec(idx1, idx2));
    }
    counter++;
  }
  
  arma::vec alpha_vec = {1 - alpha};
  arma::vec alpha_vec_eval = {alpha / 2.0, 0.5, 1 - alpha / 2.0};
  double q_alpha;
  arma::mat lower(time.n_elem, eigenvals);
  arma::mat mean(time.n_elem, eigenvals);
  arma::mat upper(time.n_elem, eigenvals);
  arma::mat eigenbands(time.n_elem, eigenvals * 3);
  arma::mat eigenval_intervals(3, eigenvals);
  arma::mat eigenval_pve_intervals(3, eigenvals);
  arma::vec magnitude_interval(3);
  for (arma::uword k = 0; k < eigenvals; k++) {
    idx1 = k * basis_dim;
    idx2 = (k + 1) * basis_dim - 1;
    q_alpha = arma::as_scalar(arma::quantile(m_alpha.col(k), alpha_vec));
    lower.col(k) = basis * (stats_vec.mean().subvec(idx1, idx2) -
      q_alpha * stats_vec.stddev().subvec(idx1, idx2));
    mean.col(k) = basis * (stats_vec.mean().subvec(idx1, idx2));
    upper.col(k) = basis * (stats_vec.mean().subvec(idx1, idx2) + 
      q_alpha * stats_vec.stddev().subvec(idx1, idx2));
    eigenval_intervals.col(k) = arma::quantile(eval_mat.col(k), alpha_vec_eval);
    eigenval_pve_intervals.col(k) = arma::quantile(eval_pve_mat.col(k), alpha_vec_eval);
  }
  magnitude_interval = arma::quantile(magnitude, alpha_vec_eval);
  arma::mat mean_cov = basis * 
    arma::reshape(stats_cov.mean(), basis_dim, basis_dim) * 
    basis.t();
  return(List::create(Named("lower_eigen", lower),
                      Named("mean_eigen", mean),
                      Named("upper_eigen", upper),
                      Named("eigenval_intervals", eigenval_intervals),
                      Named("eigenval_pve_intervals", eigenval_pve_intervals),
                      Named("surface", mean_cov),
                      Named("magnitude", magnitude_interval),
                      Named("raw_magnitude", magnitude),
                      Named("time", time)));
}

*/
// [[Rcpp::export]]
List get_posterior_eigen_cpp_correct(Rcpp::List mcmc_results,
                             arma::uword eigenvals,
                             arma::vec zi, double alpha = 0.05){
  Rcpp::List control = mcmc_results["control"];
  Rcpp::List data = mcmc_results["data"];
  Rcpp::List samples = mcmc_results["samples"];
  arma::mat basis = data["basis"];
  arma::field<arma::cube> lambda = samples["lambda"];
  arma::cube beta = samples["beta"];
  arma::vec time = data["time"];
  arma::uword iterations = control["iterations"];
  arma::uword burnin = control["burnin"];
  arma::uword num_post_iter = iterations - burnin;
  arma::uword basis_dim = basis.n_cols;
  arma::mat m_alpha(num_post_iter, eigenvals);
  arma::running_stat_vec<arma::vec> stats_vec;
  arma::running_stat_vec<arma::vec> stats_val;
  arma::running_stat_vec<arma::vec> stats_cov;
  arma::mat psi(basis_dim, basis_dim);
  arma::mat psi_sqrt;
  arma::mat psi_sqrt_inv;
  for(arma::uword j = 0; j < basis_dim; j++){
    for(arma::uword i = 0; i < basis_dim; i++){
      psi(i, j) = 
        arma::as_scalar(arma::trapz(time, basis.col(i) % basis.col(j)));
    }
  }
  psi_sqrt = arma::sqrtmat_sympd(psi);
  psi_sqrt_inv = arma::inv_sympd(psi_sqrt);
  
  arma::mat temp_evec;
  arma::mat eval_mat(num_post_iter, eigenvals);
  arma::mat eval_pve_mat(num_post_iter, eigenvals);
  List eigen_list;
  arma::vec magnitude(num_post_iter);
  arma::mat tempcov;
  arma::uword idx1, idx2;
  arma::uword counter = 0;
  for(arma::uword i = burnin; i < iterations; i++){
    eigen_list = extract_eigenfn(lambda(i),
                                 psi,
                                 psi_sqrt,
                                 psi_sqrt_inv,
                                 basis,
                                 eigenvals,
                                 zi,
                                 time);
    
    temp_evec = Rcpp::as<arma::mat>(eigen_list["eigenfn_spline"]);
    eval_mat.row(counter) = Rcpp::as<arma::rowvec>(eigen_list["eigenval"]);
    eval_pve_mat.row(counter) = Rcpp::as<arma::rowvec>(eigen_list["eigenval_pve"]);
    stats_cov(arma::vectorise(as<arma::vec>(eigen_list["cov_spline"])));
    magnitude(counter) = eigen_list["magnitude"];
    if (counter == 0) stats_vec(arma::vectorise(temp_evec));
    else {
      // align eigenvectors
      
      for (arma::uword k = 0; k < eigenvals; k++) {
        idx1 = k * basis.n_rows;
        idx2 = (k + 1) * basis.n_rows - 1;
        if (arma::sum(arma::square(temp_evec.col(k) + stats_vec.mean().subvec(idx1, idx2))) <
          arma::sum(arma::square(temp_evec.col(k) - stats_vec.mean().subvec(idx1, idx2)))) {
          temp_evec.col(k) = -temp_evec.col(k);
        }
      }
      
      stats_vec(arma::vectorise(temp_evec));
    }
    counter++;
  }
  counter = 0;
  for (arma::uword i = burnin; i < iterations; i++) {
    eigen_list = extract_eigenfn(lambda(i),
                                 psi,
                                 psi_sqrt,
                                 psi_sqrt_inv,
                                 basis,
                                 eigenvals,
                                 zi,
                                 time);
    temp_evec = Rcpp::as<arma::mat>(eigen_list["eigenfn_spline"]);
    
    for (arma::uword k = 0; k < eigenvals; k++) {
      idx1 = k * basis.n_rows;
      idx2 = (k + 1) * basis.n_rows - 1;
      if (arma::sum(arma::square(temp_evec.col(k) + stats_vec.mean().subvec(idx1, idx2))) <
        arma::sum(arma::square(temp_evec.col(k) - stats_vec.mean().subvec(idx1, idx2)))) {
        temp_evec.col(k) = -temp_evec.col(k);
      }
      m_alpha(counter, k) = arma::max((temp_evec.col(k) -
        stats_vec.mean().subvec(idx1, idx2)) / stats_vec.stddev().subvec(idx1, idx2));
    }
    counter++;
  }
  
  arma::vec alpha_vec = {1 - alpha};
  arma::vec alpha_vec_eval = {alpha / 2.0, 0.5, 1 - alpha / 2.0};
  double q_alpha;
  arma::mat lower(time.n_elem, eigenvals);
  arma::mat mean(time.n_elem, eigenvals);
  arma::mat upper(time.n_elem, eigenvals);
  arma::mat eigenbands(time.n_elem, eigenvals * 3);
  arma::mat eigenval_intervals(3, eigenvals);
  arma::mat eigenval_pve_intervals(3, eigenvals);
  arma::vec magnitude_interval(3);
  for (arma::uword k = 0; k < eigenvals; k++) {
    idx1 = k * basis.n_rows;
    idx2 = (k + 1) * basis.n_rows - 1;
    q_alpha = arma::as_scalar(arma::quantile(m_alpha.col(k), alpha_vec));
    lower.col(k) = (stats_vec.mean().subvec(idx1, idx2) -
      q_alpha * stats_vec.stddev().subvec(idx1, idx2));
    mean.col(k) = (stats_vec.mean().subvec(idx1, idx2));
    upper.col(k) = (stats_vec.mean().subvec(idx1, idx2) + 
      q_alpha * stats_vec.stddev().subvec(idx1, idx2));
    eigenval_intervals.col(k) = arma::quantile(eval_mat.col(k), alpha_vec_eval);
    eigenval_pve_intervals.col(k) = arma::quantile(eval_pve_mat.col(k), alpha_vec_eval);
  }
  magnitude_interval = arma::quantile(magnitude, alpha_vec_eval);
  arma::mat mean_cov =
    arma::reshape(stats_cov.mean(), basis.n_rows, basis.n_rows);
  return(List::create(Named("lower_eigen", lower),
                      Named("mean_eigen", mean),
                      Named("upper_eigen", upper),
                      Named("eigenval_intervals", eigenval_intervals),
                      Named("eigenval_pve_intervals", eigenval_pve_intervals),
                      Named("surface", mean_cov),
                      Named("magnitude", magnitude_interval),
                      Named("raw_magnitude", magnitude),
                      Named("time", time)));
}

// [[Rcpp::export]]
List get_variance_effects(List mod, double alpha){
  arma::mat B = mod["B"];
  arma::field<arma::cube> LambdaF = mod["Lambda"];
  arma::vec Time = mod["Time"];
  arma::mat Z = mod["Z"];
  arma::vec alpha_vec = {1 - alpha};
  arma::mat ph = arma::zeros<arma::mat>(B.n_cols, B.n_cols);
  arma::uword iter = arma::size(LambdaF)(0);
  arma::uword K = LambdaF(0, 0, 0).n_slices;
  arma::mat mean(Time.n_elem, Z.n_cols);
  arma::mat upper(Time.n_elem, Z.n_cols);
  arma::mat lower(Time.n_elem, Z.n_cols);
  arma::mat temp_lambda;
  arma::mat temp_shed;
  arma::uvec keep_cols = arma::linspace<arma::uvec>(1, Z.n_cols - 1);
  arma::vec Malpha(iter);
  for(arma::uword d = 0; d < Z.n_cols; d++){
    arma::running_stat_vec<arma::vec> stats;
    for (arma::uword i = 0; i < iter; i++) {
      for (arma::uword k = 0; k < K; k++) {
        temp_lambda = LambdaF(i, 0, 0).slice(k);
        for (arma::uword dp = 0; dp < Z.n_cols; dp++){
          if (dp == d) {
            ph = ph + temp_lambda.col(dp) * temp_lambda.col(d).t();
          }
          /*else {
            ph = ph + temp_lambda.col(d) * temp_lambda.col(dp).t() +
              temp_lambda.col(dp) * temp_lambda.col(d).t();          
          }*/
        }
      }
      stats(arma::vectorise(ph));
      ph.zeros();
    }
    for (arma::uword i = 0; i < iter; i++) {
      for (arma::uword k = 0; k < K; k++) {
        temp_lambda = LambdaF(i, 0, 0).slice(k);
        for (arma::uword dp = 0; dp < Z.n_cols; dp++){
          if (dp == d) {
            ph = ph + temp_lambda.col(dp) * temp_lambda.col(d).t();
          }
          /*else {
            ph = ph + temp_lambda.col(d) * temp_lambda.col(dp).t() +
              temp_lambda.col(dp) * temp_lambda.col(d).t();          
          }*/
        }
      }
      Malpha(i) = arma::max(arma::abs((arma::vectorise(ph) - stats.mean()) / stats.stddev()));
      ph.zeros();
    }
    double q_alpha = arma::as_scalar(arma::quantile(Malpha, alpha_vec));
    Rcout << q_alpha << std::endl;
    lower.col(d) = arma::diagvec(B * arma::reshape((stats.mean() - q_alpha * stats.stddev()), B.n_cols, B.n_cols) * B.t());
    mean.col(d) = arma::diagvec(B * arma::reshape(stats.mean(), B.n_cols, B.n_cols) * B.t());
    upper.col(d) = arma::diagvec(B * arma::reshape((stats.mean() + q_alpha * stats.stddev()), B.n_cols, B.n_cols) * B.t());
    
    if (d != Z.n_cols - 1) {
      keep_cols(d) = keep_cols(d) - 1;
    }
    
  }
  return(List::create(Named("mean", mean),
                      Named("lower", lower),
                      Named("upper", upper)));
}

