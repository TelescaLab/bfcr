#include "Sampler.h"

// [[Rcpp::depends(RcppArmadillo)]]

Rcpp::List Sampler::write_data() {
  return(dat.write_data());
}
Rcpp::List Sampler::write_control() {
  return(Rcpp::List::create(
      Rcpp::Named("iterations", dat.iter),
      Rcpp::Named("thin", dat.thin),
      Rcpp::Named("burnin", dat.burnin)
  ));
}

void SamplerUnequal::sample_parameters() {
  Progress progress_bar(dat.iter, true);
  for (arma::uword i = 0; i < dat.iter; i++) {
    for (arma::uword j = 0; j < dat.thin; j++) {
      if (Progress::check_abort() ){
        Rcpp::Rcout << "MCMC aborted" << std::endl;
        
        goto stop;
      }
      pars.update_beta(dat, transf);
      pars.update_lambda(dat, transf); 
      pars.update_eta(dat, transf);
      pars.update_tau1(dat, transf);
      pars.update_tau2(dat, transf);
      pars.update_phi(dat, transf);
      pars.update_delta(dat, transf);
      pars.update_psi(dat, transf);
      pars.update_tausq(dat, transf);
      pars.varphi = pars.psi;
      transf.complete_response(dat, pars);
    }
    progress_bar.increment();
    write_parameters();
  }
  stop:
    NULL;
}

void SamplerPooled::sample_parameters() {
  Progress progress_bar(dat.iter, true);
  for (arma::uword i = 0; i < dat.iter; i++) {
    for (arma::uword j = 0; j < dat.thin; j++) {
      if (Progress::check_abort() ){
        Rcpp::Rcout << "MCMC aborted" << std::endl;
        
        goto stop;
      }
      pars.update_beta(dat, transf);
      pars.update_lambda(dat, transf); 
      pars.update_eta(dat, transf);
      pars.update_tau1(dat, transf);
      pars.update_tau2(dat, transf);
      pars.update_phi(dat, transf);
      pars.update_delta(dat, transf);
      pars.update_varphi(dat, transf); 
      transf.complete_response(dat, pars);
    }
    progress_bar.increment();
    pars.write_parameters();
  }
  stop:
    NULL;
}

void SamplerUnequal::write_parameters() {
  pars.beta_container.slice(pars.iteration) = pars.beta;
  pars.lambda_container(pars.iteration) = pars.lambda;
  pars.eta_container.slice(pars.iteration) = pars.eta;
  pars.varphi_container.col(pars.iteration) = pars.varphi;
  pars.tau1_container.col(pars.iteration) = pars.tau1;
  pars.tau2_container.slice(pars.iteration) = pars.tau2;
  pars.delta_container.slice(pars.iteration) = pars.delta;
  pars.tausq_container(pars.iteration) = pars.tausq;
  pars.iteration++;
}

void SamplerPooled::write_parameters() {
  pars.beta_container.slice(pars.iteration) = pars.beta;
  pars.lambda_container(pars.iteration) = pars.lambda;
  pars.eta_container.slice(pars.iteration) = pars.eta;
  pars.varphi_container.col(pars.iteration) = pars.varphi;
  pars.tau1_container.col(pars.iteration) = pars.tau1;
  pars.tau2_container.slice(pars.iteration) = pars.tau2;
  pars.delta_container.slice(pars.iteration) = pars.delta;
  pars.iteration++;
}

Rcpp::List SamplerUnequal::get_samples() {
  return(Rcpp::List::create(Rcpp::Named("lambda", pars.lambda_container),
                            Rcpp::Named("beta", pars.beta_container),
                            Rcpp::Named("eta", pars.eta_container),
                            Rcpp::Named("varphi", pars.varphi_container),
                            Rcpp::Named("tau1", pars.tau1_container),
                            Rcpp::Named("tau2", pars.tau2_container),
                            Rcpp::Named("phi", pars.phi_container),
                            Rcpp::Named("delta", pars.delta_container),
                            Rcpp::Named("tausq", pars.tausq_container)));
}

Rcpp::List SamplerPooled::get_samples() {
  return(Rcpp::List::create(Rcpp::Named("lambda", pars.lambda_container),
                            Rcpp::Named("beta", pars.beta_container),
                            Rcpp::Named("eta", pars.eta_container),
                            Rcpp::Named("varphi", pars.varphi_container),
                            Rcpp::Named("tau1", pars.tau1_container),
                            Rcpp::Named("tau2", pars.tau2_container),
                            Rcpp::Named("phi", pars.phi_container),
                            Rcpp::Named("delta", pars.delta_container)));
}
