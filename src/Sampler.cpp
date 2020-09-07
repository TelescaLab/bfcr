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
      //pars.update_phi(dat, transf);
      //pars.update_delta(dat, transf);
      pars.update_psi(dat, transf);
      pars.update_tausq(dat, transf);
      pars.update_a1(dat);
      pars.update_a2(dat);
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
      pars.update_a1(dat);
      pars.update_a2(dat);
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
  pars.phi_container(pars.iteration) = pars.phi;
  pars.delta_container.slice(pars.iteration) = pars.delta;
  pars.tausq_container(pars.iteration) = pars.tausq;
  pars.a1_container.col(pars.iteration) = pars.a1_;
  pars.a2_container.col(pars.iteration) = pars.a2_;
  pars.iteration++;
}

void SamplerPooled::write_parameters() {
  pars.beta_container.slice(pars.iteration) = pars.beta;
  pars.lambda_container(pars.iteration) = pars.lambda;
  pars.eta_container.slice(pars.iteration) = pars.eta;
  pars.varphi_container.col(pars.iteration) = pars.varphi;
  pars.tau1_container.col(pars.iteration) = pars.tau1;
  pars.tau2_container.slice(pars.iteration) = pars.tau2;
  pars.phi_container(pars.iteration) = pars.phi;
  pars.delta_container.slice(pars.iteration) = pars.delta;
  pars.a1_container.col(pars.iteration) = pars.a1_;
  pars.a2_container.col(pars.iteration) = pars.a2_;
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
                            Rcpp::Named("tausq", pars.tausq_container),
                            Rcpp::Named("a1", pars.a1_container),
                            Rcpp::Named("a2", pars.a2_container),
                            Rcpp::Named("block_diag_mean", transf.blk_diag_mean_penalties),
                            Rcpp::Named("blok_diag_var", transf.blk_diag_var_penalties),
                            Rcpp::Named("block_diag_phi_delta", transf.blk_diag_phi_delta)));
}

Rcpp::List SamplerPooled::get_samples() {
  return(Rcpp::List::create(Rcpp::Named("lambda", pars.lambda_container),
                            Rcpp::Named("beta", pars.beta_container),
                            Rcpp::Named("eta", pars.eta_container),
                            Rcpp::Named("varphi", pars.varphi_container),
                            Rcpp::Named("tau1", pars.tau1_container),
                            Rcpp::Named("tau2", pars.tau2_container),
                            Rcpp::Named("phi", pars.phi_container),
                            Rcpp::Named("delta", pars.delta_container),
                            Rcpp::Named("a1", pars.a1_container),
                            Rcpp::Named("a2", pars.a2_container)));
}
