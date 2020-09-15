#ifndef PARAMETERS_H
#define PARAMETERS_H

#include <RcppArmadillo.h>
#include "Data.h"
#include "Transformations.h"
#include "Utils.h"
class Transformations;

class Parameters {
public:
  arma::cube lambda;
  arma::mat beta;
  arma::mat eta;
  arma::vec varphi;
  arma::cube phi;
  double alpha;
  double tausq;
  arma::vec psi;
  arma::vec tau1;
  arma::mat tau2;
  arma::vec tau1_delta;
  arma::mat tau2_delta;
  arma::vec tau1_nu;
  arma::mat tau2_nu;
  arma::mat delta;
  arma::uword iteration = 0;
  Parameters(Data&);
  Parameters() {};
  void update_beta(Data&, Transformations&);
  void update_lambda(Data&, Transformations&);
  void update_eta(Data&, Transformations&);
  void update_varphi(Data&, Transformations&);
  void update_tau1(Data&, Transformations&);
  void update_tau2(Data&, Transformations&);
  void update_tau1_delta(Data&, Transformations&);
  void update_tau2_delta(Data&, Transformations&);
  void update_tau1_nu(Data&, Transformations&);
  void update_tau2_nu(Data&, Transformations&);
  void update_phi(Data&, Transformations&);
  void update_delta(Data&, Transformations&);
  void update_psi(Data&, Transformations&);
  void update_tausq(Data&, Transformations&);
  void update_alpha(Data&, Transformations&);
  void update_a1(Data&);
  void update_a2(Data&);
  void write_parameters();
  double tau_a;
  double tau_b;
  double tau_cutoff;
  double varphi_a;
  double varphi_b;
  double phi_a = 2;
  double phi_b = 2;
  double delta_a1 = 2;
  double delta_a2 = 2;
  double delta_b = 1;
  double nu = 5;
  double a1, a2;
  arma::vec a1_, a2_;
  double a1_a = 2, a2_a = 2;
  arma::mat a1_container, a2_container;
  arma::cube beta_container;
  arma::field<arma::cube> lambda_container;
  arma::cube eta_container;
  arma::mat varphi_container;
  arma::mat tau1_container;
  arma::cube tau2_container;
  arma::field<arma::cube> phi_container;
  arma::cube delta_container;
  arma::vec tausq_container;
  arma::mat tau1_delta_container;
  arma::mat tau1_nu_container;
  arma::cube tau2_delta_container;
  arma::cube tau2_nu_container;
};
#endif