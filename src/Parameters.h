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
  void update_phi(Data&, Transformations&);
  void update_delta(Data&, Transformations&);
  void update_psi(Data&, Transformations&);
  void update_tausq(Data&, Transformations&);
  void update_alpha(Data&, Transformations&);
  void update_a1(Data&);
  void update_a2(Data&);
  void write_parameters();
  double varphi_a = .001;
  double varphi_b = .001;
  double tau_a = 1;
  double tau_b = 0.005;
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
};
#endif