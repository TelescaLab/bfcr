#include <RcppArmadillo.h>
List MCMC_Wrapper(arma::mat Y, arma::mat X, arma::mat B, arma::uword K, arma::uword iter, arma::uword nchains, arma::uword thin, arma::mat Theta_init, arma::cube Lambda_init, arma::mat Eta_init);
