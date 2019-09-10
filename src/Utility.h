#ifndef MODSTRING_H1
#define MODSTRING_H1

#include <RcppArmadillo.h>
using namespace Rcpp;

arma::mat getPenalty(arma::uword n);
arma::mat getPenalty2(arma::uword n, arma::uword D);
arma::vec initializeY(arma::vec y, arma::vec observedTimes, arma::vec fullTimes);
arma::mat PredictY(arma::mat& ImputedY, arma::mat X, arma::mat B, arma::mat Theta, arma::mat Eta, arma::cube Lambda, double Prec);
arma::uvec test(arma::field<arma::vec> observedTimes, arma::vec fullTimes);
arma::uvec getObservedOrder(arma::vec observedTimes, arma::vec fullTimes);
void PredictY2(arma::mat& ImputedY, arma::field<arma::uvec> observedOrder, arma::mat X,
               arma::mat B, arma::mat Theta, arma::cube Lambda, arma::mat Eta, double Prec);
double gam_trunc_left(double a, double b,  double cut);
#endif