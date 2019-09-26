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
Rcpp::List Proposal(arma::vec Theta, arma::mat Lambda, double noise = .01, arma::uword samples = 100);
double cpploglik_bayes(arma::mat &Theta, arma::cube &Lambda, double precision,
                       arma::mat &X, arma::mat &B, arma::mat &Y, int cores = 1);
#endif