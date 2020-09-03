#ifndef MODSTRING_H
#define MODSTRING_H
#include <RcppArmadillo.h>
#include "Utils.h"

void updateLambdaS(arma::field<arma::vec>& Y, arma::cube& Lambda, arma::mat& Tau, arma::mat& c, arma::mat& Gamma, arma::mat& X, arma::field<arma::mat>& B, double prec, arma::mat& Theta);
void updateThetaS(arma::field<arma::vec>& Y, arma::cube& Lambda, arma::mat& Tau, arma::mat& Gamma,
                  arma::mat& X, arma::field<arma::mat> B, double prec, arma::mat& Theta);
void updateEtaS(arma::field<arma::vec>& Y, arma::cube& Lambda, arma::vec& Sigma, arma::mat& Eta, arma::mat& X, arma::field<arma::mat>& B, double prec, arma::mat& Theta);
double updatePrecS(arma::field<arma::vec>& Y, arma::cube& Lambda, arma::mat Gamma, arma::mat& X, arma::field<arma::mat>& B, arma::mat& Theta);
#endif
