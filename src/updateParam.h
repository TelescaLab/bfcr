#ifndef MODSTRING_H
#define MODSTRING_H
#include <RcppArmadillo.h>
#include "Utility.h"

//void updateLambda2(arma::mat& Y, arma::cube& Lambda, arma::mat& Tau, arma::mat& c, arma::mat& Gamma, arma::mat& X, arma::mat& B, double prec, arma::mat& Theta);
void updateLambda2(arma::mat& Y, arma::cube& Lambda, arma::mat& Tau, arma::mat& Gamma, arma::mat& X, arma::mat& B, double prec, arma::mat& Theta);
void updateLambdaS(arma::field<arma::vec>& Y, arma::cube& Lambda, arma::mat& Tau, arma::mat& c, arma::mat& Gamma, arma::mat& X, arma::field<arma::mat>& B, double prec, arma::mat& Theta);
void updateLambda3(arma::mat& Y, arma::cube& Lambda, arma::mat& Tau, arma::mat& Eta, arma::mat& X, arma::mat& B, double prec, arma::mat& Theta);
void updateTheta(arma::mat& Y, arma::cube& Lambda, arma::mat& Tau, arma::mat& Gamma, arma::mat& X, arma::mat& B, double prec, arma::mat& Theta);
void updateThetaS(arma::field<arma::vec>& Y, arma::cube& Lambda, arma::mat& Tau, arma::mat& Gamma,
                  arma::mat& X, arma::field<arma::mat> B, double prec, arma::mat& Theta);
void updateEta(arma::mat& Y, arma::cube& Lambda, arma::vec& Sigma, arma::mat& Eta, arma::mat& X, arma::mat& B, double prec, arma::mat& Theta);
void updateEta2(arma::mat& Y, arma::cube& Lambda, arma::vec& Sigma, arma::mat& Eta, arma::mat& X, arma::mat& B, double prec, arma::mat& Theta);
void updateEta3(arma::mat& Y, arma::cube& Lambda, arma::mat& Eta, arma::mat& X, arma::mat& B, double prec, arma::mat& Theta);
void updateEtaS(arma::field<arma::vec>& Y, arma::cube& Lambda, arma::vec& Sigma, arma::mat& Eta, arma::mat& X, arma::field<arma::mat>& B, double prec, arma::mat& Theta);
double updatePrec(arma::mat& Y, arma::cube& Lambda, arma::mat Gamma, arma::mat& X, arma::mat& B, arma::mat& Theta);
double updatePrecS(arma::field<arma::vec>& Y, arma::cube& Lambda, arma::mat Gamma, arma::mat& X, arma::field<arma::mat>& B, arma::mat& Theta);
void updateTau(arma::mat& Theta, arma::cube& Lambda, arma::mat& Tau);
void updateSigma(arma::mat& Eta, arma::vec& Sigma);
void updatec(arma::cube& Lambda, arma::mat& c);
void updateTheta2(arma::mat& Y, arma::cube& Lambda, arma::mat& Tau, arma::mat& X, arma::mat& B, double prec, arma::mat& Theta);
void updateThetaLambda(arma::mat &Y, arma::cube& Lambda, arma::mat& Eta, arma::mat& Tau, arma::mat& X, arma::mat& B, double prec, arma::mat& Theta);
void updateThetaLambdaMH(arma::mat& Y, arma::mat& Theta, arma::cube& Lambda, arma::mat& Tau, double prec, arma::mat& X, arma::mat& B, double noise, arma::uword n);

void updateProj(arma::cube& Lambda, arma::mat& Theta, arma::mat& Eta, arma::vec& Delta, double Prec, arma::mat& X, arma::mat& Y, arma::mat B, arma::mat& Proj);
void updateThetaLambdaP(arma::cube& Lambda, arma::mat& Theta, arma::mat& Eta, arma::vec& Delta, arma::mat& Proj, arma::mat& Tau, arma::mat& X);
void updateEtaP(arma::cube& Lambda, arma::mat& Theta, arma::mat& Eta, arma::vec& Delta, arma::mat& Proj, arma::mat& X);
void updateDelta(arma::mat& Proj, arma::mat& Theta, arma::cube& Lambda, arma::mat& Eta, arma::vec& Delta, arma::mat& X);
double updatePrecP(arma::mat& Proj, arma::mat& Y, arma::mat& B);
#endif
