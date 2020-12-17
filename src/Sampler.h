#ifndef SAMPLER_H
#define SAMPLER_H

#include <RcppArmadillo.h>
#include "Data.h"
#include "Transformations.h"
#include "Parameters.h"
#include <progress.hpp>

class Sampler
{
public:
  Data dat;
  Parameters pars;
  Transformations transf;
  // parameterized constructor
  Sampler(Data& dat, Parameters& pars, Transformations& transf) :
    dat(dat), pars(pars), transf(transf){}
  virtual void sample_parameters() = 0;
  virtual void write_parameters() = 0;
  virtual Rcpp::List get_samples() = 0;
  Rcpp::List write_data();
  Rcpp::List write_control();
  ~Sampler(){}
};

class SamplerPooled : public Sampler
{
public:
  // parameterized constructor
  SamplerPooled(Data& dat, Parameters& pars, Transformations& transf) :
  Sampler(dat, pars, transf)
  {}
  void sample_parameters();
  void write_parameters();
  Rcpp::List get_samples();
  ~SamplerPooled(){};
};

class SamplerUnequal : public Sampler{
public:
  SamplerUnequal(Data& dat, Parameters& pars, Transformations& transf) :
  Sampler(dat, pars, transf) {}
  void sample_parameters();
  void write_parameters();
  Rcpp::List get_samples();
  ~SamplerUnequal(){}
};



class SamplerFactory {
public:
  static Sampler *new_mcmc(std::string type, Data& dat, Parameters& pars, Transformations& transf) {
    if(type == "pooled") return new SamplerPooled(dat, pars, transf);
    if(type == "unequal") return new SamplerUnequal(dat, pars, transf);
    return nullptr;
  }
};

#endif
