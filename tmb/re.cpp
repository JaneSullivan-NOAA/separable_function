// A state-space random walk model that accepts a biomass and CPUE survey index.
// Used for survey averaging

#include <TMB.hpp>
#include <iostream>

#define see(object) std::cout << #object ":\n" << object << "\n";
template <class Type>
Type square(Type x){return x * x;}


template<class Type>
Type objective_function<Type>::operator() ()
{

  // model dimensions
  DATA_IVECTOR(yrs);

  // survey biomass
  DATA_VECTOR(biomass);
  DATA_VECTOR(biomass_cv);

  // survey cpue
  DATA_VECTOR(cpue);
  DATA_VECTOR(cpue_cv);

  // parameter section

  PARAMETER(log_PE); // process errors
  PARAMETER(log_q); // scaling parameters
  PARAMETER_VECTOR(log_biomass_pred); // random effects predicted biomass

  // negative log likelihood
  vector<Type> jnll(3); // random walk, biomass obs, cpue obs
  jnll.setZero();
  Type nll = 0;

  int nyrs = yrs.size();

  // model predictions

  // predicted biomass and cpue on the natural and log scale
  vector<Type> biomass_pred(nyrs);
  biomass_pred.setZero();
  vector<Type> cpue_pred(nyrs);
  cpue_pred.setZero();
  vector<Type> log_cpue_pred(nyrs);
  log_cpue_pred.setZero();

  // derived quantities - biomass survey observations
  vector<Type> log_biomass(nyrs);
  log_biomass = log(biomass);
  vector<Type> log_biomass_sd(nyrs);
  log_biomass_sd = biomass_cv * biomass_cv + Type(1.0);
  log_biomass_sd = sqrt(log(log_biomass_sd));

  // derived quantities - cpue survey observations
  vector<Type> log_cpue(nyrs);
  log_cpue = log(cpue);
  vector<Type> log_cpue_sd(nyrs);
  log_cpue_sd = cpue_cv * cpue_cv + Type(1.0);
  log_cpue_sd = sqrt(log(log_cpue_sd));

  // random effects contribution to likelihood
  Type nll_pe = 0;

  for(int i = 1; i < nyrs; i++) {

    jnll(0) -= dnorm(log_biomass_pred(i-1), log_biomass_pred(i), exp(log_PE), 1);

    nll_pe += Type(0.5) * (log(Type(2.0) * PI * exp(Type(2.0) * log_PE)) +
      square(log_biomass_pred(i-1) - log_biomass_pred(i)) / exp(Type(2.0) * log_PE));

  }

  // likelihood for biomass survey data/observation error
  Type nll_biomass = 0;

  for(int i = 0; i < nyrs; i++) {

    if(biomass(i) > 0) {
      jnll(1) -= dnorm(log_biomass_pred(i), log_biomass(i), log_biomass_sd(i), 1);

      nll_biomass += Type(0.5) * (log(Type(2.0) * PI * square(log_biomass_sd(i))) +
        square(log_biomass_pred(i) - log_biomass(i)) / square(log_biomass_sd(i)));
    }
  }
  biomass_pred = exp(log_biomass_pred);

  // predicted values for cpue survey
  for(int i = 0; i < nyrs; i++) {
    cpue_pred(i) = exp(log_q) * biomass_pred(i);
  }
  log_cpue_pred = log(cpue_pred);

  // likelihood for cpue survey data/observation error - add 0.5 likelihood
  // weight (to get the outside admb version to converge)
  Type nll_cpue = 0;

  for(int i = 0; i < nyrs; i++) {
    if(cpue(i) > 0) {
      jnll(2) -= Type(0.5) * dnorm(log_cpue_pred(i), log_cpue(i), log_cpue_sd(i), 1);

      nll_cpue += Type(0.5) * Type(0.5) * (log(Type(2.0) * PI * square(log_cpue_sd(i))) +
        square(log_cpue_pred(i) - log_cpue(i)) / square(log_cpue_sd(i)));
    }
  }

  // report section
  ADREPORT(log_biomass_pred);
  ADREPORT(log_cpue_pred);

  REPORT(jnll);
  REPORT(nll_pe);
  REPORT(nll_biomass);
  REPORT(nll_cpue);

  nll = jnll.sum();
  return nll;
}
