// Example of SEPARABLE_FUNCTION behavior using a state-space random walk model that accepts a biomass and CPUE survey index. You get a different answer if predicted cpue (pred_cpue = q * biomass) is calcuated inside or outside the separable function (note that q = scaling parameter). The ADMB INSIDE version matches TMB; the ADMB OUTSIDE doesn't. The ADMB OUTSIDE version is what's currently being used. It's pretty sensitive; I had to play around with the likelihood weight to get it to converge.

DATA_SECTION
  // separable function flag (0 = pred cpue calculated with q inside function, 1 = predicted cpue calculated outside function)
  init_int sep_fxn_flag;

  // model start and end yrs
  init_int styr;
  init_int endyr;
  ivector yrs(styr,endyr);
  !! yrs.fill_seqadd(styr,1);

  // define biomass survey data
  init_int n_biomass
  init_ivector yrs_biomass(1,n_biomass);
  init_vector biomass(1,n_biomass);
  init_vector biomass_cv(1,n_biomass);

  vector log_biomass(1,n_biomass);
  !! log_biomass = log(biomass);

  vector log_biomass_sd(1,n_biomass);
  !! log_biomass_sd = elem_prod(biomass_cv,biomass_cv) + 1.0;
  !! log_biomass_sd = sqrt(log(log_biomass_sd));

  // define cpue survey data
  init_int n_cpue;
  init_ivector yrs_cpue(1,n_cpue);
  init_vector cpue(1,n_cpue);
  init_vector cpue_cv(1,n_cpue);

  vector log_cpue(1,n_cpue);
  !! log_cpue = log(cpue);

  vector log_cpue_sd(1,n_cpue);
  !! log_cpue_sd = elem_prod(cpue_cv,cpue_cv) + 1.0;
  !! log_cpue_sd = sqrt(log(log_cpue_sd));

PARAMETER_SECTION

  init_number log_PE;
  init_number log_q;
  random_effects_vector log_biomass_pred(styr,endyr);

  sdreport_vector biomass_pred(styr,endyr);
  sdreport_vector cpue_pred(styr,endyr);
  sdreport_vector log_cpue_pred(styr,endyr);

  // separate likelihood components
  number nll_pe;
  number nll_biomass;
  number nll_cpue;
  objective_function_value jnll;

PROCEDURE_SECTION

  jnll = 0.0;
  nll_pe = 0.0;
  nll_biomass = 0.0;
  nll_cpue = 0.0;

  // likelihood for random walk / process error
  for(int i = styr + 1; i <= endyr; ++i) {
    step(log_biomass_pred(i-1), log_biomass_pred(i), log_PE);
  }
  // calculate external likelihood for output
  for(int i = styr + 1; i <= endyr; ++i) {
    nll_pe += 0.5 * (log(2.0 * M_PI * exp(2.0 * log_PE)) + square(log_biomass_pred(i) - log_biomass_pred(i-1)) / exp(2.0 * log_PE));
  }

  // likelihood for biomass survey observations
  for(int i = 1; i <= n_biomass; ++i) {
    obs_biomass(log_biomass_pred(yrs_biomass(i)),i);
  }
  // calculate external likelihood for output
  for(int i = 1; i <= n_biomass; ++i) {
    nll_biomass += 0.5 * (log(2.0 * M_PI * square(log_biomass_sd(i))) + square(log_biomass_pred(yrs_biomass(i)) - log_biomass(i)) / square(log_biomass_sd(i)));
  }

  for(int i = styr; i <= endyr; ++i) {
    cpue_pred(i) = exp(log_q + log_biomass_pred(i));
  }
  log_cpue_pred = log(cpue_pred);

  // likelihood for cpue survey observations

  for(int i = 1; i <= n_cpue; ++i) {
    // 0 = inside
    if(sep_fxn_flag == 0) {
      // inside separable function - does not match TMB
      obs_cpue_inside(log_biomass_pred(yrs_cpue(i)),log_q,i);
    }

    // 1 = outside
    if(sep_fxn_flag == 1) {
      //obs_cpue_outside(exp(log_q + log_biomass_pred(yrs_cpue(i))),i);
  	   obs_cpue_outside(log_cpue_pred(yrs_cpue(i)),i);
    }
  }

    // calculate external likelihood for output (include small likelihood weight)
  for(int i = 1; i <= n_cpue; ++i) {
    nll_cpue += 0.5 * 0.5 * (log(2.0 * M_PI * square(log_cpue_sd(i))) + square(log_cpue_pred(yrs_cpue(i)) - log(cpue(i))) / square(log_cpue_sd(i)));
  }

 // jnll = nll_pe + nll_biomass + nll_cpue;

  if (sd_phase()) {
    biomass_pred = exp(log_biomass_pred);
  }

// separable function for random effects component of the likelihood - note dnorm() doesn't work
SEPARABLE_FUNCTION void step(const dvariable& biom1, const dvariable& biom2, const dvariable& tmp_log_PE)
   jnll += 0.5 * (log(2.0 * M_PI * exp(2.0 * tmp_log_PE)) + square(biom2 - biom1) / exp(2.0 * tmp_log_PE));

// separable function for the biomass data
SEPARABLE_FUNCTION void obs_biomass(const dvariable& tmp_log_biomass_pred, int i)
  jnll += 0.5 * (log(2.0 * M_PI * square(log_biomass_sd(i))) + square(tmp_log_biomass_pred - log_biomass(i)) / square(log_biomass_sd(i)));

// separable function for the cpue data when predicted log cpue is calculated INSIDE of the function. add 0.5 likelihood weight to cpue component of the likelihood 
SEPARABLE_FUNCTION void obs_cpue_inside(const dvariable& tmp_log_biomass_pred, const dvariable& tmp_log_q, int i)
  jnll += 0.5 * 0.5 * (log(2.0 * M_PI * square(log_cpue_sd(i))) + square((tmp_log_biomass_pred + tmp_log_q) - log(cpue(i))) / square(log_cpue_sd(i)));

// separable function for the cpue data when predicted log cpue is calculated OUTSIDE of the function - note that there is a likelihood weight of 0.5 on here.. doesn't want to converge without it
SEPARABLE_FUNCTION void obs_cpue_outside(const dvariable& tmp_log_cpue_pred, int i)
  jnll += 0.5 * 0.5 * (log(2.0 * M_PI * square(log_cpue_sd(i))) + square(tmp_log_cpue_pred - log_cpue(i)) / square(log_cpue_sd(i)));

REPORT_SECTION
  report << "nll_pe \t" << nll_pe << endl;
  report << "nll_biomass \t" << nll_biomass << endl;
  report << "nll_cpue \t" << nll_cpue << endl;
  report << "jnll \t" << jnll << endl;