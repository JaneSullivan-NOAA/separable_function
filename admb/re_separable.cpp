#ifdef DEBUG
  #ifndef __SUNPRO_C
    #include <cfenv>
    #include <cstdlib>
  #endif
#endif
#ifdef DEBUG
  #include <chrono>
#endif
#include <admodel.h>
#ifdef USE_ADMB_CONTRIBS
#include <contrib.h>

#endif
#include <df1b2fun.h>

#include <adrndeff.h>

  extern "C"  {
    void ad_boundf(int i);
  }
#include <re_separable.htp>

  df1b2_parameters * df1b2_parameters::df1b2_parameters_ptr=0;
  model_parameters * model_parameters::model_parameters_ptr=0;
model_data::model_data(int argc,char * argv[]) : ad_comm(argc,argv)
{
  adstring tmpstring;
  tmpstring=adprogram_name + adstring(".dat");
  if (argc > 1)
  {
    int on=0;
    if ( (on=option_match(argc,argv,"-ind"))>-1)
    {
      if (on>argc-2 || argv[on+1][0] == '-')
      {
        cerr << "Invalid input data command line option"
                " -- ignored" << endl;
      }
      else
      {
        tmpstring = adstring(argv[on+1]);
      }
    }
  }
  global_datafile = new cifstream(tmpstring);
  if (!global_datafile)
  {
    cerr << "Error: Unable to allocate global_datafile in model_data constructor.";
    ad_exit(1);
  }
  if (!(*global_datafile))
  {
    delete global_datafile;
    global_datafile=NULL;
  }
  sep_fxn_flag.allocate("sep_fxn_flag");
  styr.allocate("styr");
  endyr.allocate("endyr");
  yrs.allocate(styr,endyr);
 yrs.fill_seqadd(styr,1);
  n_biomass.allocate("n_biomass");
  yrs_biomass.allocate(1,n_biomass,"yrs_biomass");
  biomass.allocate(1,n_biomass,"biomass");
  biomass_cv.allocate(1,n_biomass,"biomass_cv");
  log_biomass.allocate(1,n_biomass);
 log_biomass = log(biomass);
  log_biomass_sd.allocate(1,n_biomass);
 log_biomass_sd = elem_prod(biomass_cv,biomass_cv) + 1.0;
 log_biomass_sd = sqrt(log(log_biomass_sd));
  n_cpue.allocate("n_cpue");
  yrs_cpue.allocate(1,n_cpue,"yrs_cpue");
  cpue.allocate(1,n_cpue,"cpue");
  cpue_cv.allocate(1,n_cpue,"cpue_cv");
  log_cpue.allocate(1,n_cpue);
 log_cpue = log(cpue);
  log_cpue_sd.allocate(1,n_cpue);
 log_cpue_sd = elem_prod(cpue_cv,cpue_cv) + 1.0;
 log_cpue_sd = sqrt(log(log_cpue_sd));
  if (global_datafile)
  {
    delete global_datafile;
    global_datafile = NULL;
  }
}

model_parameters::model_parameters(int sz,int argc,char * argv[]) : 
 model_data(argc,argv) , function_minimizer(sz)
{
  model_parameters_ptr=this;
  initializationfunction();
  log_PE.allocate("log_PE");
  log_q.allocate("log_q");
  log_biomass_pred.allocate(styr,endyr,"log_biomass_pred");
  biomass_pred.allocate(styr,endyr,"biomass_pred");
  cpue_pred.allocate(styr,endyr,"cpue_pred");
  log_cpue_pred.allocate(styr,endyr,"log_cpue_pred");
  prior_function_value.allocate("prior_function_value");
  likelihood_function_value.allocate("likelihood_function_value");
  jnll.allocate("jnll");  /* ADOBJECTIVEFUNCTION */
}
void model_parameters::userfunction(void)
{
  jnll =0.0;
  jnll = 0.0;
  //nll_pe = 0.0;
  //nll_biomass = 0.0;
  //nll_cpue = 0.0;
  // likelihood for random walk / process error
  for(int i = styr + 1; i <= endyr; ++i) {
    step(log_biomass_pred(i-1), log_biomass_pred(i), log_PE);
  }
  // likelihood for biomass survey observations
  for(int i = 1; i <= n_biomass; ++i) {
    obs_biomass(log_biomass_pred(yrs_biomass(i)), i);
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
      obs_cpue_inside(log_biomass_pred(yrs_cpue(i)), log_q, i);
    }
    // 1 = outside
    if(sep_fxn_flag == 1) {
      // obs_cpue_outside(exp(log_q + log_biomass_pred(yrs_cpue(i))),i);
  	  obs_cpue_outside(log_cpue_pred(yrs_cpue(i)), i);
    }
  }
 //jnll = nll_pe + nll_biomass + nll_cpue;
  if (sd_phase()) {
    biomass_pred = exp(log_biomass_pred);
    //log_biomass_pred2 = log_biomass_pred;
  }
}

void SEPFUN1  model_parameters::step(const dvariable& biom1, const dvariable& biom2, const dvariable& tmp_log_PE)
{
  begin_df1b2_funnel();
  // dvariable tmp_PE = exp(tmp_log_PE);
  // jnll += dnorm(biom1, biom2, tmp_PE, 1);
  jnll += 0.5 * (log(2.0 * M_PI * mfexp(2.0 * tmp_log_PE))) + square(biom2 - biom1) / (mfexp(2.0 * tmp_log_PE));
  end_df1b2_funnel();
}

void SEPFUN1  model_parameters::obs_biomass(const dvariable& tmp_log_biomass_pred, int i)
{
  begin_df1b2_funnel();
  // jnll += dnorm(tmp_log_biomass_pred, log_biomass(i), log_biomass_sd(i), true);
  jnll += 0.5 * log(2.0 * M_PI * square(log_biomass_sd(i))) + square(tmp_log_biomass_pred - log_biomass(i)) / square(log_biomass_sd(i));
  end_df1b2_funnel();
}

void SEPFUN1  model_parameters::obs_cpue_inside(const dvariable& tmp_log_biomass_pred, const dvariable& tmp_log_q, int i)
{
  begin_df1b2_funnel();
  jnll += 0.5 * log(2.0 * M_PI * square(log_cpue_sd(i))) + square((tmp_log_biomass_pred + tmp_log_q) - log_cpue(i)) / square(log_cpue_sd(i));
  end_df1b2_funnel();
}

void SEPFUN1  model_parameters::obs_cpue_outside(const dvariable& tmp_log_cpue_pred, int i)
{
  begin_df1b2_funnel();
  //nll_cpue += dnorm(tmp_log_cpue_pred, log_cpue(i), log_cpue_sd(i), true);
  jnll += 0.5 * log(2.0 * M_PI * square(log_cpue_sd(i))) + square(tmp_log_cpue_pred - log_cpue(i)) / square(log_cpue_sd(i));
  end_df1b2_funnel();
}

void model_parameters::report(const dvector& gradients)
{
 adstring ad_tmp=initial_params::get_reportfile_name();
  ofstream report((char*)(adprogram_name + ad_tmp));
  if (!report)
  {
    cerr << "error trying to open report file"  << adprogram_name << ".rep";
    return;
  }
  //report << "log_biomass_pred" << log_biomass_pred2 << endl;
  //report << "log_cpue_pred" << log_cpue_pred << endl;
  //report << "Like components" << endl;
  //report << "pe " << nll_pe << endl;
  //report << "biomass survey  " << nll_biomass << endl;
  //report << "cpue survey " << nll_cpue << endl;
  report << "total nll " << jnll << endl;
}
  long int arrmblsize=0;

int main(int argc,char * argv[])
{
#ifdef DEBUG
  auto start = std::chrono::high_resolution_clock::now();
  #ifndef __SUNPRO_C
std::feclearexcept(FE_ALL_EXCEPT);
  #endif
#endif
  ad_set_new_handler();
  ad_exit=&ad_boundf;
    gradient_structure::set_NO_DERIVATIVES();
    gradient_structure::set_YES_SAVE_VARIABLES_VALUES();
      if (!arrmblsize) arrmblsize=150000;
    df1b2variable::noallocate=1;
df1b2variable::pool = new adpool();
initial_df1b2params::varsptr = new P_INITIAL_DF1B2PARAMS[1000];
{
    df1b2_parameters mp(arrmblsize,argc,argv);
    mp.iprint=10;

    function_minimizer::random_effects_flag=1;
    df1b2variable::noallocate=0;
    mp.preliminary_calculations();
    initial_df1b2params::separable_flag=1;
    mp.computations(argc,argv);
}
delete [] init_df1b2variable::list;
init_df1b2variable::list = NULL;
delete [] initial_df1b2params::varsptr;
initial_df1b2params::varsptr = NULL;
delete df1b2variable::pool;
df1b2variable::pool = NULL;
#ifdef DEBUG
  std::cout << endl << argv[0] << " elapsed time is " << std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - start).count() << " microseconds." << endl;
  #ifndef __SUNPRO_C
bool failedtest = false;
if (std::fetestexcept(FE_DIVBYZERO))
  { cerr << "Error: Detected division by zero." << endl; failedtest = true; }
if (std::fetestexcept(FE_INVALID))
  { cerr << "Error: Detected invalid argument." << endl; failedtest = true; }
if (std::fetestexcept(FE_OVERFLOW))
  { cerr << "Error: Detected overflow." << endl; failedtest = true; }
if (std::fetestexcept(FE_UNDERFLOW))
  { cerr << "Error: Detected underflow." << endl; }
if (failedtest) { std::abort(); } 
  #endif
#endif
    return 0;
}

extern "C"  {
  void ad_boundf(int i)
  {
    /* so we can stop here */
    exit(i);
  }
}

void model_parameters::preliminary_calculations(void){
  #if defined(USE_ADPVM)

  admaster_slave_variable_interface(*this);

  #endif

}

model_data::~model_data()
{}

model_parameters::~model_parameters()
{}

void model_parameters::final_calcs(void){}

void model_parameters::set_runtime(void){}

#ifdef _BORLANDC_
  extern unsigned _stklen=10000U;
#endif


#ifdef __ZTC__
  extern unsigned int _stack=10000U;
#endif

void df1b2_parameters::user_function(void)
{
  jnll =0.0;
  jnll = 0.0;
  //nll_pe = 0.0;
  //nll_biomass = 0.0;
  //nll_cpue = 0.0;
  // likelihood for random walk / process error
  for(int i = styr + 1; i <= endyr; ++i) {
    step(log_biomass_pred(i-1), log_biomass_pred(i), log_PE);
  }
  // likelihood for biomass survey observations
  for(int i = 1; i <= n_biomass; ++i) {
    obs_biomass(log_biomass_pred(yrs_biomass(i)), i);
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
      obs_cpue_inside(log_biomass_pred(yrs_cpue(i)), log_q, i);
    }
    // 1 = outside
    if(sep_fxn_flag == 1) {
      // obs_cpue_outside(exp(log_q + log_biomass_pred(yrs_cpue(i))),i);
  	  obs_cpue_outside(log_cpue_pred(yrs_cpue(i)), i);
    }
  }
 //jnll = nll_pe + nll_biomass + nll_cpue;
  if (sd_phase()) {
    biomass_pred = exp(log_biomass_pred);
    //log_biomass_pred2 = log_biomass_pred;
  }
}

void   df1b2_pre_parameters::step(const funnel_init_df1b2variable& biom1, const funnel_init_df1b2variable& biom2, const funnel_init_df1b2variable& tmp_log_PE)
{
  begin_df1b2_funnel();
  // df1b2variable tmp_PE = exp(tmp_log_PE);
  // jnll += dnorm(biom1, biom2, tmp_PE, 1);
  jnll += 0.5 * (log(2.0 * M_PI * mfexp(2.0 * tmp_log_PE))) + square(biom2 - biom1) / (mfexp(2.0 * tmp_log_PE));
  end_df1b2_funnel();
}

void   df1b2_pre_parameters::obs_biomass(const funnel_init_df1b2variable& tmp_log_biomass_pred, int i)
{
  begin_df1b2_funnel();
  // jnll += dnorm(tmp_log_biomass_pred, log_biomass(i), log_biomass_sd(i), true);
  jnll += 0.5 * log(2.0 * M_PI * square(log_biomass_sd(i))) + square(tmp_log_biomass_pred - log_biomass(i)) / square(log_biomass_sd(i));
  end_df1b2_funnel();
}

void   df1b2_pre_parameters::obs_cpue_inside(const funnel_init_df1b2variable& tmp_log_biomass_pred, const funnel_init_df1b2variable& tmp_log_q, int i)
{
  begin_df1b2_funnel();
  jnll += 0.5 * log(2.0 * M_PI * square(log_cpue_sd(i))) + square((tmp_log_biomass_pred + tmp_log_q) - log_cpue(i)) / square(log_cpue_sd(i));
  end_df1b2_funnel();
}

void   df1b2_pre_parameters::obs_cpue_outside(const funnel_init_df1b2variable& tmp_log_cpue_pred, int i)
{
  begin_df1b2_funnel();
  //nll_cpue += dnorm(tmp_log_cpue_pred, log_cpue(i), log_cpue_sd(i), true);
  jnll += 0.5 * log(2.0 * M_PI * square(log_cpue_sd(i))) + square(tmp_log_cpue_pred - log_cpue(i)) / square(log_cpue_sd(i));
  end_df1b2_funnel();
}
   
void df1b2_pre_parameters::setup_quadprior_calcs(void) 
{ 
  df1b2_gradlist::set_no_derivatives(); 
  quadratic_prior::in_qp_calculations=1; 
}  
  
void df1b2_pre_parameters::begin_df1b2_funnel(void) 
{ 
  (*re_objective_function_value::pobjfun)=0; 
  other_separable_stuff_begin(); 
  f1b2gradlist->reset();  
  if (!quadratic_prior::in_qp_calculations) 
  { 
    df1b2_gradlist::set_yes_derivatives();  
  } 
  funnel_init_var::allocate_all();  
}  
 
void df1b2_pre_parameters::end_df1b2_funnel(void) 
{  
  lapprox->do_separable_stuff(); 
  other_separable_stuff_end(); 
  funnel_init_var::deallocate_all(); 
} 
  
void model_parameters::begin_df1b2_funnel(void) 
{ 
  if (lapprox)  
  {  
    {  
      begin_funnel_stuff();  
    }  
  }  
}  
 
void model_parameters::end_df1b2_funnel(void) 
{  
  if (lapprox)  
  {  
    end_df1b2_funnel_stuff();  
  }  
} 

void df1b2_parameters::deallocate() 
{
  log_PE.deallocate();
  log_q.deallocate();
  log_biomass_pred.deallocate();
  biomass_pred.deallocate();
  cpue_pred.deallocate();
  log_cpue_pred.deallocate();
  prior_function_value.deallocate();
  likelihood_function_value.deallocate();
  jnll.deallocate();
} 
void df1b2_parameters::allocate(void) 
{
  log_PE.allocate("log_PE");
  log_q.allocate("log_q");
  log_biomass_pred.allocate(styr,endyr,"log_biomass_pred");
  biomass_pred.allocate(styr,endyr,"biomass_pred");
  cpue_pred.allocate(styr,endyr,"cpue_pred");
  log_cpue_pred.allocate(styr,endyr,"log_cpue_pred");
  prior_function_value.allocate("prior_function_value");
  likelihood_function_value.allocate("likelihood_function_value");
  jnll.allocate("jnll");  /* ADOBJECTIVEFUNCTION */
}
