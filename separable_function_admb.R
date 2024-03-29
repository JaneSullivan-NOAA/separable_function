# compare TMB model results to two ADMB versions (one with statement outside of
# the separable function, the other inside)

# set up ----
library(dplyr)
library(tidyr)
library(ggplot2)
library(TMB)
# devtools::install_github('Cole-Monnahan-NOAA/adnuts', ref='dev')
library(adnuts)

# function to write ADMB dat files
write_dat <- function(df = NULL,
                      sep_fxn_flag,
                      fn) {
  cat("# Separable function flag (0 = inside, 1 = outside)","\n",
      sep_fxn_flag,"\n","\n",
      "# Model start and end years","\n",
      min(df$year), max(df$year), "\n","\n",
      "# Number of biomass surveys", "\n",
      length(which(!is.na(df$biomass))),"\n","\n",
      "# Biomass survey years", "\n",
      df %>% filter(!is.na(biomass)) %>% pull(year),"\n","\n",
      "# Biomass observations", "\n",
      df %>% filter(!is.na(biomass)) %>% pull(biomass),"\n","\n",
      "# Biomass CV", "\n",
      df %>% filter(!is.na(biomass)) %>% pull(biomass_cv),"\n","\n",
      "# Number of CPUE surveys", "\n",
      length(which(!is.na(df$cpue))),"\n","\n",
      "# CPUE survey years", "\n",
      df %>% filter(!is.na(cpue)) %>% pull(year),"\n","\n",
      "# CPUE observations", "\n",
      df %>% filter(!is.na(cpue)) %>% pull(cpue),"\n","\n",
      "# CPUE CV", "\n",
      df %>% filter(!is.na(cpue)) %>% pull(cpue_cv),
      "\n", sep = " ", file = paste0(fn, ".dat"))
}

# Data ----

# ye example
# df <- read.csv('ye_cseo.csv')

# sst cgoa
df <- read.csv('sst_cgoa.csv')

# TMB ----

tmbdf <- df %>%
  tidyr::expand(year = c(min(df$year): max(df$year))) %>%
  left_join(df) %>%
  arrange(year)
tmbdat <- list(yrs = tmbdf$year,
               biomass = tmbdf$biomass,
               biomass_cv = tmbdf$biomass_cv,
               cpue = tmbdf$cpue,
               cpue_cv = tmbdf$cpue_cv)
tmbpar <- list(log_PE = 1,
               log_q = 1,
               log_biomass_pred = rep(log(mean(df$biomass, na.rm = TRUE)),
                                      length(tmbdat$yrs)))

# dyn.unload(dynlib('tmb/re'))
# TMB::compile('tmb/re.cpp')
dyn.load(dynlib('tmb/re'))
obj <- MakeADFun(tmbdat, tmbpar, random = 'log_biomass_pred', dll = 're')
opt <- with(obj, nlminb(par, fn, gr))
rep <- sdreport(obj)

pars <- as.data.frame(do.call('rbind', as.list(rep, 'Est')[c(1, 2)]))
names(pars) <- c('log_est')
sd <- as.data.frame(do.call('rbind', as.list(rep, 'Std')[c(1, 2)]))
names(sd) <- c('sd')
tmb_pars <- pars %>%
  bind_cols(sd)
tmb_pars <- tmb_pars %>%
  mutate(version = 'TMB',
         variable = c('PE', 'q'),
         pred = exp(log_est),
         pred_lci = exp(log_est - 1.96 * sd),
         pred_uci = exp(log_est + 1.96 * sd)) %>%
  select(version, variable, pred, pred_lci, pred_uci)

tmb_nll <- data.frame(version = 'TMB',
                      nll = c('nll_pe', 'nll_biomass',
                      'nll_cpue', 'jnll', 'mnll'),
                      value = c(obj$report()$jnll[1],
                                 obj$report()$jnll[2],
                                obj$report()$jnll[3],
                                sum(obj$report()$jnll),
                                opt$objective))

tmb_results <- tidyr::expand_grid(variable = unique(names(rep$value)),
                                  year = tmbdat$yrs) %>%
  mutate(version = 'TMB',
         estimate = rep$value,
         sd = rep$sd,
         variable = ifelse(grepl('biomass', variable), 'biomass_pred', 'cpue_pred')) %>%
  mutate(pred = exp(estimate),
         pred_lci = exp(estimate - 1.96 * sd),
         pred_uci = exp(estimate + 1.96 * sd)) %>%
  dplyr::select(version, variable, year, pred, pred_lci, pred_uci)

# ADMB inside ----

getwd()
unlink(file.path('admb', 'inside'), recursive = T)
dir.create(file.path('admb', 'inside'))
file.copy(from = file.path('admb', 're_separable.exe'),
          to = file.path('admb/inside', 're_separable.exe'),
          overwrite = TRUE)
write_dat(df = df,
          sep_fxn_flag = 0,
          fn = 'admb/inside/inside')
setwd('admb/inside/')
write.table(as.numeric(obj$env$last.par.best), file='init.pin', row.names=FALSE, col.names=FALSE)
system('re_separable -ind inside.dat -ainp init.pin -maxfn 0 -nohess')
nll_inside_init<- read.table('RE_SEP~1.rep')
names(nll_inside_init) <- c('nll', 'value')
system('re_separable -ind inside.dat -ainp init.pin')
results_inside <- read.table('RE_SEP~1.std',header = TRUE) %>%
  as.data.frame()
nll_inside <- read.table('RE_SEP~1.rep')
names(nll_inside) <- c('nll', 'value')
setwd('../..')

# ADMB outside ----
unlink(file.path('admb', 'outside'), recursive = T)
dir.create(file.path('admb', 'outside'))
file.copy(from = file.path('admb', 're_separable.exe'),
          to = file.path('admb/outside', 're_separable.exe'),
          overwrite = TRUE)
write_dat(df = df,
          sep_fxn_flag = 1,
          fn = file.path('admb/outside/','outside'))
setwd('admb/outside/')
write.table(as.numeric(obj$env$last.par.best), file='init.pin', row.names=FALSE, col.names=FALSE)
system('re_separable -ind outside.dat -ainp init.pin -maxfn 0 -nohess')
nll_outside_init<- read.table('RE_SEP~1.rep')
names(nll_outside_init) <- c('nll', 'value')
system('re_separable -ind outside.dat -ainp init.pin')
results_outside <- read.table('RE_SEP~1.std',header = TRUE) %>%
  as.data.frame()
nll_outside <- read.table('RE_SEP~1.rep')
names(nll_outside) <- c('nll', 'value')
setwd('../..')

admb_results <- results_inside %>%
  mutate(version = 'ADMB inside') %>%
  bind_rows(results_outside %>%
              mutate(version = 'ADMB outside (status quo)')) %>%
  filter(name %in% c('log_biomass_pred', 'log_cpue_pred', 'log_q', 'log_PE')) %>%
  mutate(variable = dplyr::case_when(name == 'log_biomass_pred' ~ 'biomass_pred',
                                     name == 'log_cpue_pred' ~ 'cpue_pred',
                                     name == 'log_q' ~ 'q',
                                     name == 'log_PE' ~ 'PE'),
         pred = exp(value),
         pred_lci = exp(value - 1.96 * std.dev),
         pred_uci = exp(value + 1.96 * std.dev)) %>%
  dplyr::select(version, variable, pred, pred_lci, pred_uci)

admb_pars <- admb_results %>% filter(variable %in% c('PE', 'q'))
admb_results <- admb_results %>%
  filter(!variable %in% c('PE', 'q')) %>%
  mutate(year = rep(unique(tmb_results$year), 4))

admb_results <- results_inside %>%
  mutate(version = 'ADMB inside') %>%
  bind_rows(results_outside %>%
              mutate(version = 'ADMB outside (status quo)')) %>%
  filter(name %in% c('log_biomass_pred', 'log_cpue_pred')) %>%
  mutate(variable = ifelse(grepl('biomass', name), 'biomass_pred', 'cpue_pred'),
         pred = exp(value),
         pred_lci = exp(value - 1.96 * std.dev),
         pred_uci = exp(value + 1.96 * std.dev),
         year = rep(unique(tmb_results$year), 4)) %>%
  dplyr::select(version, variable, year, pred, pred_lci, pred_uci)

# compare results ----

# biomass survey data and model fits
obs <- df %>%
  dplyr::rename(obs = biomass, obs_cv = biomass_cv) %>%
  mutate(log_obs = ifelse(obs > 0, log(obs), NA),
                sd_log_obs = ifelse(obs > 0, sqrt(log(obs_cv^2 + 1)), NA),
                obs_lci = exp(log_obs - 1.96 * sd_log_obs),
                obs_uci = exp(log_obs + 1.96 * sd_log_obs))

ggplot(bind_rows(tmb_results, admb_results) %>%
         filter(variable == 'biomass_pred')) +
  geom_ribbon(aes(x = year, col = version, fill = version, ymin = pred_lci, ymax = pred_uci),
              alpha = 0.2, col = NA) +
  geom_line(aes(x = year, y = pred, lty = version, col = version, fill = version), size = 1) +
  geom_point(data = obs, aes(x = year, y = obs)) +
  geom_errorbar(data = obs, aes(x = year, ymin = obs_lci, ymax = obs_uci))  +
  labs(title = 'REMA model fits to biomass survey data', y = 'Biomass (t)', x = NULL) +
  scale_y_continuous(labels = scales::comma) +
  ggplot2::scale_fill_viridis_d(direction = 1) + ggplot2::scale_colour_viridis_d(direction = 1)

ggsave('biomass_survey_fits.png')

# cpue survey data and model fits
obs <- df %>%
  dplyr::rename(obs = cpue, obs_cv = cpue_cv) %>%
  mutate(log_obs = ifelse(obs > 0, log(obs), NA),
         sd_log_obs = ifelse(obs > 0, sqrt(log(obs_cv^2 + 1)), NA),
         obs_lci = exp(log_obs - 1.96 * sd_log_obs),
         obs_uci = exp(log_obs + 1.96 * sd_log_obs))

ggplot(bind_rows(tmb_results, admb_results) %>%
         filter(variable == 'cpue_pred')) +
  geom_ribbon(aes(x = year, col = version, fill = version, ymin = pred_lci, ymax = pred_uci), alpha = 0.2, col = NA) +
  geom_line(aes(x = year, y = pred, lty = version, col = version, fill = version), size = 1) +
  geom_point(data = obs, aes(x = year, y = obs)) +
  geom_errorbar(data = obs, aes(x = year, ymin = obs_lci, ymax = obs_uci))  +
  labs(title = 'REMA model fits to CPUE survey data', y = 'CPUE', x = NULL) +
  scale_y_continuous(labels = scales::comma) +
  ggplot2::scale_fill_viridis_d(direction = 1) + ggplot2::scale_colour_viridis_d(direction = 1)

ggsave('cpue_survey_fits.png')

## get the marginal NLL too
minside <- adnuts:::.read_mle_fit('re_separable', path='admb/inside')$nll
moutside <- adnuts:::.read_mle_fit('re_separable', path='admb/outside')$nll
minside <- adnuts:::.read_mle_fit('RE_SEP~1', path='admb/inside')$nll
moutside <- adnuts:::.read_mle_fit('RE_SEP~1', path='admb/outside')$nll

nlls <- cbind(tmb_nll[,-1], inside=c(nll_inside_init[,2],minside), outside=c(nll_outside_init[,2],moutside))
write.csv(nlls, 'nll_components.csv')
## # nlls
## nlls <- bind_rows(tmb_nll,
##           bind_rows(nll_inside %>%
##                       mutate(version = 'ADMB inside'),
##                     nll_outside %>%
##                       mutate(version = 'ADMB outside (status quo)'))) %>%
##   mutate(nll = ifelse(nll == 'jnll', 'total_nll', nll))

## nlls %>%
##   ggplot(aes(x = nll, y = value, col = version, fill = version)) +
##   geom_bar(stat = 'identity', position = 'dodge')

## ggsave('likelihood_components.png')

# parameter estimates
pars <- bind_rows(admb_pars, tmb_pars)
pars %>%
  ggplot(aes(x = variable, y = pred, col = version, fill = version)) +
  geom_bar(stat = 'identity', position = 'dodge') +
  geom_errorbar(aes(x = variable, ymin = pred_lci, ymax = pred_uci, group = version),
                position = 'dodge', col = 'black') +
  facet_wrap(~variable, scales = 'free_y') +
  labs(x = 'parameter', y = 'estimate')
ggsave('parameter_estimates.png')

## compare to MCMC ----
unlink(file.path('admb', 'inside'), recursive = T)
dir.create(file.path('admb', 'inside'))
file.copy(from = file.path('admb', 're_sep_mcmc.exe'),
          to = file.path('admb/inside', 're_sep_mcmc.exe'),
          overwrite = TRUE)
write_dat(df = df,
          sep_fxn_flag = 0,
          fn = 'admb/inside/inside')
init.fn <- function() {rnorm(39)}
fit_inside <- sample_nuts('re_sep_mcmc', path='admb/inside',
                          cores=1, chains=3,
                          control=list(metric=NULL),
                          skip_optimization=TRUE,
                          init=init.fn,
                          admb_args='-ind inside.dat')
## Need to run mceval manually from command line
setwd('admb/inside')
system("re_sep_mcmc -ind inside.dat -mceval")
setwd("../..")
post.inside <- read.table('admb/inside/post.dat', sep=' ',
                          header=FALSE)[,-1] %>%
  setNames(1984:2020) %>% mutate(iter=1:n()) %>%
  pivot_longer(-iter) %>%
  transmute(iter=iter, year=as.numeric(name), biomass=exp(value))


## repeat with outside
unlink(file.path('admb', 'outside'), recursive = T)
dir.create(file.path('admb', 'outside'))
file.copy(from = file.path('admb', 're_sep_mcmc.exe'),
          to = file.path('admb/outside', 're_sep_mcmc.exe'),
          overwrite = TRUE)
write_dat(df = df,
          sep_fxn_flag = 0,
          fn = 'admb/outside/outside')
init.fn <- function() {rnorm(39)}
fit_outside <- sample_nuts('re_sep_mcmc', path='admb/outside',
                          cores=1, chains=3,
                          control=list(metric=NULL),
                          skip_optimization=TRUE,
                          init=init.fn,
                          admb_args='-ind outside.dat')

## Need to run mceval manually from command line
setwd('admb/outside')
system("re_sep_mcmc -ind outside.dat -mceval")
setwd("../..")
post.outside <- read.table('admb/outside/post.dat', sep=' ',
                          header=FALSE)[,-1] %>%
  setNames(1984:2020) %>% mutate(iter=1:n()) %>%
  pivot_longer(-iter) %>%
  transmute(iter=iter, year=as.numeric(name), biomass=exp(value))




posts <- bind_rows(cbind(version='ADMB inside', post.inside),
                   cbind(version='ADMB outside (status quo)',
                         post.outside)) %>%
  group_by(year,version) %>%
  summarize(pred=median(biomass),
            pred_lci=quantile(biomass, 0.025),
            pred_uci=quantile(biomass, .975)) %>%
    mutate(type='MCMC')
results <- bind_rows(posts, cbind(filter(admb_results, variable=='biomass_pred'), type='MLE'))


ggplot(results, aes(year, pred, ymin=pred_lci, ymax=pred_uci,
  fill=type,  color=type)) +
  geom_line() + geom_ribbon(alpha=.5) + facet_wrap('version', ncol=1)
ggsave("mcmc_comparison.png")

ggplot(results,# %>%
         # bind_rows(tmb_results %>%
         #             filter(variable == 'biomass_pred') %>%
         #             select(-variable) %>%
         #             mutate(type = 'MLE')),
       aes(year, pred, ymin=pred_lci, ymax=pred_uci,
  fill=version,  color=version)) +
  geom_ribbon(alpha=.2, col = NA) +
  geom_line(size = 1) +
  facet_wrap(~type, ncol=1) +
  labs(title = 'REMA model predictions', y = 'Biomass (t)', x = NULL) +
  scale_y_continuous(labels = scales::comma) +
  ggplot2::scale_fill_viridis_d(direction = 1) + ggplot2::scale_colour_viridis_d(direction = 1)

ggsave("mcmc_comparison2.png")


