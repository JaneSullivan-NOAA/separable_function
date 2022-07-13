# separable_function
The goal of this repo is to better understand the behavior of the `SEPARABLE_FUNCTION` in ADMB. We get different parameter estimates and model results dependening on if the calculations for a predicted value are specified [inside](https://github.com/JaneSullivan-NOAA/separable_function/blob/f917b8c6a12796d3605392b92e939212e404bb26/admb/re_separable.tpl#L121) or [outside](https://github.com/JaneSullivan-NOAA/separable_function/blob/f917b8c6a12796d3605392b92e939212e404bb26/admb/re_separable.tpl#L97) the `SEPARABLE_FUNCTION`. We demonstrate these differences using a random walk model used for survey averaging that accepts a biomass index and an additional CPUE index. Biomass (`log_biomass_pred`) is estimated as a vector of random effects, and the process error (`log_PE`) and scaling parameter (`log_q`) are estimated as fixed effects. We find that the ADMB version with the predicted values calculated inside the `SEPARABLE_FUNCTION` yields the same results as TMB. We assume the TMB/ADMB inside version of the model is "correct" but do not understand why. The full example with ADMB inside, ADMB outside, and TMB versions can be reproduced using the [separable_function_admb.R](https://github.com/JaneSullivan-NOAA/separable_function/blob/main/separable_function_admb.R) script file.

Other oddities:
In the two data examples provided, the ADMB outside version of the model is sensitive and only converges when a likelihood weight is applied. For consistency, this likelihood weight (0.5) is added to the ADMB inside and TMB versions as well. Additionally we found that we were unable to use `dnorm()` inside the `SEPARABLE_FUNCTION`. 

# session info
```
─ Session info ────────────────────────────────────────────────────────────────────
 setting  value
 version  R version 4.2.0 (2022-04-22 ucrt)
 os       Windows 10 x64 (build 19044)
 system   x86_64, mingw32
 ui       RStudio
 language (EN)
 collate  English_United States.utf8
 ctype    English_United States.utf8
 tz       America/Anchorage
 date     2022-07-13
 rstudio  2022.02.3+492 Prairie Trillium (desktop)
 pandoc   NA

─ Packages ────────────────────────────────────────────────────────────────────────
 ! package     * version date (UTC) lib source
   assertthat    0.2.1   2019-03-21 [1] CRAN (R 4.2.0)
   brio          1.1.3   2021-11-30 [1] CRAN (R 4.2.0)
   cachem        1.0.6   2021-08-19 [1] CRAN (R 4.2.0)
   callr         3.7.0   2021-04-20 [1] CRAN (R 4.2.0)
   cli           3.3.0   2022-04-25 [1] CRAN (R 4.2.0)
   colorspace    2.0-3   2022-02-21 [1] CRAN (R 4.2.0)
   crayon        1.5.1   2022-03-26 [1] CRAN (R 4.2.0)
   DBI           1.1.3   2022-06-18 [1] CRAN (R 4.2.0)
   desc          1.4.1   2022-03-06 [1] CRAN (R 4.2.0)
   devtools    * 2.4.3   2021-11-30 [1] CRAN (R 4.2.0)
   digest        0.6.29  2021-12-01 [1] CRAN (R 4.2.0)
   dplyr       * 1.0.9   2022-04-28 [1] CRAN (R 4.2.0)
   ellipsis      0.3.2   2021-04-29 [1] CRAN (R 4.2.0)
   fansi         1.0.3   2022-03-24 [1] CRAN (R 4.2.0)
   farver        2.1.0   2021-02-28 [1] CRAN (R 4.2.0)
   fastmap       1.1.0   2021-01-25 [1] CRAN (R 4.2.0)
   fs            1.5.2   2021-12-08 [1] CRAN (R 4.2.0)
   generics      0.1.2   2022-01-31 [1] CRAN (R 4.2.0)
   ggplot2     * 3.3.6   2022-05-03 [1] CRAN (R 4.2.0)
   glue          1.6.2   2022-02-24 [1] CRAN (R 4.2.0)
   gtable        0.3.0   2019-03-25 [1] CRAN (R 4.2.0)
   labeling      0.4.2   2020-10-20 [1] CRAN (R 4.2.0)
   lattice       0.20-45 2021-09-22 [1] CRAN (R 4.2.0)
   lifecycle     1.0.1   2021-09-24 [1] CRAN (R 4.2.0)
   magrittr      2.0.3   2022-03-30 [1] CRAN (R 4.2.0)
   Matrix        1.4-1   2022-03-23 [1] CRAN (R 4.2.0)
   memoise       2.0.1   2021-11-26 [1] CRAN (R 4.2.0)
   munsell       0.5.0   2018-06-12 [1] CRAN (R 4.2.0)
   pillar        1.7.0   2022-02-01 [1] CRAN (R 4.2.0)
   pkgbuild      1.3.1   2021-12-20 [1] CRAN (R 4.2.0)
   pkgconfig     2.0.3   2019-09-22 [1] CRAN (R 4.2.0)
   pkgload       1.2.4   2021-11-30 [1] CRAN (R 4.2.0)
   prettyunits   1.1.1   2020-01-24 [1] CRAN (R 4.2.0)
   processx      3.5.3   2022-03-25 [1] CRAN (R 4.2.0)
   ps            1.7.0   2022-04-23 [1] CRAN (R 4.2.0)
   purrr         0.3.4   2020-04-17 [1] CRAN (R 4.2.0)
   R6            2.5.1   2021-08-19 [1] CRAN (R 4.2.0)
   Rcpp          1.0.8.3 2022-03-17 [1] CRAN (R 4.2.0)
   remotes       2.4.2   2021-11-30 [1] CRAN (R 4.2.0)
   rlang         1.0.2   2022-03-04 [1] CRAN (R 4.2.0)
   rprojroot     2.0.3   2022-04-02 [1] CRAN (R 4.2.0)
   rstudioapi    0.13    2020-11-12 [1] CRAN (R 4.2.0)
   scales        1.2.0   2022-04-13 [1] CRAN (R 4.2.0)
   sessioninfo   1.2.2   2021-12-06 [1] CRAN (R 4.2.0)
   testthat    * 3.1.4   2022-04-26 [1] CRAN (R 4.2.0)
   tibble        3.1.7   2022-05-03 [1] CRAN (R 4.2.0)
   tidyr       * 1.2.0   2022-02-01 [1] CRAN (R 4.2.0)
   tidyselect    1.1.2   2022-02-21 [1] CRAN (R 4.2.0)
   TMB         * 1.9.0   2022-05-27 [1] CRAN (R 4.2.0)
   usethis     * 2.1.6   2022-05-25 [1] CRAN (R 4.2.0)
   utf8          1.2.2   2021-07-24 [1] CRAN (R 4.2.0)
   vctrs         0.4.1   2022-04-13 [1] CRAN (R 4.2.0)
   withr         2.5.0   2022-03-03 [1] CRAN (R 4.2.0)
```
