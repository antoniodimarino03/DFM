# DFM
Matlab and R functions to estimate the dynamic factor models and to perform forecasting as in Stock and Watson 2002, Forni Hallin Lippi Reichlin 2000, Forni Hallin Lippi Reichlin 2005, Forni Hallin Lippi Zaffaroni 2017.
1 gdfm_twosided.m estimates the Generalized Dynamic Factor Model as in Forni, Hallin, Lippi, and Reichlin (2000);
2 gdfm_onesided.m estimates and forecasts the Generalized Dynamic Factor Model as in Forni, Hallin, Lippi, and Reichlin (2005); 
3 gdfm_unrestricted.m estimates and forecasts the Generalized Dynamic Factor Model as in Forni, Hallin, Lippi, and Zaffaroni (2017);
4
5 numfactors.m estimates the number of factors as in Hallin and Liška (2007);
6 spectral.m computes the spectral density decomposition used by all other functions;
7 prepare_missing.m  transforms the raw data into stationary form
8 mrsq.m  computes the R-squared and marginal R-squared  values from estimated factors and factor loadings
9
10
11

The codes for 1 2 3 6 were written by Matteo Barigozzi, Mario Forni, Roman Liška, and Matteo
Luciani, were commented and debugged by Matteo Barigozzi and are available from
www.barigozzi.eu/codes.html.

The codes fore 4 5 7 8 ...
