# DFM
Matlab and R functions to estimate the dynamic factor models and to perform forecasting as in Stock and Watson (2002); Forni, Hallin, Lippi and Reichlin (2000); Forni, Hallin, Lippi and Reichlin (2005); Forni, Hallin, Lippi and Zaffaroni (2017).

1) gdfm_twosided.m estimates the Generalized Dynamic Factor Model as in Forni, Hallin, Lippi, and Reichlin (2000);
2) gdfm_onesided.m estimates and forecasts the Generalized Dynamic Factor Model as in Forni, Hallin, Lippi, and Reichlin (2005); 
3) gdfm_unrestricted.m estimates the Generalized Dynamic Factor Model as in Forni, Hallin, Lippi, and Zaffaroni (2017);
4) factors_SW.m estimates a set of factors for a given dataset as in Stock and Watson (2002);
5) numfactors.m estimates the number of dynamic factors as in Hallin and Liška (2007);
6) spectral.m computes the spectral density decomposition used by all other functions;
7) prepare_missing.m  transforms the raw data into stationary form;
8) mrsq.m  computes the R-squared and marginal R-squared values from estimated factors and factor loadings;
9) mrsq_1F.m computes R-squared and marginal R-squared for a SINGLE dynamic factor model estimated as in Forni, Hallin, Lippi and Reichlin (2000);
10) gdfm17_forecast.m performs forecasting as in Forni, Hallin, Lippi and Zaffaroni (2017);
11) forecast_SW.m  performs forecasting as in Stock and Watson (2002);
12) Cleaning_1.R performs the first cleaning operations of the dataset and split the data into training set and test set, it also estimates the number of static factors as in Bai and Ng (2003)
13) Cleaning_2.m performs advanced cleaning operations of the dataset ending up with raw data;
14) MAPE and PLOTS.R computes MAPE for the performed forecasts for each method, variable and horizon. It also displays the plots of the absolute percentage errors for each variable. 

Codes 1, 2, 3, 6, 13 were written by Matteo Barigozzi, Mario Forni, Roman Liška, and Matteo Luciani, were commented and debugged by Matteo Barigozzi and are available from www.barigozzi.eu/codes.html.
Codes 4, 7, 8 have been taken by the codes associated with FRED-QD at the folllowing link: https://www.stlouisfed.org/research/economists/mccracken/fred-databases.
Codes 2, 3, 7 have been slightly modified by the undersigned. 
Codes 9, 10, 11, 12, 14 have been written by the undersigned. These scripts incorporates elements generated or refined using AI tools.
Code 9 has been adpated to a frequency domain approach by the undersigned from code 8.
