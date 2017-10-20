# basal-area-application

This folder contains data and code for modeling zero-inflated Black Cherry basal areas via Yang and Tokdar's 2016 method of quantile regression.

* Data is in main folder.

* Code is the "code" folder:
   - Files starting with "6Eco" include basal-are, application-specific code for reading in data, running models,
   assessing convergence, creating coefficient plots, performing diagnostics, and running k-fold comparisons
   to Tobit Regression and independently-estimated quantile regressions.
   - Files starting with "0" include general functions, e.g. performing training and test calculations for
   k-fold split, summarizing posteriors, finding zero-probability estimates, performing convergence diagnostics, etc.
