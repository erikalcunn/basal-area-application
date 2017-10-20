# basal-area-application

Contains data and code for modeling zero-inflated Black Cherry basal areas via Yang and Tokdar's 2016 method of quantile regression.

* Data is in main repo directory.

* Code is the "code" directory:
   - "6Eco*.R" files include basal-area, application-specific code for reading in data, running models,
   assessing convergence, creating coefficient plots, performing diagnostics, and running k-fold comparisons
   to Tobit Regression and independently-estimated quantile regressions.
   - "0*.R" files include general or helper functions, e.g. performing training and test calculations for
   k-fold split, summarizing posteriors, finding zero-probability estimates, performing convergence diagnostics, etc.
