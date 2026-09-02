
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                        		   %
%       Panel Cointegrating Polynomial Regressions:			   % 
%	Group-Mean Fully Modified OLS Estimation and Inference           % 
%						    			   %
%       by Martin Wagner and Karsten Reichold				   %
%									   %
%       DOI: https://doi.org/10.1080/07474938.2023.2178141              %
%									   % 
%	GitHub: https://github.com/kreichold/GroupMeanFM	           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Key function:
-------------------------------

-	GroupMeanFMOLS.m
	(This is the only function that needs to be executed by the practitioner. It returns group-mean fully modified OLS estimation results.)

(For more details see the descriptions in the Matlab-file.)


 
List of all functions included:
-------------------------------

-	And_HAC91.m (This function implements the data-dependent bandwidth selection rule of Andrews (1991) and is used in lr_varmod.m.)
-	demean_detrend.m (This function allows to demean and detrend time series.)
-	lr_varmod.m (This function estimates long-run covariance matrices.)
-	lr_weights.m (This function computes kernel weights used in lr_varmod.m.)

(More information is given in the corresponding Matlab-files.)


