# Simulated-QRJM

This code is in reference to the paper:
Kundu, D. and Das, K., 2024. "A quantile-regression approach to bivariate longitudinal joint modeling." Journal of Statistical Research 58(1) , pp : 111-130 , doi: 10.3329/jsr.v58i1.75417

link : https://www.researchgate.net/publication/383265217_A_quantile-regression_approach_to_bivariate_longitudinal_joint_modeling

Simulation code in : Simulation QRJM.R
Plots : 

(1) Data contour plot 

(2) Convergence plot (to show parameter convergence) 

(3) Longitudinal plot for the estimated qunatiles

(4) Estimated survival plots

(5) Parameter significance plots (based on 95% credible set).

#Plots and their meaning

Quantile_Vs_tau : shows no quantile crossing

Multi-normal_QQ_plot : serves as a motivation for joint modeling, citing deviation from bivariate normality

Contour : Shows the bivariate -density of mean responses across subjects.

Medicine : shows the estimate and 95% CI for the "betas" corresponding to the medicine in the longitudinal submodel across quantiles

Fixed_Var : Tile plot for showing significance of coefficients for each "fixed variable" across (submodel X quantile)

Association_Convergence : Plot showing convergence of the association parameter estimates (note: on an average Gelman-Rubin convergence stat < 1.1)

Association: shows the estimate and 95% CI for the "psis" corresponding to the medicine in the survival submodel across quantiles

Longitudinal : Estimated bi-variate longitudinal quantiles across time for various quantiles

Survival : Estimated non relapse probabilities for various quantiles






