# ------------------------------------------------------------------------------
# CBF_LMM.R
# Paper-facing entry point.
#
# The manuscript ("Conditional Bayes-Factor Stepwise Refinement for Fine-Mapping
# in Gaussian Linear Mixed Models") calls the primary procedure CBF-LMM:
# the marginal conditional Bayes factor with numerical marginalization of the
# variance-component ratio delta (MBF) and eBIC stopping.
#
# In this code base that procedure is MS_L_LMM_stepwise_fast() with its default
# arguments (delta_eval = "marginal", criterion = "eBIC").  CBF_LMM_stepwise()
# below is a thin alias so that code and paper share one name.
#
# Sensitivity options (as in the paper):
#   delta_eval = "reml"            plug-in REML evaluation of delta
#   criterion  = "JointPosterior"  posterior-score stopping
#
# JS_L_LMM_stepwise_fast() is the exploratory joint Schur-complement variant
# discussed as future work in the manuscript; it is retained for exploratory
# use and is not part of the primary method specification.
# ------------------------------------------------------------------------------

CBF_LMM_stepwise <- function(y, X, tau2 = 0.04, K_max = 10L,
                             criterion = "eBIC",
                             delta_eval = "marginal", ...) {
  MS_L_LMM_stepwise_fast(y, X, tau2 = tau2, K_max = K_max,
                         criterion = criterion,
                         delta_eval = delta_eval, ...)
}
