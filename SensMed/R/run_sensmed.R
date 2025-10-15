#' Sensitivity analysis for causal mediation analysis
#'
#' Implements a sensitivity analysis method for unmeasured pretreatment confounders in assessing causal mediation effects based on the debiased machine learning approach.
#'
#' @param data A \code{data.frame} containing all analysis variables.
#' @param treatment The column name of the treatment variable. Currently supports binary treatment.
#' @param outcome The column name of the outcome variable.
#' @param mediators A character vector containing the column names of the mediator variables.
#' @param covariates A character vector containing the column names of pretreatment covariates.
#' @param learners A character vector of \code{mlr3superlearner} methods for estimating models. Default is \code{learners = c("glm")}. All available methods could be found at \link{https://mlr-org.com/learners.html}; examples include "ranger" (random forest), "nnet" (neural network), "lightgbm" (light implementation of gradient boosting), "xgboost" (extreme gradient boosting), "rpart" (regression tree), "glmnet" (GLM with elastic net regularization), "earth" (multivariate adaptive splines).
#' @param cf.y (optional) A \code{numeric} vector of the R-squared sensitivity parameters for the regressions, that is, R^2_\{Y-theta_\{2,s\} ~ theta_2 -theta_\{2,s\} \} and R^2_\{bM-theta_\{1,s\} ~ theta_\{1,u,s\} -theta_\{1,s\} \}.
#' @param cf.m (optional) A \code{numeric} vector of the R-squared sensitivity parameters for the weights, that is, 1-R^2_\{alpha_2 ~ alpha_\{2,s\}\} and 1-R^2_\{alpha_1 ~ alpha_\{1,s\}\}.
#' @param rho (optional) A \code{numeric} vector of the correlation sensitivity parameters, that is, rho_\{1\} and rho_\{2\}.
#' @param RV.options (optional) A list of the options for computing the robustness value specified following the template \code{RV.options = list(theta.null = 0, cf = seq(0, 0.99, by = 0.01), rho = c(1, 1))}, where
#' "theta.null" is the null value,
#' "cf" is the vector of R-squared sensitivity parameter values to search over, and
#' "rho" is a vector of the correlation sensitivity parameters rho_\{1\} and rho_\{2\}.
# @param benchmark_covariates (optional) A character vector containing the column names of pretreatment covariates to be used as benchmarks for the strengths of the unmeasured confounders.
# @param if_contour (optional) Whether to return results for making contour plots.
# @param contour_grid (optional) A list of two \code{numeric} vectors that contain the x- and y-coordinates of grid lines. By default, \code{contour_grid = list(x = seq(0, 1, by = 0.03), y = seq(0, 1, by = 0.03))}.
#' @param control (optional) A list of control parameters specified via \code{estimate_control()}. The default setting is based on R package \code{crumble}.
#'
#' @export
#'

run.sensmed <- function(data,
                        treatment,
                        outcome,
                        mediators,
                        covariates,
                        learners = c("glm"),
                        cf.y = c(0.26, 0.26), # large R2
                        cf.m = c(0.26, 0.26),
                        rho = c(1, 1),
                        # other settings
                        RV.options = list(theta.null = 0, cf = seq(0, 0.99, by = 0.01), rho = c(1, 1)),
                        control = estimate_control()
) {
    rho2 <- rho^2

    sens_M1 <- Sens.Mediate(
        data = data,
        treatment = treatment,
        outcome = outcome,
        mediators = mediators,
        covariates = covariates,
        cluster_id = NA_character_,
        learners_outcome = learners,
        learners_weights = learners,
        learners_else = learners,
        cf.y = cf.y,
        cf.m = cf.m,
        rho2 = rho2,
        # known long regressions and long weights
        long_mus_IIE = NULL,
        long_alphas_IIE = NULL,
        long_mus_IDE = NULL,
        long_alphas_IDE = NULL,
        # other settings
        conf.level = 0.95, # significance level 0.05 for testing short parameter; 90% CI for Confidence Bounds, with one-sided significance level 0.05
        benchmark_covariates = NULL,
        control = control,
        RV.options = RV.options
    )

    list(
        tab_sensitivity_reporting = sens_M1$tab_sensitivity_reporting,
        sens_res = sens_M1$sens_res,
        short_list = sens_M1$short_list)
}

