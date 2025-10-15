\donttest{
    library(glue)
    library(purrr)
    library(tidyverse)
    library(origami)
    library(mlr3superlearner)
    library(mlr3extralearners)


    data("simulated_data")
    dat <- aa

    # RV and sensitivity results at specified sensitivity parameter values

    sens_M1 <- SensMed::run.sensmed(
        data = dat,
        treatment = "A",
        outcome = "Y",
        mediators = "M",
        covariates = c("X1obs", "X2obs"),
        learners = c("mean", "glm", "ranger"),
        cf.y = c(0.26, 0.26),
        cf.m = c(0.26, 0.26),
        rho = c(1, 1),
        RV.options = list(theta.null = 0, cf = seq(0, 0.99, by = 0.01), rho = c(1, 1))
    )

    # Robustness values, and estimates, standard errors, and confidence intervals for the short parameter
    sens_M1$tab_sensitivity_reporting
    # Sensitivity results at the specified sensitivity parameter values: estimates, standard errors, and confidence intervals for the short parameter (theta_s), lower bound (theta_m), and upper bound (theta_p)
    sens_M1$sens_res

}
