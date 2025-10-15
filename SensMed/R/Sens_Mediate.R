
Sens.Mediate <- function(
    data,
    treatment,
    outcome,
    mediators,
    covariates,
    cluster_id = NA_character_,
    learners_outcome = c("glm"),
    learners_weights = c("glm"),
    learners_else = c("glm"),
    cf.y = 0.13,
    cf.m = 0.13,
    rho2 = 1,
    # known long regressions and long weights
    long_mus_IIE = NULL,
    long_alphas_IIE = NULL,
    long_mus_IDE = NULL,
    long_alphas_IDE = NULL,
    # other settings
    conf.level = 0.95, # significance level 0.05 for testing the short parameter theta_s
    benchmark_covariates = NULL,
    control = estimate_control(),
    RV.options = NULL # RV.options = list(theta.null = 0, cf = seq(0, 0.99, by = 0.05), rho = 1)
) {

    # get short components, empirical bias
    short_list <- estimate.short(data,
                                 treatment,
                                 outcome,
                                 mediators,
                                 covariates,
                                 cluster_id,
                                 learners_outcome,
                                 learners_weights,
                                 learners_else,
                                 # known long regressions and long weights
                                 long_mus_IIE,
                                 long_alphas_IIE,
                                 long_mus_IDE,
                                 long_alphas_IDE,
                                 # other settings
                                 conf.level,
                                 benchmark_covariates,
                                 control)


    list2env(short_list, envir = environment())


    ## confidence bounds at given sens params ----
    # significance.level <- 1-conf.level

    sens_IIE <- sens_Mjo(med, eif_IIE, cal_nu2s_IIE, cal_sigma2s_IIE, cf.y, cf.m, rho2, only_YM_confounded = FALSE, estimate = estimates["eif_IIE", ])

    sens_IDE <- sens_Mjo(med, eif_IDE, cal_nu2s_IDE, cal_sigma2s_IDE, cf.y, cf.m, rho2, only_YM_confounded = FALSE, estimate = estimates["eif_IDE", ])

    # RV.options = list(theta.null = 0, cf = seq(0, 0.99, by = 0.05), rho2 = 1)
    if (!is.null(RV.options)) {
        RV.options$rho2 <- RV.options$rho^2

        rv_IIE <- rv.value(med, eif_IIE, cal_nu2s_IIE, cal_sigma2s_IIE, RV.options$theta.null, RV.options$cf, RV.options$rho2, only_YM_confounded = FALSE, estimates["eif_IIE", ])

        rv_IDE <- rv.value(med, eif_IDE, cal_nu2s_IDE, cal_sigma2s_IDE, RV.options$theta.null, RV.options$cf, RV.options$rho2, only_YM_confounded = FALSE, estimates["eif_IDE", ])

        ## Minimal sensitivity reporting ----
        rv_values <- map_dbl(list(IIE=rv_IIE, IDE=rv_IDE), \(x) x$RV_null0_sig0.05)

        tab_sensitivity_reporting <- estimates %>%
            rownames_to_column(var = "estimand") %>%
            mutate(estimand = gsub("eif_", "", estimand)) %>%
            dplyr::filter(str_detect(estimand, "I")) %>%
            # mutate(CI = glue("[{round(`95% CI_low`, 2)}, {round(`95% CI_high`, 2)}]")) %>%
            mutate(RV = rv_values[estimand])

    }

    sens_res <- rbind(IIE = sens_IIE$sens_estimates,
                      IDE = sens_IDE$sens_estimates) %>%
        rownames_to_column(var = "estimand")



    out <- mget(ls(envir = environment()))
    out
}
