
cal.nu2s.IIE_Mjo <- function(med, alphas_IIE, folds, learners_regressions, control, covariates_cl = FALSE, benchmark_covar = NULL) {

    q1 <- estimate_q1.IIE_Mjo(med, alphas_IIE, folds, learners_regressions, control, covariates_cl = FALSE, benchmark_covar = NULL)

    nu2s_list <- nu2s.IIE_Mjo(med, alphas_IIE, q1)
    nu2s_list
}

cal.nu2s.IDE <- function(med, alphas_IDE, folds, learners_regressions, control, covariates_cl = FALSE, benchmark_covar = NULL) {

    q1 <- estimate_q1.IDE(med, alphas_IDE, folds, learners_regressions, control, covariates_cl = FALSE, benchmark_covar = NULL)

    nu2s_list <- nu2s.IDE(med, alphas_IDE, q1)
    nu2s_list
}


nu2s.IIE_Mjo <- function(med, alphas_IIE, q1 = q1_IIE_Mjo) {
    # for IIE, ay is the same for both potential outcomes, am differs
    # alphas[[2]][["alpha2_target"]]=alphas[[1]][["alpha2_target"]]

    natural <- list()
    natural[["alpha1"]] <- alphas_IIE[["alpha1"]]
    natural[["alpha2"]] <- alphas_IIE[["alpha2"]]

    shifted <- list()
    shifted[["alpha1_target"]] <- alphas_IIE[["alpha1_target"]]
    # q1[["b2_alpha2"]] == alphas_IIE[["alpha2_target"]]
    shifted[["alpha2_target"]] <- natural[["alpha1"]] * (q1[["b2_alpha2"]] - q1[["q1_fit"]]) + q1[["b1_q1"]]

    eif_nu2 <- vector("list", length = 2)
    names(eif_nu2) <- glue("nu{1:2}")
    # nu2
    eif_nu2[["nu1"]] <- 2 * shifted$alpha1_target - natural$alpha1^2
    eif_nu2[["nu2"]] <- 2 * shifted$alpha2_target - natural$alpha2^2

    nu2 <- map(eif_nu2, \(x) mean(x))
    nu2plugin <- map(natural, \(x) (mean(x^2)))

    cluster_id <- med@vars@cluster_id
    if (is.na(cluster_id)) {
        S_cluster <- NA
    } else {
        S_cluster <- med@data[[cluster_id]]
    }
    se_nu2 <- inference.cl(
        data.frame(
            nu1 = eif_nu2[["nu1"]],
            nu2 = eif_nu2[["nu2"]]
        ), S = S_cluster, conf.level = 0.95)
    
    list(nu2 = nu2, eif_nu2 = eif_nu2, nu2plugin = nu2plugin,
         se_nu2 = c(nu1 = se_nu2$std_error[1], nu2 = se_nu2$std_error[2]))
}

nu2s.IDE <- function(med, alphas_IDE, q1 = q1_IDE) {
    # for IDE, am is the same for both potential outcomes,  ay differs
    natural <- list()
    natural[["alpha1"]] <- alphas_IDE[["alpha1"]]
    natural[["alpha2"]] <- alphas_IDE[["alpha2"]]

    shifted <- list()
    shifted[["alpha1_target"]] <- alphas_IDE[["alpha1_target"]]
    shifted[["alpha2_target"]] <- natural[["alpha1"]] * (q1[["b2_alpha2"]] - q1[["q1_fit"]]) + q1[["b1_q1"]]

    eif_nu2 <- vector("list", length = 2)
    names(eif_nu2) <- glue("nu{1:2}")
    # nu2
    eif_nu2[["nu1"]] <- 2 * shifted$alpha1_target - natural$alpha1^2
    eif_nu2[["nu2"]] <- 2 * shifted$alpha2_target - natural$alpha2^2

    nu2 <- map(eif_nu2, \(x) mean(x))
    nu2plugin <- map(natural, \(x) (mean(x^2)))
    
    cluster_id <- med@vars@cluster_id
    if (is.na(cluster_id)) {
        S_cluster <- NA
    } else {
        S_cluster <- med@data[[cluster_id]]
    }
    se_nu2 <- inference.cl(
        data.frame(
            nu1 = eif_nu2[["nu1"]],
            nu2 = eif_nu2[["nu2"]]
        ), S = S_cluster, conf.level = 0.95)
    
    list(nu2 = nu2, eif_nu2 = eif_nu2, nu2plugin = nu2plugin,
         se_nu2 = c(nu1 = se_nu2$std_error[1], nu2 = se_nu2$std_error[2]))
}

