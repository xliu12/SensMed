
cal.nu2s.IIE_Mjo <- function(med, alphas_IIE, folds, learners_regressions, control, covariates_cl = FALSE, benchmark_covar = NULL) {

    q1 <- estimate_q1.IIE_Mjo(med, alphas_IIE, folds, learners_regressions, control, covariates_cl = FALSE, benchmark_covar = NULL)

    nu2s_list <- nu2s.IIE_Mjo(med, alphas_IIE, q1)
    nu2s_list
}

cal.nu2s.IDE <- function(med, alphas_IDE, folds, learners_regressions, control, covariates_cl = FALSE, benchmark_covar = NULL) {

    # m_each = list(IDE_m = c(az = "", am = "data_0", ay = "data_0"),
    #        IDE_p = c(az = "", am = "data_0", ay = "data_1"))

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
    # debiased pseudo outcome / doubly robust transformation
    # q1[["b2_alpha2"]] == alphas_IIE[["alpha2_target"]]
    shifted[["alpha2_target"]] <- natural[["alpha1"]] * (q1[["b2_alpha2"]] - q1[["q1_fit"]]) + q1[["b1_q1"]]

    # verify with data generation
    # mean(2*gen_data$alphas_IIE$alpha1_target - gen_data$alphas_IIE$alpha1^2)
    # mean(gen_data$alphas_IIE$alpha1^2)
    # mean(2*alphas_IIE$alpha1_target - alphas_IIE$alpha1^2)

    # mean(gen_data$alphas_IIE$alpha2^2)
    # mean(2*gen_data$alphas_IIE$alpha1*gen_data$alphas_IIE$alpha2_target - gen_data$alphas_IIE$alpha2^2)
    # mean(2*alphas_IIE$alpha1*alphas_IIE$alpha2_target - alphas_IIE$alpha2^2)

    eif_nu2 <- vector("list", length = 2)
    names(eif_nu2) <- glue("nu{1:2}")
    # nu2
    eif_nu2[["nu1"]] <- 2 * shifted$alpha1_target - natural$alpha1^2
    eif_nu2[["nu2"]] <- 2 * shifted$alpha2_target - natural$alpha2^2

    nu2 <- map(eif_nu2, \(x) mean(x))
    nu2plugin <- map(natural, \(x) (mean(x^2)))
    # nu2[["nu2"]] <- ifelse(nu2[["nu2"]] < 0, mean(natural[["alpha2"]]^2), nu2[["nu2"]])
    # nu2[["nu1"]] <- ifelse(nu2[["nu1"]] < 0, mean(natural[["alpha1"]]^2), nu2[["nu1"]])

    # check nu2 estimator SE
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
    # nu2[[1]] <- se_nu2$estimate[1]
    # nu2[[2]] <- se_nu2$estimate[2] # can differ if cluster-average


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
    # debiased pseudo outcome / doubly robust transformation
    shifted[["alpha2_target"]] <- natural[["alpha1"]] * (q1[["b2_alpha2"]] - q1[["q1_fit"]]) + q1[["b1_q1"]]

    eif_nu2 <- vector("list", length = 2)
    names(eif_nu2) <- glue("nu{1:2}")
    # nu2
    eif_nu2[["nu1"]] <- 2 * shifted$alpha1_target - natural$alpha1^2
    eif_nu2[["nu2"]] <- 2 * shifted$alpha2_target - natural$alpha2^2

    nu2 <- map(eif_nu2, \(x) mean(x))
    nu2plugin <- map(natural, \(x) (mean(x^2)))
    # nu2[["nu2"]] <- ifelse(nu2[["nu2"]] < 0, mean(natural[["alpha2"]]^2), nu2[["nu2"]])
    # nu2[["nu1"]] <- ifelse(nu2[["nu1"]] < 0, mean(natural[["alpha1"]]^2), nu2[["nu1"]])

    # check nu2 estimator SE
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
    # nu2[[1]] <- se_nu2$estimate[1]
    # nu2[[2]] <- se_nu2$estimate[2] # can differ if cluster-average

    list(nu2 = nu2, eif_nu2 = eif_nu2, nu2plugin = nu2plugin,
         se_nu2 = c(nu1 = se_nu2$std_error[1], nu2 = se_nu2$std_error[2]))
}


# cal.nu2.IIE_M1 <- function(alphas) {
#     # alphas=alphas_each_IIE_M1_r
#     natural <- list()
#     if (length(alphas) == 2) {
#         natural <- map(1:3, \(k) alphas[[2]][[glue("alpha{k}")]] - alphas[[1]][[glue("alpha{k}")]] )
#
#         shifted <- list()
#         # shifted[["alpha1_target"]] <- alphas[[2]][["alpha1_target"]] - alphas[[1]][["alpha1_target"]]
#         # shifted[["alpha2_target"]] <- alphas[[2]][["alpha2_target"]] * alphas[[2]][["alpha1"]]  - alphas[[1]][["alpha2_target"]] * alphas[[1]][["alpha1"]]
#         # shifted[["alpha3_target"]] <- alphas[[2]][["alpha3_target"]] * alphas[[2]][["alpha2"]]  - alphas[[1]][["alpha3_target"]] * alphas[[1]][["alpha2"]]
#
#         shifted[["alpha1_target"]] <- alphas[[2]][["alpha1_target"]] + alphas[[1]][["alpha1_target"]]
#         shifted[["alpha2_target"]] <- alphas[[2]][["alpha2_target"]] * alphas[[2]][["alpha1"]]  + alphas[[1]][["alpha2_target"]] * alphas[[1]][["alpha1"]]
#         shifted[["alpha3_target"]] <- alphas[[2]][["alpha3_target"]] * alphas[[2]][["alpha2"]]  + alphas[[1]][["alpha3_target"]] * alphas[[1]][["alpha2"]]
#
#     }
#     # alphas=alphas_each_IIE_M1_r[1]
#     if (length(alphas) == 1) {
#         natural <- map(1:3, \(k) alphas[[1]][[glue("alpha{k}")]] )
#         # shifted <- map(1:3, \(k) alphas[[1]][[glue("alpha{k}_target")]] )
#         shifted <- list()
#         shifted[["alpha1_target"]] <- alphas[[1]][["alpha1_target"]]
#         shifted[["alpha2_target"]] <- alphas[[1]][["alpha2_target"]] * alphas[[1]][["alpha1"]]
#         shifted[["alpha3_target"]] <- alphas[[1]][["alpha3_target"]] * alphas[[1]][["alpha2"]]
#     }
#
#     names(natural) <- glue("alpha{1:3}")
#     names(shifted) <- glue("alpha{1:3}_target")
#
#
#     eif_nu2 <- vector("list", length = 3)
#     names(eif_nu2) <- glue("nu{1:3}") # make sure all things are ordered as 1,2,3
#     # nu2
#     # eif_nu2[["nu1"]] <- 2 * shifted$alpha1_target * 1 - natural$alpha1^2
#     # eif_nu2[["nu2"]] <- 2 * shifted$alpha2_target * natural$alpha1 - natural$alpha2^2
#     # eif_nu2[["nu3"]] <- 2 * shifted$alpha3_target * natural$alpha2 - natural$alpha3^2
#
#     eif_nu2[["nu1"]] <- 2 * shifted$alpha1_target - natural$alpha1^2
#     eif_nu2[["nu2"]] <- 2 * shifted$alpha2_target - natural$alpha2^2
#     eif_nu2[["nu3"]] <- 2 * shifted$alpha3_target - natural$alpha3^2
#
#     ## currently, focusing on the individual-average effects
#     # individual average for the nuisance parameters
#
#     nu2 <- map(eif_nu2, \(x) {
#         mean(x)
#     })
#     # nu2 <- vector("list", length = 3)
#     names(nu2) <- glue("nu{1:3}")
#     # nu2[["nu1"]] <- mean(natural$alpha1^2)
#     # nu2[["nu2"]] <- mean(natural$alpha2^2)
#     # nu2[["nu3"]] <- mean(natural$alpha3^2)
#
#
#     list(nu2 = nu2, eif_nu2 = eif_nu2)
# }

# alphas <- alphas_each_IIE_Mjo[1:2]

