estimate.short <- function(
        data,
        treatment,
        outcome,
        mediators,
        covariates,
        cluster_id = NA_character_,
        learners_outcome = c("glm"),
        learners_weights = c("glm"),
        learners_else = c("glm"),
        # known long regressions and long weights for evaluating the bias formula
        long_mus_IIE = NULL,
        long_alphas_IIE = NULL,
        long_mus_IDE = NULL,
        long_alphas_IDE = NULL,
        # other settings
        conf.level = 0.95,
        # run the long short regression theta_{1,u,s}
        benchmark_covariates = "U",
        control = estimate_control()
) {
    data$id <- 1:nrow(data)
    med <- med_data(
        data = data,
        vars = med_vars(
            A = treatment,
            Y = outcome,
            M = mediators,
            C = covariates,
            cluster_id = cluster_id,
            id = "id"
        ),
        weights = rep(1, nrow(data)),
        a0 = 0,
        a1 = 1
    )
    # if (!is.na(cluster_id)) {
    #     # dummy indicators
    #     Sdumm <- dummy_cols(data[[cluster_id]], remove_first_dummy = FALSE, remove_selected_columns = TRUE)
    #     colnames(Sdumm) <- paste0("S", 1:ncol(Sdumm))
    #     Sname_dummies <- colnames(Sdumm)
    #
    #     # cluster means
    #     data <- data %>%
    #         mutate(
    #             id = 1:n(),
    #             S = data[[cluster_id]]) %>%
    #         group_by(S) %>%
    #         mutate(across(
    #             all_of(c(treatment, mediators, outcome, covariates)),
    #             list(clmean = ~ mean(.), cwc = ~ . - mean(.))
    #         )) %>%
    #         bind_cols(Sdumm) %>%
    #         ungroup()
    #
    #     med <- med_data(
    #         data = data,
    #         vars = med_vars(
    #             A = treatment,
    #             Y = outcome,
    #             M = mediators,
    #             C = covariates,
    #             # C = c(glue("{covariates}_cwc"), Sname_dummies),
    #             cluster_id = cluster_id,
    #             id = "id"
    #         ),
    #         weights = rep(1, nrow(data)),
    #         a0 = 0,
    #         a1 = 1
    #     )
    # }

    m_each_Mjo <- list(
        y01 = c(am = "data_0", ay = "data_1"),
        y11 = c(am = "data_1", ay = "data_1"),
        y00 = c(am = "data_0", ay = "data_0"),
        y10 = c(am = "data_1", ay = "data_0")
    )
    # m_each_Mjo[c(1,2)]
    m_each_IIE <- list(
        IIE_Mjo_m = c(am = "data_0", ay = "data_1"),
        IIE_Mjo_p = c(am = "data_1", ay = "data_1")
    )
    # m_each_Mjo[c(3,1)]
    m_each_IDE <- list(
        IDE_m = c(am = "data_0", ay = "data_0"),
        IDE_p = c(am = "data_0", ay = "data_1")
    )

    if (is.na(cluster_id)) {
        folds <- origami::make_folds(med@data, V = control$crossfit_folds)
    }
    if (!is.na(cluster_id)) {
        folds <- origami::make_folds(med@data, cluster_ids = med@data[[cluster_id]], V = control$crossfit_folds)
        # folds <- make.fold_K(med@data, Snames = "S", control$crossfit_folds)
    }
    if (control$crossfit_folds == 1) {
        folds[[1]]$training_set <- folds[[1]]$validation_set
    }

    # estimate models -----------------------

    # if (!is.na(cluster_id)) {
    #     med_cwc <- med
    #     med_cwc@vars@A <- glue("{med@vars@A}_cwc")
    #     med_cwc@vars@M <- glue("{med@vars@M}_cwc")
    #     med_cwc@vars@Y <- glue("{med@vars@Y}_cwc")
    #     # med_cwc@vars@C <- glue("{med@vars@C}_cwc")
    #     med_cwc@vars@C <- c(glue("{med@vars@C}_cwc"), Sname_dummies)
    # }

    # for each potential outcome
    mus_Mjo <- estimate_mu_cl(med, m_each_Mjo, folds, learners_outcome, control, covariates_cl = FALSE)

    a_c <- a.c(med, folds, learners_weights, control)
    a_mc <- a.mc(med, folds, learners_weights, control)
    if (length(med@vars@M) > 1) {
        m_ac <- NULL
    } else {
        if(length(unique(med@data[, med@vars@M])) == 2) {
            m_ac <- m.ac(med, folds, learners_weights, control)
        } else {
            m_ac <- NULL
        }
    }
    
    # m_ac <- NULL
    alphas_Mjo <- fit_alpha_each(med, a_c, a_mc, m_ac, m_each_Mjo, folds, learners_weights, control)
    # alphas_Mjo <- estimate_alpha_each(med, m_each_Mjo, folds, nn_module, control)

    # for IIE
    mus_IIE <- with(mus_Mjo, {data.frame(
        mu_fit_am_jo = y11$mu_fit_am_jo,
        mu_fit_ay_jo = y11$mu_fit_ay_jo,
        b_am_jo = y11$b_am_jo - y01$b_am_jo,
        b_ay_jo = y11$b_ay_jo)})

    alphas_IIE <- with(alphas_Mjo, {data.frame(
        alpha1 = y11$alpha1 - y01$alpha1,
        alpha2 = y11$alpha2 - y01$alpha2,
        # long.short_alpha2 = y11$long.short_alpha2 - y01$long.short_alpha2,
        # alpha1_IIE(1,X)-alpha1_IIE(0,X) = (alpha1(1, X)-0) - (0-alpha1(0, X)) = alpha1(1, X) + alpha1(0, X)
        alpha1_target = y11$alpha1_target + y01$alpha1_target,
        # alpha2_IIE(1,M,X) = alpha2(1, M, X, am=1, ay=1) - alpha2(1, M, X, am=0, ay=1)
        alpha2_target = y11$alpha2_target - y01$alpha2_target)})


    mus_IDE <- with(mus_Mjo, {data.frame(
        mu_fit_am_jo = y01$mu_fit_am_jo - y00$mu_fit_am_jo,
        mu_fit_ay_jo = y01$mu_fit_ay_jo,
        b_am_jo = y01$b_am_jo - y00$b_am_jo,
        b_ay_jo = y01$b_ay_jo - y00$b_ay_jo)})

    alphas_IDE <- with(alphas_Mjo, {data.frame(
        alpha1 = y01$alpha1,
        alpha2 = y01$alpha2 - y00$alpha2,
        alpha1_target = y01$alpha1_target,
        # (alpha2_IDE(1,M,X,am=0,ay=1)-alpha2_IDE(1,M,X,am=0,ay=0)) - (alpha2_IDE(0,M,X,am=0,ay=1)-alpha2_IDE(0,M,X,am=0,ay=0)) = alpha2_IDE(1,M,X,am=0,ay=1) + alpha2_IDE(0,M,X,am=0,ay=0)
        alpha2_target = y01$alpha2_target + y00$alpha2_target)})

    # mean((alphas_Mjo$y01$alpha2 - gen_data$y01$alpha_each$alpha2)^2)
    # mean(alphas_Mjo$y01$alpha2 * data$Y)
    # mean(gen_data$y01$alpha_each$alpha2 * data$Y)
    # mean(gen_data$mus_IIE$b_am_jo)
    # mean(mus_IIE$b_am_jo)

    ## short EIF -----
    eifs_each <- eif_each_cal(med, mus_Mjo, alphas_Mjo, m_each_Mjo)
    eif_IDE <- eif_each_cal(med, mus_IDE, alphas_IDE, m_each = NULL)
    eif_IIE <- eif_each_cal(med, mus_IIE, alphas_IIE, m_each = NULL)
    eifs <- data.frame(Reduce(cbind, eifs_each), eif_IIE, eif_IDE)
    names(eifs) <- c(names(eifs_each), "eif_IIE", "eif_IDE")

    if (is.null(cluster_id) | is.na(cluster_id)) {
        estimates <- inference.iid(eifs)
    }
    if (!is.null(cluster_id) & !is.na(cluster_id)) {
        estimates <- inference.cl(eifs, med@data[[cluster_id]], average = "individual", conf.level = 0.95)
    }
    #
    # estimates <- estimates %>%
    #     rownames_to_column(var = "estimand") %>%
        # arrange(estimand)

    # Sensitivity analysis ------------
    ## nu2s, sigma2s -------
    cal_nu2s_IIE <- cal.nu2s.IIE_Mjo(med, alphas_IIE, folds, learners_else, control, covariates_cl = FALSE, benchmark_covar = NULL)
    cal_nu2s_IDE <- cal.nu2s.IDE(med, alphas_IDE, folds, learners_else, control, covariates_cl = FALSE, benchmark_covar = NULL)

    cal_sigma2s_IIE <- cal.sigma2s.IIE_Mjo(med@data[[med@vars@Y]], mus_IIE)
    # debias hat.sigma1sq for the estimation of theta2
    # cal_sigma1sq_IIE <- cal.sigma1sq.IIE(med, a_mc, mus_IIE, folds, learners_outcome, control, covariates_cl = FALSE, benchmark_covar = NULL)
    cal_sigma1sq_IIE <- cal.sigma1sq.IIE(med, a_mc, mus_IIE, mus_Mjo)
    cal_sigma2s_IIE$eif_sigma2s$sigma1 <- cal_sigma1sq_IIE$eif_sigma1sq
    if (cal_sigma1sq_IIE$sigma1sq > 0) {
        cal_sigma2s_IIE$sigma2s$sigma1 <- cal_sigma1sq_IIE$sigma1sq
    }

    cal_sigma2s_IDE <- cal.sigma2s.IDE(med@data[[med@vars@Y]], mus_IDE)
    # cal_sigma1sq_IDE <- cal.sigma1sq.IDE(med, a_mc, mus_IDE, folds, learners_outcome, control, covariates_cl = FALSE, benchmark_covar = NULL)
    cal_sigma1sq_IDE <- cal.sigma1sq.IDE(med, a_mc, mus_IDE, mus_Mjo)
    cal_sigma2s_IDE$eif_sigma2s$sigma1 <- cal_sigma1sq_IDE$eif_sigma1sq
    if (cal_sigma1sq_IDE$sigma1sq > 0) {
        cal_sigma2s_IDE$sigma2s$sigma1 <- cal_sigma1sq_IDE$sigma1sq
    }

    # evaluate
    # unlist(cal_nu2s_IIE$nu2) - gen_largeN$nu2s_IIE
    # unlist(cal_nu2s_IDE$nu2) - gen_largeN$nu2s_IDE
    # unlist(cal_sigma2s_IIE$sigma2s) - gen_largeN$sigma2s_IIE
    # unlist(cal_sigma2s_IDE$sigma2s) - gen_largeN$sigma2s_IDE

    ## long.short regression --------

    if (!is.null(benchmark_covariates)){
        # get theta1_u,s(A,X)
        med_U <- med # benchmark_covariates = "U"
        med_U@vars@C <- base::union(med@vars@C, benchmark_covariates) # the long covariates
        long.short_mus_Mjo <- estimate_mu.long.short(benchmark_covariates, med_U, m_each_Mjo, folds, learners_outcome, control)

        long.short_mus_IIE <- with(long.short_mus_Mjo, {data.frame(
            mu_fit_am_jo = y11$mu_fit_am_jo,
            mu_fit_ay_jo = y11$mu_fit_ay_jo,
            b_am_jo = y11$b_am_jo - y01$b_am_jo,
            b_ay_jo = y11$b_ay_jo)})
        long.short_mus_IDE <- with(long.short_mus_Mjo, {data.frame(
            mu_fit_am_jo = y01$mu_fit_am_jo - y00$mu_fit_am_jo,
            mu_fit_ay_jo = y01$mu_fit_ay_jo,
            b_am_jo = y01$b_am_jo - y01$b_am_jo,
            b_ay_jo = y01$b_ay_jo - y00$b_ay_jo)})

        # bias formula empirical -----
        if (!is.null(long_mus_IIE)) {
            # bias formula (covariance of regression error and weighting error)
            # IIE
            est_bias_IIE <- -mean((long_mus_IIE$mu_fit_ay_jo - mus_IIE$mu_fit_ay_jo) * (long_alphas_IIE$alpha2 - alphas_IIE$alpha2))  -mean(
                (long.short_mus_IIE$mu_fit_am_jo - mus_IIE$mu_fit_am_jo) * (long_alphas_IIE$alpha1 - alphas_IIE$alpha1))

            # # empirical bias
            # # gen_data$IIE is expected to eqaul mean(gen_data$eif_IIEeta);
            # emp_bias_IIE <- mean(eif_IIE) - mean(gen_data$eif_IIEeta)
            # est_bias_IIE - emp_bias_IIE

            # IDE
            est_bias_IDE <- -mean((long_mus_IDE$mu_fit_ay_jo - mus_IDE$mu_fit_ay_jo) * (long_alphas_IDE$alpha2 - alphas_IDE$alpha2)) -mean(
                (long.short_mus_IDE$mu_fit_am_jo - mus_IDE$mu_fit_am_jo) * (long_alphas_IDE$alpha1 - alphas_IDE$alpha1))
            # # empirical bias
            # emp_bias_IDE <- mean(eif_IDE) - mean(gen_data$IDE)
            # est_bias_IDE - emp_bias_IDE

        }
    }



    mget(ls(), envir = environment())
}
