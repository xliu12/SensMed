

sens_Mjo <- function(med, eif_short, cal_nu2s, cal_sigma2s, cf.y = c(0.04, 0.04), cf.m = c(0.03, 0.03), rho2 = 1, only_YM_confounded = FALSE, estimate) {
    # alphas <- alphas_IIE; mus <- mus_IIE
    # alphas <- alphas_IDE; mus <- mus_IDE
    Y <- med@data[[med@vars@Y]]

    # if (params == "IIE_Mjo") {
    #     cal_nu2s <- cal.nu2s.IIE_Mjo(med, alphas, folds, learners_regressions, control, covariates_cl = FALSE, benchmark_covar = NULL)
    #     cal_sigma2s <- cal.sigma2s.IIE_Mjo(Y, mus)
    # }
    # if (params == "IDE") {
    #     cal_nu2s <- cal.nu2s.IDE(med, alphas, folds, learners_regressions, control, covariates_cl = FALSE, benchmark_covar = NULL)
    #     cal_sigma2s <- cal.sigma2s.IDE(Y, mus)
    # }

    nu2 <- cal_nu2s$nu2
    nu2[nu2<0] <- unlist(cal_nu2s$nu2plugin)[nu2<0]
    eif_nu2 <- cal_nu2s$eif_nu2

    psi_nu2 <- map(1:2, \(y) eif_nu2[[y]] - nu2[[y]])

    resY <- cal_sigma2s$resY
    R2_Y <- cal_sigma2s$R2_Y
    eif_sigma2 <- cal_sigma2s$eif_sigma2
    sigma2 <- cal_sigma2s$sigma2s

    psi_sigma2 <- map(1:2, \(y) eif_sigma2[[y]] - sigma2[[y]]) # eif_sigma2 - sigma2

    # S2
    S2 <- map(1:2, \(x) sigma2[[x]] * nu2[[x]])

    # eif
    # eif_short <- eif_each(med, mus_1, alpha_each1, m_each1)
    # eif_short <- eif_each_cal(med, mus, alphas, m_each)
    # names(eif_short) <- names(m_each)

    param_short <- mean(eif_short)
    psi_param_short <- eif_short - param_short

    # bias bound

    if (length(cf.y) < 2) {
        cf.y = c(cf.y, cf.y)
    }
    if (length(cf.m) < 2) {
        cf.m = c(cf.m, cf.m)
    }
    if (length(rho2) < 2) {
        rho2 = c(rho2, rho2)
    }
    bf <- lapply(1:2, \(k) bias.factor(cf.y = cf.y[k], cf.m = cf.m[k], rho2 = rho2[k]))
    # bf <- bias.factor(cf.y = cf.y, cf.m = cf.m, rho2 = rho2)
    # if (length(bf) < 2) {
    #     if (only_YM_confounded) {
    #         bf <- c(0, bf)
    #     } else {
    #         bf <- c(bf, bf)
    #     }
    # }
    if (only_YM_confounded) {
        bf[[1]] <- 0
    }
    bias <- map(1:2, \(k) bf[[k]] * sqrt(S2[[k]]))

    psi_bias <- map(1:2, \(k) (bf[[k]] / 2) * (1 / sqrt(S2[[k]])) * (sigma2[[k]] * psi_nu2[[k]] + nu2[[k]] * psi_sigma2[[k]]))


    bias_sum <- reduce(bias, `+`)

    # psi_S2[[1]][1:3, ] + psi_S2[[2]][1:3, ]
    psi_bias_sum <- reduce(psi_bias, `+`)

    # currently, focusing on individual average effect
    param_bound_p <- param_short + bias_sum
    param_bound_m <- param_short - bias_sum


    psi_param_bound_p <- psi_param_short + psi_bias_sum
    psi_param_bound_m <- psi_param_short - psi_bias_sum

    cluster_id <- med@vars@cluster_id

    if (is.na(cluster_id)) {
        S_cluster <- NA
    } else {
        S_cluster <- med@data[[cluster_id]]
    }

    se_bounds <- inference.cl(
        data.frame(
            # psi_param_bound_m, psi_param_bound_p
            theta_m = psi_param_bound_m + param_bound_m,
            theta_p = psi_param_bound_p + param_bound_p
        ),
        S = S_cluster,
        conf.level = 0.9
    )

    se_theta_s <- estimate
    if (is.na(cluster_id)) {
        se_theta_m <- se_bounds["theta_m", c("estimate", "std_error", "90% CI_low", "90% CI_high")]
        se_theta_p <- se_bounds["theta_p", c("estimate", "std_error", "90% CI_low", "90% CI_high")]
    }
    if (!is.na(cluster_id)) {
        se_theta_m <- se_bounds["theta_m", ]
        se_theta_p <- se_bounds["theta_p", ]
    }

    sens_estimates <- data.frame(
        theta_s = se_theta_s,
        theta_m = se_theta_m,
        theta_p = se_theta_p
    ) %>%
        rename_with(., .fn = ~ gsub("..CI", "%CI", .)) %>%
        mutate(cf.y1 = cf.y[1], cf.y2 = cf.y[2],
               cf.m1 = cf.m[1], cf.m2 = cf.m[2],
               rho21 = rho2[1], rho22 = rho2[2])

    mget(ls(), envir = environment())
}


# Bias factor sensitivity parameter ------------
# based on R package dml.sensemakr and Chernozhukov et al. (2024) Long story short
bias.factor <- function(cf.y = 0.03, cf.m = 0.04, rho2 = 1) {
    sqrt(rho2 * cf.y * (cf.m / (1 - cf.m)))
}



# eif_short=eif_IIE; cal_nu2s=cal_nu2s_IIE; cal_sigma2s=cal_sigma2s_IIE; estimate=estimates["eif_IIE", ]

rv.value <- function(med, eif_short, cal_nu2s, cal_sigma2s, theta.null = 0, cf = seq(0, 0.99, by = 0.05), rho2 = 1, only_YM_confounded = FALSE, estimate) {


    values <- map2(
        .x = cf, .y = cf,
        .f = ~ sens_Mjo(med, eif_short, cal_nu2s, cal_sigma2s, cf.y = .x, cf.m = .y, rho2, only_YM_confounded, estimate)[["sens_estimates"]]
    )

    df_values <- do.call(rbind, values)
    # view(df_values)

    df_values1 <- df_values %>%
        rename(conf_bound.l = `theta_m.90%CI_low`, conf_bound.u = `theta_p.90%CI_high`) %>%
        mutate(
            # Rsq = cf,
            cf = cf, theta.null = theta.null,
            if_cover = (conf_bound.l - theta.null) * (conf_bound.u - theta.null) < 0
        )

    # original nonsignificant, reverse to significant
    if (df_values1[df_values1$cf.m1 == 0 & df_values1$cf.y1 == 0, "if_cover"] == TRUE) {
        if (df_values1$theta_s.estimate[1] < 0) {
            rv <- df_values1$cf[which(df_values1$`theta_m.90%CI_high` <= 0)][1]
        }
        if (df_values1$theta_s.estimate[1] >= 0) {
            rv <- df_values1$cf[which(df_values1$`theta_p.90%CI_low` >= 0)][1]
        }
    }
    # original significant, reverse to nonsignificant
    if (df_values1[df_values1$cf.m1 == 0 & df_values1$cf.y1 == 0, "if_cover"] == FALSE) {
        rv <- df_values1$cf[df_values1$if_cover == TRUE][1]
    }
    rv[is.na(rv)] <- max(df_values1)


    rv_res <- list(rv, df_values1)
    names(rv_res) <- c(glue("RV_null{theta.null}_sig{0.05}"), glue("RVsens_res"))

    rv_res
}

