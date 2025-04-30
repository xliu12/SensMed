eif_each <- function(med, mus_1, alpha_each1, m_each1) {
    Y <- med@data[[med@vars@Y]]
    if (m_each1[["az"]] != "") {
        eif_each <- alpha_each1$alpha3 * (Y - mus_1$mu_fit3) +
            alpha_each1$alpha2 * (mus_1$b3 - mus_1$mu_fit2) +
            alpha_each1$alpha1 * (mus_1$b2 - mus_1$mu_fit1) +
            mus_1$b1
    }
    if (m_each1[["az"]] == "") {
        eif_each <- alpha_each1$alpha3 * (Y - mus_1$mu_fit_ay_jo) +
            alpha_each1$alpha2 * (mus_1$b_ay_jo - mus_1$mu_fit_am_jo) +
            mus_1$b_am_jo
    }

    eif_each1 <- data.frame(eif_each)
    names(eif_each1) <- paste0("zmy", paste(substr(m_each1, 6, 6), collapse = ""))
    eif_each1
}


eif_each_cal <- function(med, mus, alphas_each, m_each) {
    eifs <- lapply(1:length(m_each), \(x = 1) eif_each(med, mus[[x]], alphas_each[[x]], m_each[[x]]))

    eif_df <- data.frame(do.call(cbind, eifs))
    eif_df
}


# based on R package dml.sensemakr and Chernozhukov et al. (2024) Long story short

sens_IIE_M1 <- function(med, mus, alphas_each, m_each, params, cf.y = 0.04, cf.m = 0.03, rho2 = 1, only_YM_confounded = TRUE, mediatedY = c(1)) {
    # alphas_each <- alphas_each_IIE_M2_r
    # cal_nu2 <- cal.nu2.IIE_M1(alphas_each)
    cal_nu2 <- cal.nu2.IIE_M1(alphas_each[mediatedY]) # alpha's bias component is bounded by the plus-signed mean potential outcome
    nu2 <- cal_nu2$nu2
    eif_nu2 <- cal_nu2$eif_nu2

    psi_nu2 <- map(1:3, \(y) eif_nu2[[y]] - nu2[[y]])

    # sigma2
    # med <- med2_cwc
    # mus <- mus_IIE_M2
    # m_each <- m_each_IIE_M2

    Y <- med@data[[med@vars@Y]]
    resY <- cal.resY(Y, m_each, mus[mediatedY])
    # resY <- resY[, 1] # x[, 1] = x[, 2]
    eif_sigma2 <- resY^2
    # eif_sigma2 <- (Y - mus[[1]]$mu_fit3)^2 # mus[[1]]$mu_fit3 = mus[[2]]$mu_fit3 when m_each[[1]]$ay = m_each[[2]]$ay
    sigma2 <- colMeans(eif_sigma2)

    psi_sigma2 <- map(1:3, \(y) eif_sigma2[[y]] - sigma2[[y]]) # eif_sigma2 - sigma2

    # S2
    S2 <- map(1:3, \(x) sigma2[[x]] * nu2[[x]])



    # eif
    # eif_short <- eif_each(med, mus_1, alpha_each1, m_each1)
    eif_short <- eif_each_cal(med, mus, alphas_each, m_each)
    names(eif_short) <- names(m_each)

    param <- params
    if (length(m_each) == 2) {
        eif_short_each <- eif_short[, paste0(param, c("_m", "_p"))]
        param_short <- mean(eif_short_each[[2]]) - mean(eif_short_each[[1]])
        psi_param_short <- eif_short_each[[2]] - eif_short_each[[1]] - param_short
    }
    if (length(m_each) == 1) {
        eif_short_each <- eif_short
        param_short <- mean(eif_short_each[[1]])
        psi_param_short <- eif_short_each[[1]] - param_short
    }

    # bias bound

    bf <- bias.factor(cf.y = cf.y, cf.m = cf.m, rho2 = rho2)
    if (length(bf) < 3) {
        if (only_YM_confounded) {
            bf <- c(0, 0, bf)
        } else {
            bf <- c(bf, bf, bf)
        }
    }

    bias <- map(1:3, \(k) bf[[k]] * sqrt(S2[[k]]))

    psi_bias <- map(1:3, \(k) (bf[[k]] / 2) * (1 / sqrt(S2[[k]])) * (sigma2[[k]] * psi_nu2[[k]] + nu2[[k]] * psi_sigma2[[k]]))


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
        S_cluster <- med@data[["S"]]
    }
    se_bounds <- inference.cl(
        data.frame(
            theta_m = psi_param_bound_m + param_bound_m,
            theta_p = psi_param_bound_p + param_bound_p
        ),
        S = S_cluster,
        conf.level = 0.9
    )

    se_theta_s <- inference.cl(
        data.frame(
            theta_s = psi_param_short + param_short
        ),
        S = S_cluster,
        conf.level = 0.95
    )


    out <- data.frame(
        theta_s = se_theta_s,
        theta_m = se_bounds["theta_m", c("estimate", "std_error", "90% CI_low")],
        theta_p = se_bounds["theta_p", c("estimate", "std_error", "90% CI_high")]
    ) %>%
        rename_with(., .fn = ~ gsub("..CI", "%CI", .)) %>%
        mutate(cf.y = cf.y, cf.m = cf.m, rho2 = rho2)
    # fit_short <- data.frame(nu2_Ma1 = nu2[2], nu2_Ma0 = nu2[1], yfit_short =  mus[[1]]$mu_fit3)
    out
}

sens_IIE_Mjo <- function(med, mus, alphas_each, m_each, params, cf.y = 0.04, cf.m = 0.03, rho2 = 1, only_YM_confounded = TRUE, mediatedY = c(1)) {
    # cal_nu2 <- cal.nu2.IIE_Mjo(alphas_each)
    # cal_nu2 <- cal.nu2.IIE_Mjo(alphas_each[mediatedY]) # alpha's bias component is bounded by the plus-signed mean potential outcome
    if (params != "IDE") {
        cal_nu2 <- cal.nu2.IIE_Mjo(alphas_each[ mediatedY ])
    }
    if (params == "IDE") {
        cal_nu2 <- cal.nu2.IDE(alphas_each[ mediatedY ])
    }
    nu2 <- cal_nu2$nu2
    eif_nu2 <- cal_nu2$eif_nu2

    psi_nu2 <- map(1:2, \(y) eif_nu2[[y]] - nu2[[y]])

    # sigma2
    # mus <- mus_IIE_M1
    # m_each <- m_each_IIE_M1

    Y <- med@data[[med@vars@Y]]
    resY <- cal.resY(Y, m_each, mus[mediatedY])
    eif_sigma2 <- resY^2
    # eif_sigma2 <- (Y - mus[[1]]$mu_fit3)^2 # mus[[1]]$mu_fit3 = mus[[2]]$mu_fit3 when m_each[[1]]$ay = m_each[[2]]$ay
    sigma2 <- map(eif_sigma2, \(x) mean(x)) # mean(eif_sigma2)

    psi_sigma2 <- map(1:2, \(y) eif_sigma2[[y]] - sigma2[[y]]) # eif_sigma2 - sigma2

    # S2
    S2 <- map(1:2, \(x) sigma2[[x]] * nu2[[x]])



    # eif
    # eif_short <- eif_each(med, mus_1, alpha_each1, m_each1)
    eif_short <- eif_each_cal(med, mus, alphas_each, m_each)
    names(eif_short) <- names(m_each)

    param <- params

    if (length(m_each) == 2) {
        eif_short_each <- eif_short[, paste0(param, c("_m", "_p"))]
        param_short <- mean(eif_short_each[[2]]) - mean(eif_short_each[[1]])
        psi_param_short <- eif_short_each[[2]] - eif_short_each[[1]] - param_short
    }
    if (length(m_each) == 1) {
        eif_short_each <- eif_short
        param_short <- mean(eif_short_each[[1]])
        psi_param_short <- eif_short_each[[1]] - param_short
    }

    # bias bound

    bf <- bias.factor(cf.y = cf.y, cf.m = cf.m, rho2 = rho2)
    if (length(bf) < 2) {
        if (only_YM_confounded) {
            bf <- c(0, bf)
        } else {
            bf <- c(bf, bf)
        }
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
        S_cluster <- med@data[["S"]]
    }
    se_bounds <- inference.cl(
        data.frame(
            theta_m = psi_param_bound_m + param_bound_m,
            theta_p = psi_param_bound_p + param_bound_p
        ),
        S = S_cluster,
        conf.level = 0.9
    )

    se_theta_s <- inference.cl(
        data.frame(
            theta_s = psi_param_short + param_short
        ),
        S = S_cluster,
        conf.level = 0.95
    )


    out <- data.frame(
        theta_s = se_theta_s,
        theta_m = se_bounds["theta_m", c("estimate", "std_error", "90% CI_low")],
        theta_p = se_bounds["theta_p", c("estimate", "std_error", "90% CI_high")]
    ) %>%
        rename_with(., .fn = ~ gsub("..CI", "%CI", .)) %>%
        mutate(cf.y = cf.y, cf.m = cf.m, rho2 = rho2)
    # fit_short <- data.frame(nu2_Ma1 = nu2[2], nu2_Ma0 = nu2[1], yfit_short =  mus[[1]]$mu_fit3)
    out
}

# Mediation Potential Outcome --------------


sens_Y.Mjo <- function(med, mus, alphas_each, m_each, params, cf.y = 0.04, cf.m = 0.03, rho2 = 1) {
    mediatedY <- grep("am0_ay1", names(alphas_each))
    cal_nu2 <- cal.nu2.IIE_Mjo(alphas_each[mediatedY])

    nu2 <- cal_nu2$nu2
    eif_nu2 <- cal_nu2$eif_nu2

    psi_nu2 <- map(1:length(nu2), \(y) eif_nu2[[y]] - nu2[[y]])

    # sigma2
    # mus <- mus_Mjo[c(3:4)]; m_each <- m_each_Mjo[3:4]
    # med <- med_cwc

    Y <- med@data[[med@vars@Y]]
    resY <- cal.resY(Y, m_each, mus)
    resY <- map(resY, \(x) x[, 1]) # x[, 1] = x[, 2]
    eif_sigma2 <- map(resY, \(x) x^2)
    # eif_sigma2 <- (Y - mus[[1]]$mu_fit3)^2 # mus[[1]]$mu_fit3 = mus[[2]]$mu_fit3 when m_each[[1]]$ay = m_each[[2]]$ay
    sigma2 <- map(eif_sigma2, \(x) mean(x)) # mean(eif_sigma2)

    psi_sigma2 <- map(1:2, \(y) eif_sigma2[[y]] - sigma2[[y]]) # eif_sigma2 - sigma2

    # S2
    S2 <- map(1:2, \(x) sigma2[[x]] * nu2[[x]])



    # eif
    # eif_short <- eif_each(med, mus_1, alpha_each1, m_each1)
    eif_short <- eif_each_cal(med, mus, alphas_each, m_each)
    names(eif_short) <- names(m_each)

    # params = "IDE"
    param <- params

    if (length(m_each) == 2) {
        eif_short_each <- eif_short[, paste0(param, c("_m", "_p"))]
        param_short <- mean(eif_short_each[[2]]) - mean(eif_short_each[[1]])
        psi_param_short <- eif_short_each[[2]] - eif_short_each[[1]] - param_short
    }
    if (length(m_each) == 1) {
        eif_short_each <- eif_short
        param_short <- mean(eif_short_each[[1]])
        psi_param_short <- eif_short_each[[1]] - param_short
    }

    # bias bound

    bf <- bias.factor(cf.y = cf.y, cf.m = cf.m, rho2 = rho2)
    if (length(cf.y) < 2) {
        # only Y-M confounding
        bf <- c(0, bf)
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
        S_cluster <- med@data[["S"]]
    }
    se_bounds <- inference.cl(
        data.frame(
            theta_m = psi_param_bound_m + param_bound_m,
            theta_p = psi_param_bound_p + param_bound_p
        ),
        S = S_cluster,
        conf.level = 0.9
    )

    se_theta_s <- inference.cl(
        data.frame(
            theta_s = psi_param_short + param_short
        ),
        S = S_cluster,
        conf.level = 0.95
    )


    out <- data.frame(
        theta_s = se_theta_s,
        theta_m = se_bounds["theta_m", c("estimate", "std_error", "90% CI_low")],
        theta_p = se_bounds["theta_p", c("estimate", "std_error", "90% CI_high")]
    ) %>%
        rename_with(., .fn = ~ gsub("..CI", "%CI", .)) %>%
        mutate(cf.y = cf.y, cf.m = cf.m, rho2 = rho2)
    # fit_short <- data.frame(nu2_Ma1 = nu2[2], nu2_Ma0 = nu2[1], yfit_short =  mus[[1]]$mu_fit3)
    out
}

# Potential outcomes combined --------------



sens_Y2 <- function(med, mus, alphas_each, m_each, params, cf.y = 0.04, cf.m = 0.03, rho2 = 1) {
    # params <- "IIE_Mjo"
    # alphas_each = alphas_each_IIE_Mjo_r
    # cal_nu2 <- cal.nu2.IIE_M1(alphas_each)
    if (params %in% c("IDE", "IIE_Mjo")) {
        cal_nu2 <- map(1:2, \(x) cal.nu2.IIE_Mjo(alphas_each[x]))
    } else {
        cal_nu2 <- map(1:2, \(x) cal.nu2.IIE_M1(alphas_each[x]))
    }

    nu2 <- map(1:2, ~ unlist(cal_nu2[[.]]$nu2))
    eif_nu2 <- map(1:2, ~ cal_nu2[[.]]$eif_nu2)

    psi_nu2 <- map(1:length(nu2), \(y) map(1:2, \(x) eif_nu2[[x]][[y]] - nu2[[x]][[y]]))

    # sigma2
    # m_each = m_each_Mjo[1:2]
    # mus = mus_Mjo[1:2]
    # mus <- mus_IIE_M1
    # m_each <- m_each_IIE_M1
    # med <- med1_cwc
    Y <- med@data[[med@vars@Y]]
    # map(mus, ~names(.))

    resY <- map(1:2, \(x) cal.resY(Y, m_each[x], mus[x]))

    eif_sigma2 <- map(1:2, \(x) resY[[x]]^2)
    # eif_sigma2 <- (Y - mus[[1]]$mu_fit3)^2 # mus[[1]]$mu_fit3 = mus[[2]]$mu_fit3 when m_each[[1]]$ay = m_each[[2]]$ay
    sigma2 <- map(1:2, \(x) apply(eif_sigma2[[x]], 2, mean))

    psi_sigma2 <- map(1:length(sigma2), \(y) eif_sigma2[[y]] - sigma2[[y]]) # eif_sigma2 - sigma2

    # S2
    S2 <- map(1:length(sigma2), \(x) map_dbl(1:2, \(i) sigma2[[i]][[x]] * nu2[[i]][[x]]))



    # eif
    # eif_short <- eif_each(med, mus_1, alpha_each1, m_each1)
    eif_short <- eif_each_cal(med, mus, alphas_each, m_each)
    names(eif_short) <- names(m_each)

    # params = "IIE_M1"
    param <- params
    if (length(m_each) == 2) {
        eif_short_each <- eif_short[, paste0(param, c("_m", "_p"))]
        param_short <- mean(eif_short_each[[2]]) - mean(eif_short_each[[1]])
        psi_param_short <- eif_short_each[[2]] - eif_short_each[[1]] - param_short
    }
    if (length(m_each) == 1) {
        eif_short_each <- eif_short
        param_short <- mean(eif_short_each[[1]])
        psi_param_short <- eif_short_each[[1]] - param_short
    }

    # bias bound

    bf <- bias.factor(cf.y = cf.y, cf.m = cf.m, rho2 = rho2)

    if (length(cf.y) == 1) {
        bf <- rep(bf, length())
    }

    bias <- map(1:length(S2), \(k) bf[[k]] * (sqrt(S2[[k]][[1]]) + sqrt(S2[[k]][[2]])))

    psi_bias <- map(1:3, \(k)
    (bf[[k]] / 2) * (1 / sqrt(S2[[k]][[1]])) * (sigma2[[k]] * psi_nu2[[k]][[1]] + nu2[[1]][[k]] * psi_sigma2[[k]][[1]]) +
        (bf[[k]] / 2) * (1 / sqrt(S2[[k]][[2]])) * (sigma2[[k]] * psi_nu2[[k]][[2]] + nu2[[2]][[k]] * psi_sigma2[[k]][[1]]))


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
        S_cluster <- med@data[["S"]]
    }
    se_bounds <- inference.cl(
        data.frame(
            theta_m = psi_param_bound_m + param_bound_m,
            theta_p = psi_param_bound_p + param_bound_p
        ),
        S = S_cluster,
        conf.level = 0.9
    )

    se_theta_s <- inference.cl(
        data.frame(
            theta_s = psi_param_short + param_short
        ),
        S = S_cluster,
        conf.level = 0.95
    )


    out <- data.frame(
        theta_s = se_theta_s,
        theta_m = se_bounds["theta_m", c("estimate", "std_error", "90% CI_low")],
        theta_p = se_bounds["theta_p", c("estimate", "std_error", "90% CI_high")]
    ) %>%
        rename_with(., .fn = ~ gsub("..CI", "%CI", .)) %>%
        mutate(cf.y = cf.y, cf.m = cf.m, rho2 = rho2)
    # fit_short <- data.frame(nu2_Ma1 = nu2[2], nu2_Ma0 = nu2[1], yfit_short =  mus[[1]]$mu_fit3)
    out
}

# Bias factor sensitivity parameter ------------
bias.factor <- function(cf.y = 0.03, cf.m = 0.04, rho2 = 1) {
    sqrt(rho2 * cf.y * (cf.m / (1 - cf.m)))
}




rv.value <- function(med, mus, alphas_each, m_each, params, theta.null = 0, cf = seq(0, 0.99, by = 0.05), only_YM_confounded = TRUE, mediatedY = c(1)) {
    # cf <- seq(0, 0.99,by = 0.05)
    sens_fun <- ifelse(m_each[[1]][["az"]] == "", sens_IIE_Mjo, sens_IIE_M1)

    values <- map2(
        .x = cf, .y = cf,
        .f = ~ sens_fun(med, mus, alphas_each, m_each, params, cf.y = .x, cf.m = .y, rho2 = 1, only_YM_confounded, mediatedY)
    )
    # values <- map2(.x = cf, .y = cf, .f = ~ sens_fun(med, mus_IIE_M1, alphas_each=alphas_each_IIE_M1_r, m_each=m_each_IIE_M1, params = "IIE_M1", cf.y = .x, cf.m = .y, rho2 = 1) )

    df_values <- do.call(rbind, values)

    df_values1 <- df_values %>%
        rename(conf_bound.l = `theta_m.90%CI_low`, conf_bound.u = `theta_p.90%CI_high`) %>%
        mutate(
            # Rsq = cf,
            cf = cf, theta.null = theta.null,
            if_cover = (conf_bound.l - theta.null) * (conf_bound.u - theta.null) < 0
        )

    if (df_values1[df_values1$cf.m == 0 & df_values1$cf.y == 0, "if_cover"] == TRUE) {
        rv <- min(df_values1$cf[df_values1$if_cover == FALSE])
    } else {
        rv <- min(df_values1$cf[df_values1$if_cover == TRUE])
    }
    

    rv_res <- list(rv, df_values1)
    names(rv_res) <- c(glue("RV_null{theta.null}_sig{0.05}"), glue("RVsens_res"))

    rv_res
}


# NOT USED --------------------


sens_YMA <- function(med, mus, alphas_each, m_each, params, cf.y = 0.04, cf.m = 0.03, rho2 = 1) {
    cal_nu2 <- cal.nu2(alphas_each)
    nu2 <- cal_nu2$nu2
    eif_nu2 <- cal_nu2$eif_nu2

    psi_nu2 <- map(1:3, \(y) sapply(1:2, \(x) eif_nu2[[y]][, x] - nu2[[y]][x]))

    # sigma2
    # mus <- mus_IIE_M1
    # m_each <- m_each_IIE_M1

    Y <- med@data[[med@vars@Y]]
    resY <- cal.resY(Y, m_each, mus)
    eif_sigma2 <- map(resY, \(x) x^2)
    # eif_sigma2 <- (Y - mus[[1]]$mu_fit3)^2 # mus[[1]]$mu_fit3 = mus[[2]]$mu_fit3 when m_each[[1]]$ay = m_each[[2]]$ay
    sigma2 <- map(eif_sigma2, \(x) colMeans(x)) # mean(eif_sigma2)

    psi_sigma2 <- map(1:3, \(y) sapply(1:2, \(x) eif_sigma2[[y]][, x] - sigma2[[y]][x])) # eif_sigma2 - sigma2

    # S2
    S2 <- map(1:3, \(x) sigma2[[x]] * nu2[[x]])
    S2_sum <- reduce(S2, `+`)
    psi_S2 <- map(1:3, \(y) sapply(1:2, \(x) sigma2[[y]][x] * psi_nu2[[y]][, x] + nu2[[y]][x] * psi_sigma2[[y]][, x]))
    # psi_S2[[1]][1:3, ] + psi_S2[[2]][1:3, ]
    psi_S2_sum <- reduce(psi_S2, `+`)

    # eif
    # eif_short <- eif_each(med, mus_1, alpha_each1, m_each1)
    eif_short <- eif_each_cal(med, mus, alphas_each, m_each)
    names(eif_short) <- names(m_each)

    # one.param ----
    one.param <- function(param = "IIE_M1", cf.y = 0.03, cf.m = 0.03, rho2 = 1) {
        eif_short_each <- eif_short[, paste0(param, c("_m", "_p"))]
        param_short <- mean(eif_short_each[[2]]) - mean(eif_short_each[[1]])
        psi_param_short <- eif_short_each[[2]] - eif_short_each[[1]]

        # bias bound
        bf <- bias.factor(cf.y = cf.y, cf.m = cf.m, rho2 = rho2)
        bias_bound <- sqrt(S2_sum) * bf
        # bounds IFs
        psi_bias_bound <- lapply(1:2, \(x) (bf / 2) * (1 / sqrt(S2_sum[x])) * psi_S2_sum[, x])

        # thetas_short <- colMeans(eif_short_each)
        # theta_short <- mean(eif_short)
        # psi_theta_short <- eif_short - theta_short
        # currently, focusing on individual average effect
        param_bound_p <- (mean(eif_short_each[[2]]) + bias_bound[2]) - (mean(eif_short_each[[1]]) + bias_bound[1])
        param_bound_m <- (mean(eif_short_each[[2]]) - bias_bound[2]) - (mean(eif_short_each[[1]]) - bias_bound[1])


        psi_param_bound_p <- (eif_short_each[[2]] + psi_bias_bound[[2]]) - (eif_short_each[[1]] + psi_bias_bound[[1]])
        psi_param_bound_m <- (eif_short_each[[2]] - psi_bias_bound[[2]]) - (eif_short_each[[1]] - psi_bias_bound[[1]])


        # sens_res <- data.frame(
        #     param_short = param_short,
        #     param_bound_m = param_bound_m,
        #     param_bound_p = param_bound_p,
        #     psi_param_short = eif_short_each[[2]] - eif_short_each[[1]],
        #     psi_param_bound_m = psi_param_bound_m,
        #     psi_param_bound_p = psi_param_bound_p #,
        #     # psi_bias_bound_Ma = psi_bias_bound[,1],
        #     # bias_bound = bias_bound
        # )
        # ci_bounds <- inference.cl(sens_res[, c("psi_param_bound_m", "psi_param_bound_p")], med@data[["S"]], average = "individual", conf.level = 0.9) # one-sided
        # sens_effect1 <- sens_res %>%
        #     mutate(
        #         ci_long_l = ci_bounds["psi_param_bound_m", "ci_low"],
        #         ci_long_u = ci_bounds["psi_param_bound_p", "ci_high"]
        #     ) %>%
        #     summarise(across(everything(), ~mean(.))) %>%
        #     mutate(cf.y = cf.y, cf.m = cf.m, rho2 = rho2)
        #
        # sens_effect2 <- sens_effect1 %>%
        #     # mutate(cf.y, cf.m, rho2) %>%
        #     select(param_short, param_bound_m, param_bound_p, ci_long_l, ci_long_u, cf.y, cf.m, rho2)
        se_bounds <- inference.cl(
            data.frame(
                theta_m = psi_param_bound_m + param_bound_m,
                theta_p = psi_param_bound_p + param_bound_p
            ),
            S = med@data[["S"]], average = "individual",
            conf.level = 0.9
        )

        se_theta_s <- inference.cl(
            data.frame(
                theta_s = psi_param_short + param_short
            ),
            S = med@data[["S"]], average = "individual",
            conf.level = 0.95
        )


        out <- data.frame(
            theta_s = se_theta_s,
            theta_m = se_bounds["theta_m", c("estimate", "std_error", "90% CI_low")],
            theta_p = se_bounds["theta_p", c("estimate", "std_error", "90% CI_high")]
        ) %>%
            rename_with(., .fn = ~ gsub("..CI", "%CI", .)) %>%
            mutate(cf.y = cf.y, cf.m = cf.m, rho2 = rho2)



        # names(sens_effect2) <- c("estimate", "bound.l", "bound.u", "conf_bound.l", "conf_bound.u", "cf.y", "cf.m", "rho2")

        # sens_effect2
        out
    }
    # fit_short <- data.frame(nu2_Ma1 = nu2[2], nu2_Ma0 = nu2[1], yfit_short =  mus[[1]]$mu_fit3)
    sens_out <- lapply(params, \(param) one.param(param, cf.y, cf.m, rho2))
    names(sens_out) <- params
    # list(sens_res = sens_res, fit_short = fit_short)
    sens_out
}

sens_cal <- function(med, mus, alphas_each, m_each, cf.y = 0.04, cf.m = 0.03, rho2 = 1, one_sided.level = 0.95) {
    # sens_list <- lapply(1:length(m_each), \(x=1) sens_each(med, mus[[x]], alphas_each[[x]], m_each[[x]], cf.y = cf.y, cf.m = cf.m, rho2 = rho2) )
    # sens_effect <- sens_list[[2]] - sens_list[[1]]
    sens_fit_res <- sens_YM(med, mus, alphas_each, m_each, cf.y, cf.m, rho2)
    sens_res <- sens_fit_res[["sens_res"]]
    n <- nrow(sens_res)
    ci_bounds <- inference.cl(sens_res[, c("psi_IIE_bound_m", "psi_IIE_bound_p")], med@data[[med@vars@cluster_id]], average = "individual", conf.level = 1 - (1 - one_sided.level) * 2) # one-sided
    sens_effect1 <- sens_res %>%
        mutate(
            ci_long_l = ci_bounds["psi_IIE_bound_m", "ci_low"],
            ci_long_u = ci_bounds["psi_IIE_bound_p", "ci_high"]
        ) %>%
        summarise(across(everything(), ~ mean(.))) %>%
        mutate(cf.y = cf.y, cf.m = cf.m, rho2 = rho2)

    sens_effect2 <- sens_effect1 %>%
        # mutate(cf.y, cf.m, rho2) %>%
        select(IIE_short, IIE_bound_m, IIE_bound_p, ci_long_l, ci_long_u, cf.y, cf.m, rho2)

    names(sens_effect2) <- c("estimate", "bound.l", "bound.u", "conf_bound.l", "conf_bound.u", "cf.y", "cf.m", "rho2")

    list(
        sens_effect = sens_effect2,
        sens_res = sens_fit_res[["sens_res"]],
        fit_short = sens_fit_res[["fit_short"]]
    )
}
