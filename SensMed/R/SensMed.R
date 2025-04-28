#' @export
#'
Sens.Med <- function(
        data,
        treatment,
        outcome,
        mediators,
        covariates,
        a1 = 1, a0 = 0,
        id,
        cluster_id = NA_character_,
        weights = rep(1, nrow(data)),
        learners_outcome = c("glm", "xgboost"),
        learners_treatment = c("glm", "mean", "lightgbm"),
        cf.y = 0.13, cf.m = 0.13, rho2 = 1, conf.level = 0.95,
        control = estimate_control(crossfit_folds = 1)
) {
    data_in <- data
    
    if (!is.na(cluster_id)) {
        # dummy indicators
        Sdumm <- dummy_cols(data[[cluster_id]], remove_first_dummy = FALSE, remove_selected_columns = TRUE)
        colnames(Sdumm) <- paste0("S", 1:ncol(Sdumm))
        Sname_dummies <- colnames(Sdumm)
        
        # cluster means
        data <- data %>%
            mutate(S = data_in[[cluster_id]]) %>%
            group_by(S) %>%
            mutate(across(
                all_of(c(treatment, mediators, outcome, covariates)),
                list(clmean = ~ mean(.), cwc = ~ . - mean(.))
            )) %>%
            bind_cols(Sdumm) %>%
            ungroup()
        
        med <- med_data(
            data = data,
            vars = med_vars(
                A = treatment,
                Y = outcome,
                M = mediators,
                # Z = moc %??% NA_character_,
                C = covariates,
                # C = c(glue("{covariates}_cwc"), Sname_dummies),
                cluster_id = cluster_id,
                Wj = NA_character_,
                S = "S",
                id = id
            ),
            weights = rep(1, nrow(data)),
            a0 = a0,
            a1 = a1
        )
    }
    
    if (is.na(cluster_id)) {
        med <- med_data(
            data = data,
            vars = med_vars(
                A = treatment,
                Y = outcome,
                M = mediators,
                # Z = moc %??% NA_character_,
                C = covariates,
                # C = c(glue("{covariates}_cwc"), Sname_dummies),
                cluster_id = NA_character_,
                Wj = NA_character_,
                S = NA_character_,
                id = id
            ),
            weights = rep(1, nrow(data)),
            a0 = a0,
            a1 = a1
        )
    }
    
    
    # ?crumble::crumble
    dat_crumble <- data %>% ungroup() 
    
    Mname <- mediators
    crumble_res <- crumble::crumble(
        data = dat_crumble,
        trt = treatment,
        outcome = outcome, mediators = Mname[[2]], moc = Mname[[1]],
        covar = covariates,
        d0 = \(data, trt) rep(0, nrow(data)),
        d1 = \(data, trt) rep(1, nrow(data)),
        effect = c("RI"),
        learners = learners_outcome,
        nn_module = crumble::sequential_module(dropout = 0.05),
        id = "id",
        control = crumble::crumble_control(crossfit_folds = 1, mlr3superlearner_folds = 1, device = "cpu", epochs = 10)
    )
    
    # library(doFuture)
    # HDmed_res <- HDmediation::mediation(
    #   data = data,
    #   A = treatment, W = covariates, Z = mediators[2], M = mediators[1],
    #   Y = outcome, family = "gaussian",
    #   folds = 1, partial_tmle = F,
    #   learners_g = c("SL.glm", "SL.ranger", "SL.mean"),
    #   learners_e = c("SL.glm", "SL.ranger", "SL.mean"),
    #   learners_c = c("SL.glm", "SL.ranger", "SL.mean"),
    #   learners_b = c("SL.glm", "SL.ranger", "SL.mean"),
    #   learners_hz = c("SL.glm", "SL.ranger", "SL.nnet"),
    #   learners_u = c("SL.glm", "SL.ranger", "SL.mean"),
    #   learners_ubar = c("SL.glm", "SL.ranger", "SL.mean"),
    #   learners_v = c("SL.glm", "SL.ranger", "SL.mean"),
    #   learners_vbar = c("SL.glm", "SL.ranger", "SL.mean")
    # )
    
    n <- nrow(data)
    mus_IE <- list() 
    mus_DE <- list() 
    
    mus_IE[["az_am0_ay1"]] <- data.frame(
        mu_fit_ay_jo =  crumble_res$outcome_reg$theta_n$natural$fit3_natural[, "100"],
        mu_fit_am_jo = crumble_res$outcome_reg$theta_n$natural$fit2_natural[, "100"],
        b_ay_jo = crumble_res$outcome_reg$theta_n$bs$b3[, "100"],
        b_am_jo = crumble_res$outcome_reg$theta_n$bs$b2[, "100"])
    mus_IE[["az_am1_ay1"]] <- data.frame(
        mu_fit_ay_jo =  crumble_res$outcome_reg$theta_n$natural$fit3_natural[, "111"],
        mu_fit_am_jo = crumble_res$outcome_reg$theta_n$natural$fit2_natural[, "111"],
        b_ay_jo = crumble_res$outcome_reg$theta_n$bs$b3[, "111"],
        b_am_jo = crumble_res$outcome_reg$theta_n$bs$b2[, "111"])
    
    mus_DE[["az_am0_ay0"]] <- data.frame(
        mu_fit_ay_jo =  crumble_res$outcome_reg$theta_n$natural$fit3_natural[, "000"],
        mu_fit_am_jo = crumble_res$outcome_reg$theta_n$natural$fit2_natural[, "000"],
        b_ay_jo = crumble_res$outcome_reg$theta_n$bs$b3[, "000"],
        b_am_jo = crumble_res$outcome_reg$theta_n$bs$b2[, "000"])
    mus_DE[["az_am0_ay1"]] <- mus_IE[["az_am0_ay1"]]
    
    
    alphas_IE <- list()
    alphas_DE <- list()
    
    alphas_IE[["az_am0_ay1"]] <- data.frame(
        alpha1 =  crumble_res$alpha_n$alpha1[, "100"],
        alpha2 =  crumble_res$alpha_n$alpha2[, "100"],
        alpha3 =  crumble_res$alpha_n$alpha3[, "100"])
    alphas_IE[["az_am1_ay1"]] <- data.frame(
        alpha1 =  crumble_res$alpha_n$alpha1[, "111"],
        alpha2 =  crumble_res$alpha_n$alpha2[, "111"],
        alpha3 =  crumble_res$alpha_n$alpha3[, "111"])
    
    
    # EIF -----
    m_each_IE <- list(IE_m = c(az = "", am = "data_0", ay = "data_1"),
                   IE_p = c(az = "", am = "data_1", ay = "data_1"))
    eifs_IE <- eif_each_cal(med, mus_IE, alphas_IE, m_each_IE)
    eifs_IE[["IE"]] <- eifs_IE[[2]] - eifs_IE[[1]]
    
    eifs <- eifs_IE
    if (is.null(cluster_id) | is.na(cluster_id)) {
        estimates <- inference.iid(eifs)
    }
    if (!is.null(cluster_id) & !is.na(cluster_id)) {
        data$S <- data_in[[cluster_id]]
        estimates <- inference.cl(eifs, data$S, average = "individual", conf.level = 0.95)
    }
    
    eifs_IIE_M1 <- eif_each_cal(med1, mus_IIE_M1, alphas_each_IIE_M1_r, m_each_IIE_M1)
    eifs_IIE_M1[["IIE_M1"]] <- eifs_IIE_M1[[2]] - eifs_IIE_M1[[1]]
    
    eifs_IIE_M2 <- eif_each_cal(med2_cwc, mus_IIE_M2, alphas_each_IIE_M2_r, m_each_IIE_M2)
    eifs_IIE_M2[["IIE_M2"]] <- eifs_IIE_M2[[2]] - eifs_IIE_M2[[1]]
    
    eifs_IDE <- eif_each_cal(med_cwc, mus_Mjo[c(3:4)], alphas_each_IDE_r, m_each_Mjo[3:4])
    eifs_IDE[["IDE"]] <- eifs_IDE[[2]] - eifs_IDE[[1]]
    
    
    eifs_IIE_Mjo <- eif_each_cal(med_cwc, mus_Mjo[c(1:2)], alphas_each_IIE_Mjo_r, m_each_Mjo[1:2])
    eifs_IIE_Mjo[["IIE_Mjo"]] <- eifs_IIE_Mjo[[2]] - eifs_IIE_Mjo[[1]]
    
    eifs_IIE_Mjo[["IIE_mu"]] <- eifs_IIE_Mjo[["IIE_Mjo"]] - eifs_IIE_M1[["IIE_M1"]] - eifs_IIE_M2[["IIE_M2"]]
    # eifs_IDE <- eif_each_cal(med, mus_Mjo[c(3:4)], alphas_each_IDE_r, m_each_Mjo[3:4])
    #     eifs_IDE[["IDE"]] <- eifs_IDE[[2]] - eifs_IDE[[1]]
    #
    #     eifs_IIE_Mjo <- eif_each_cal(med, mus_Mjo[c(1:2)], alphas_each_IIE_Mjo_r, m_each_Mjo[1:2])
    #     eifs_IIE_Mjo[["IIE_Mjo"]] <- eifs_IIE_Mjo[[2]] - eifs_IIE_Mjo[[1]]
    
    # eifs <- c(eifs_IDE, eifs_IIE_Mjo) %>% do.call(cbind, .) %>% data.frame()
    
    eifs <- c(eifs_IIE_M1, eifs_IIE_M2, eifs_IDE, eifs_IIE_Mjo) %>%
        do.call(cbind, .) %>%
        data.frame()
    
    if (is.null(cluster_id) | is.na(cluster_id)) {
        estimates <- inference.iid(eifs)
    }
    if (!is.null(cluster_id) & !is.na(cluster_id)) {
        data$S <- data_in[[cluster_id]]
        estimates <- inference.cl(eifs, data$S, average = "individual", conf.level = 0.95)
    }
    
    # estimates_C <- inference.cl(eifs_each, data_in[["S"]], average = "cluster")
    estimates <- estimates %>%
        rownames_to_column(var = "estimand") %>%
        arrange(estimand)
    # estimates_C <- inference.cl(eifs_each, data_in[["S"]], average = "cluster")
    
    # Sensitivity analysis ------------
    
    
    ## RV values ----
    sens_fit_res <- list()
    rv_res <- list()
    # IIE_M1
    sens_fit_res[["IIE_M1"]] <- sens_IIE_M1(med_cwc, mus_IIE_M1, alphas_each_IIE_M1_r, m_each_IIE_M1, params = "IIE_M1", cf.y, cf.m, rho2)
    
    rv_res[["IIE_M1"]] <- rv.value(med_cwc, mus_IIE_M1, alphas_each_IIE_M1_r, m_each_IIE_M1, params = "IIE_M1", theta.null = 0, cf = seq(0, 0.99, by = 0.01))
    
    # IIE_M2
    sens_fit_res[["IIE_M2"]] <- sens_IIE_M1(med_cwc, mus_IIE_M2, alphas_each_IIE_M2_r, m_each_IIE_M2, params = "IIE_M2", cf.y, cf.m, rho2)
    
    rv_res[["IIE_M2"]] <- rv.value(med_cwc, mus_IIE_M2, alphas_each_IIE_M2_r, m_each_IIE_M2, params = "IIE_M2", theta.null = 0, cf = seq(0, 0.99, by = 0.01))
    
    # IDE
    sens_fit_res[["IDE"]] <- sens_Y.M(med_cwc, mus_Mjo[3:4], alphas_each_IDE_r, m_each_Mjo[3:4], params = "IDE", cf.y, cf.m, rho2)
    
    rv_res[["IDE"]] <- rv.value(med_cwc, mus_Mjo[3:4], alphas_each_IDE_r, m_each_Mjo[3:4], params = "IDE", theta.null = 0, cf = seq(0, 0.99, by = 0.01))
    
    # IIE_Mjo
    sens_fit_res[["IIE_Mjo"]] <- sens_IIE_Mjo(med_cwc, mus_Mjo[1:2], alphas_each_IIE_Mjo_r, m_each_Mjo[1:2], params = "IIE_Mjo", cf.y, cf.m, rho2)
    
    rv_res[["IIE_Mjo"]] <- rv.value(med_cwc, mus_Mjo[1:2], alphas_each_IIE_Mjo_r, m_each_Mjo[1:2], params = "IIE_Mjo", theta.null = 0, cf = seq(0, 0.99, by = 0.01))
    
    ## Minimal sensitivity reporting ----
    rv_values <- map_dbl(rv_res, \(x) x$RV_null0_sig0.05)
    tab_sensitivity_reporting <- estimates %>%
        dplyr::filter(estimand %in% c("IDE", "IIE_M1", "IIE_M2")) %>%
        mutate(CI = glue("[{round(`95% CI_low`, 2)}, {round(`95% CI_high`, 2)}]")) %>% 
        mutate(RV_null0_sig0.05 = rv_values[estimand])
    
    tab_sensitivity_reporting %>% write_csv(file = "Example/Results/tab_sensitivity_reporting.csv")
    
    ## Contour plots -----
    sens_grid <- list()
    
    cf_grid <- expand.grid(cf.y = seq(0.01, 0.99, by = 0.01), cf.m = seq(0.01, 0.99, by = 0.01))
    
    ### IIE_M1 contour ----
    sens_grid[["IIE_M1"]] <- map2_df(
        .x = cf_grid$cf.m, .y = cf_grid$cf.y,
        .f = ~ sens_IIE_M1(med, mus_IIE_M1, alphas_each_IIE_M1_r, m_each_IIE_M1,
                           params = "IIE_M1",
                           cf.y = .y, cf.m = .x, rho2 = 1
        )
    )
    
    
    z_grid <- sens_grid[["IIE_M1"]] %>% 
        select(cf.y, cf.m, `theta_m.90%CI_low`) %>% 
        pivot_wider(names_from = cf.m, values_from = `theta_m.90%CI_low`) %>% 
        as.matrix(.)
    z_grid <- z_grid[, -1]
    
    
    
    pdf(file = "Example/Results/contour_IIE_M1.pdf", width = 5, height = 5)
    contour_plot(grid_values.x = unique(cf_grid$cf.m),
                 grid_values.y = unique(cf_grid$cf.y),
                 z_axis = z_grid,
                 levels = c(-3.5, -2, -1.5, -1, -0.8, -0.5, -0.3, -0.1, 0, 0.1, 0.2),
                 labels = NULL,
                 threshold = 0,
                 col.contour = "blue4")
    
    dev.off()
    
    
    ### IDE contour -----
    sens_grid[["IDE"]] <- map2_df(
        .x = cf_grid$cf.m, .y = cf_grid$cf.y,
        .f = ~ sens_IIE_Mjo(med, mus_Mjo[3:4], alphas_each_IDE_r, m_each_Mjo[3:4],
                            params = "IDE",
                            cf.y = .y, cf.m = .x, rho2 = 1
        )
    )
    
    
    z_grid <- sens_grid[["IDE"]] %>% 
        select(cf.y, cf.m, `theta_m.90%CI_low`) %>% 
        pivot_wider(names_from = cf.m, values_from = `theta_m.90%CI_low`) %>% 
        as.matrix(.)
    z_grid <- z_grid[, -1]
    
    
    quantile(z_grid)
    pdf(file = "Example/Results/contour_IDE.pdf", width = 5, height = 5)
    contour_plot(grid_values.x = unique(cf_grid$cf.m),
                 grid_values.y = unique(cf_grid$cf.y),
                 z_axis = z_grid,
                 levels = c( -4,  -2, -1.5, -1,  -0.5, 0, 0.3, 0.6, 1, 1.3, 1.4),
                 labels = NULL,
                 threshold = 0,
                 col.contour = "blue4")
    
    dev.off()
    
    
    # Benchmark ------------------------
    Benchmark_res <- list()
    
    Benchmark_res[["IIE_M1"]] <- benchmark.IIE_M1(
        med1_cwc,
        benchmark_covariates = c("X_BYSES1"), param = "IIE.M1", mus_IIE_M1, alphas_each_IIE_M1_r, m_each_IIE_M1, estimates, folds, control
    )
    
    med_sens <- list(estimates, rv_res, sens_fit_res)
    names(med_sens) <- c("estimates_short", "robustness_value", "sens_fit_res")
    out <- mget(ls(envir = environment()))
    
    out
    # med_sens
}
