

mediate <- function(
    data,
    treatment, 
    outcome,
    mediators,
    covariates,
    a1 = 1, a0 = 0,
    id,
    cluster_id,
    weights = rep(1, nrow(data)),
    learners = c("glm", "nnet"),
    cf.y = 0.13, cf.m = 0.13, rho2 = 1, conf.level = 0.95,
    control = estimate_control()
) {
    data_in <- data
    Sname <- cluster_id
    # AM <- data_in[[Aname]] * data_in[, Mnames]
    # colnames(AM) <- c("AM.M1", "AM.M2")
    # Mint <- data.frame(data_in[, Mnames[1]] * data_in[, Mnames[2]])
    # colnames(Mint) <- "Mint"
    # AMint <- data.frame(data_in[[Aname]] * data_in[, Mnames[1]] * data_in[, Mnames[2]])
    # colnames(AMint) <- "AMint"
    
    # dummy indicators
    Sdumm <- dummy_cols(data[[Sname]], remove_first_dummy = FALSE, remove_selected_columns = TRUE)
    colnames(Sdumm) <- paste0("S", 1:ncol(Sdumm))
    Sname_dummies <- colnames(Sdumm)
    
    # cluster means
    data <- data %>% # data.frame(data, AM, Mint, AMint) %>%
        rename(A = all_of(treatment), S = all_of(cluster_id), Y = all_of(outcome), M = all_of(mediators)) %>% 
        # mutate(A0 = 1 - A) %>% 
        group_by(S) %>%
        mutate(across(c(!id), list(clmean = ~mean(.), cwc = ~.-mean(.)))) %>%
        bind_cols(Sdumm)
    
    data[["S"]] <- match(data[["S"]], unique(data[["S"]])) %>% as.factor()
    
    quantile(data$A_clmean)
    # cluster level covariates
    Wj <- data %>% group_by(S) %>% summarise(across(all_of(covariates), sd)) %>% 
        ungroup() %>% summarise(across(all_of(covariates), ~max(.))) %>% 
        select(where(~.==0)) %>% 
        names(.)
    
    covariates_c <- c(Wj, Sname_dummies)
    # individual level covariates
    covariates_i <- setdiff(covariates, Wj)
    covariates_cl <- list(Wj = Wj, 
                          # Sdumm = NULL,
                          Sdumm = Sname_dummies, 
                          Xij = covariates_i, covariates_c = covariates_c)
    # names(covariates_cl)
    varnames <- list("A" = "A", "M" = glue("M{1:length(mediators)}"), "Y" = "Y",
                     # "AM" = paste0("AM.", Mnames),
                     "S" = "S", "Sdumm" = Sname_dummies,
                     "X" = covariates, "W" = Wj)
    
    
    med <- med_data(
        data = data,
        vars = med_vars(
            A = "A",
            Y = "Y",
            M = glue("M{1:length(mediators)}"),
            # Z = moc %??% NA_character_,
            C = covariates,
            # C = c(glue("{covariates}_cwc"), Sname_dummies),
            cluster_id = "S",
            id = id 
        ),
        weights = weights,
        a0 = a0,
        a1 = a1
    )
    
    # Multiple mediators
    # Create independent M -----
    # this estimates mediation for the first mediator M1
    if (length(med@vars@M) >= 2) {
        m_each_Mjo <- list(
            IIE_Mjo_m = c(az = "", am = "data_0", ay = "data_1"),
            IIE_Mjo_p = c(az = "", am = "data_1", ay = "data_1"),
            IDE_m = c(az = "", am = "data_0", ay = "data_0"),
            IDE_p = c(az = "", am = "data_0", ay = "data_1")
        )
        # IIE_M1, independently draw M2 from its marginal distribution
        med <- add_indepM(med, control, which.M = 2, cwc = TRUE)
        med1 <- med
        # med <- add_indepM.within(med, control, which.M = which.M, if_indepM = TRUE)
        m_each_IIE_M1 <- list(
            IIE_M1_m = c(az = "data_0_indepM", am = "data_0_indepM", ay = "data_1_indepM") 
            ,
            IIE_M1_p = c(az = "data_0_indepM", am = "data_1_indepM", ay = "data_1_indepM"))
        
        # IIE_M2, independently draw M1 from its marginal distribution
        med2 <- add_indepM(med, control, which.M = 1, cwc = TRUE)
        m_each_IIE_M2 <- list(
            IIE_M2_m = c(az = "data_1_indepM", am = "data_0_indepM", ay = "data_1_indepM") 
            ,
            IIE_M2_p = c(az = "data_1_indepM", am = "data_1_indepM", ay = "data_1_indepM"))

    }
    
    # folds <- make.fold_K(data, Snames = "S", control$crossfit_folds)
    folds <- origami::make_folds(med@data, V = control$crossfit_folds)
    if (control$crossfit_folds == 1) {
        folds[[1]]$training_set <- folds[[1]]$validation_set
    }
    
    # Fit models ------------------------
    ## IIE_Mjo and IDE -----------------------
    # Estimate mu 
    # med_cwc <- med
    # med_cwc@vars <- med_vars(A = "A_cwc", 
    #                          Y = "Y_cwc", 
    #                          M = glue("M{1:length(mediators)}_cwc"), 
    #                          C = glue("{covariates}_cwc"), 
    #                          cluster_id = "S", id = id)
    # 
    # med_cwcSdumm <- med_cwc
    # med_cwcSdumm@vars@C <- c(glue("{covariates}_cwc"), Sname_dummies)
    
    #  mus <- estimate_mu_cl(med, m_each, folds, learners, control, which.M, covariates_cl)
    mus_Mjo <- estimate_mu_cl(med, m_each_Mjo, folds, learners, control, which.M = NULL)
  
    # mean(mus_Mjo[[2]]$b_am_jo - mus_Mjo[[1]]$b_am_jo)
    # mean(mus_Mjo[[4]]$b_am_jo - mus_Mjo[[3]]$b_am_jo)

    # Estimate Riesz representer alphas
    # a_c <- a.c(med, varnames, cluster_opt = "noncluster.mlr", 
    #            # cluster_opt = "cwc", 
    #            folds, glue("SL.{learners}"), bounded = TRUE)
    a_c <- a.c(med, folds, learners, control)
    
    # learners <- c("glm", "xgboost")
    # a_mc <- a.mc(med, # varnames, cluster_opt = "noncluster.mlr", bounded = TRUE,
    #              # cluster_opt = "cwc", 
    #              folds, glue("SL.{learners}"), )
    a_mc <- a.mc(med, folds, learners, control)
    
    
    alpha_ac <- alpha.ac(A = med@data$A, a_c)
    # mean((alpha_ac$`A=0`) * data$Y)
    
    # control$learners_riesz <- list("nn" = list(.key = "nn", verbose = FALSE, interactions = 1, learning_rate = 0.001, hidden = 30, layers = 1))
    
    alphas_each_IIE_Mjo <- riesz_crossfit(med, m_each = m_each_Mjo[c(1:2)], folds, which.M = NULL, learners_riesz = control$learners_riesz, control, riesz_alpha1 = F, riesz_alpha3 = F, alpha_ac, a_c, a_mc)
    
    alphas_each_IDE <- riesz_crossfit(med, m_each = m_each_Mjo[c(3:4)], folds, which.M = NULL, learners_riesz = control$learners_riesz, control, riesz_alpha1 = F, riesz_alpha3 = F, alpha_ac, a_c, a_mc)
    
    # mean(alphas_each_IDE[[2]]$alpha3 * data$Y) - mean(alphas_each_IDE[[1]]$alpha3 * data$Y)
    # mean(alphas_each_IIE_Mjo[[2]]$alpha3 * data$Y) - mean(alphas_each_IIE_Mjo[[1]]$alpha3 * data$Y)
    
    ## IIE_M1 --------------
    # med1_cwc <- med1
    # med1_cwc@vars <- med_vars(A = "A_cwc", 
    #                          Y = "Y_cwc", 
    #                          M = glue("M{1:length(mediators)}_cwc"), 
    #                          C = glue("{covariates}_cwc"), 
    #                          cluster_id = "S", id = id)
    # med1_cwcSdumm <- med1_cwc
    # med1_cwcSdumm@vars@C <- c(glue("{covariates}_cwc"), Sname_dummies)
    mus_IIE_M1 <- estimate_mu_cl(med1, m_each_IIE_M1, folds, learners, control, which.M = 2)
    
    # mean(mus_IIE_M1[[2]]$b1 - mus_IIE_M1[[1]]$b1)
    
    # nn_module <- sequential_module()
    # alphas_each <- estimate_alpha_cl(med, m_each, folds, nn_module, control, which.M, covariates_cl)
    # alphas_each <- estimate_alpha_each(med, m_each[c(5:6)], folds, nn_module, control, which.M)
    
    # aa=lrn("riesz.nn"); aa$param_set; mlr_learners
    control$learners_riesz <- list("nn" = list(.key = "nn", verbose = FALSE, interactions = 1, learning_rate = 0.001/2, hidden = 30, layers = 1))
    
    alphas_each_IIE_M1 <- riesz_crossfit(med1, m_each = m_each_IIE_M1, folds, which.M = 2, learners_riesz = control$learners_riesz, control, riesz_alpha1 = F, riesz_alpha3 = TRUE, alpha_ac, a_c, a_mc)
    
    # mean(alphas_each_IIE_M1[[2]]$alpha3 * data$Y) - mean( alphas_each_IIE_M1[[1]]$alpha3 * data$Y)
    
    
    ## IIE_M2 --------------
    # med2_cwc <- med2
    # med2_cwc@vars <- med_vars(A = "A_cwc", 
    #                           Y = "Y_cwc", 
    #                           M = glue("M{1:length(mediators)}_cwc"), 
    #                           C = glue("{covariates}_cwc"), 
    #                           cluster_id = "S", id = id)
    # med2_cwcSdumm <- med2_cwc
    # med2_cwcSdumm@vars@C <- c(glue("{covariates}_cwc"), Sname_dummies)
    
    #  mus <- estimate_mu_cl(med, m_each, folds, learners, control, which.M, covariates_cl)
    mus_IIE_M2 <- estimate_mu_cl(med2, m_each_IIE_M2, folds, learners, control, which.M = 1)
    # mean(mus_IIE_M2[[2]]$b1 - mus_IIE_M2[[1]]$b1)
    
    control$learners_riesz <- list("nn" = list(.key = "nn", verbose = FALSE, interactions = 1, learning_rate = 0.001, hidden = 30, layers = 1))
    
    alphas_each_IIE_M2 <- riesz_crossfit(med2, m_each = m_each_IIE_M2, folds, which.M = 1, learners_riesz = control$learners_riesz, control, riesz_alpha1 = F, riesz_alpha3 = TRUE, alpha_ac, a_c, a_mc)
    
    # mean(alphas_each_IIE_M2[[2]]$alpha3 * data$Y) - mean( alphas_each_IIE_M2[[1]]$alpha3 * data$Y)
    
    # EIF -----
    
    eifs_IDE <- eif_each_cal(med, mus_Mjo[c(3:4)], alphas_each_IDE, m_each_Mjo[3:4])
    eifs_IDE[["IDE"]] <- eifs_IDE[[2]] - eifs_IDE[[1]] 
    
    eifs_IIE_Mjo <- eif_each_cal(med, mus_Mjo[c(1:2)], alphas_each_IIE_Mjo, m_each_Mjo[1:2])
    eifs_IIE_Mjo[["IIE_Mjo"]] <- eifs_IIE_Mjo[[2]] - eifs_IIE_Mjo[[1]] 
    
    eifs_IIE_M1 <- eif_each_cal(med1, mus_IIE_M1, alphas_each_IIE_M1, m_each_IIE_M1)
    eifs_IIE_M1[["IIE_M1"]] <- eifs_IIE_M1[[2]] - eifs_IIE_M1[[1]]
    
    eifs_IIE_M2 <- eif_each_cal(med2, mus_IIE_M2, alphas_each_IIE_M2, m_each_IIE_M2)
    eifs_IIE_M2[["IIE_M2"]] <- eifs_IIE_M2[[2]] - eifs_IIE_M2[[1]]
    
    eifs_IIE_Mjo[["IIE_mu"]] <- eifs_IIE_Mjo[["IIE_Mjo"]] - eifs_IIE_M1[["IIE_M1"]] - eifs_IIE_M2[["IIE_M2"]]
    
    eifs <- c(eifs_IIE_M1, eifs_IIE_M2, eifs_IDE, eifs_IIE_Mjo) %>% do.call(cbind, .) %>% data.frame()
    estimates <- inference.iid(eifs)
    
    estimates <- estimates %>% rownames_to_column(var = "estimand") %>% 
      arrange(estimand)
    # estimates_I <- inference.cl(eifs_IIE_Mjo, data[["id"]], average = "individual", conf.level = 0.95)
    # estimates_C <- inference.cl(eifs_each, data_in[["S"]], average = "cluster")
    
    # sensitivity analysis and short parameters
    sens_fit_res <- sens_YM(med, mus, alphas_each, m_each, params, cf.y, cf.m, rho2)
    
    rv_res <- rv.value(med, mus, alphas_each, m_each, params, theta.null = 0, cf = seq(0, 0.99,by = 0.05))
    
    
    med_sens <- list(estimates_I, rv_res, sens_fit_res)
    names(med_sens) <- c("estimates_short", "robustness_value", "sens_fit_res")
    
    med_sens
}
