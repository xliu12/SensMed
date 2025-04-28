

mediate_s <- function(
        data,
        treatment, 
        outcome,
        mediators,
        covariates,
        a1 = 1, a0 = 0,
        id,
        cluster_id,
        weights = rep(1, nrow(data)),
        learners_outcome = c("glm", "xgboost"),
        learners_treatment = c("glm", "mean", "lightgbm"),
        cf.y = 0.13, cf.m = 0.13, rho2 = 1, conf.level = 0.95,
        control = estimate_control()
) {
    # data_in <- data
    # Sname <- cluster_id
    # AM <- data_in[[Aname]] * data_in[, Mnames]
    # colnames(AM) <- c("AM.M1", "AM.M2")
    # Mint <- data.frame(data_in[, Mnames[1]] * data_in[, Mnames[2]])
    # colnames(Mint) <- "Mint"
    # AMint <- data.frame(data_in[[Aname]] * data_in[, Mnames[1]] * data_in[, Mnames[2]])
    # colnames(AMint) <- "AMint"
    
    # dummy indicators
    # Sdumm <- dummy_cols(data[[cluster_id]], remove_first_dummy = FALSE, remove_selected_columns = TRUE)
    # colnames(Sdumm) <- paste0("S", 1:ncol(Sdumm))
    # Sname_dummies <- colnames(Sdumm)
    # 
    # cluster means
    data <- data %>% # data.frame(data, AM, Mint, AMint) %>%
        # #     # rename(A = all_of(treatment), S = all_of(cluster_id), Y = all_of(outcome), M = all_of(mediators)) %>%
        rename( S = all_of(cluster_id) ) %>%
        # #     # mutate(A0 = 1 - A) %>%
        group_by(S) %>%
        mutate(across(all_of(c(treatment, mediators, outcome, covariates)), 
                      list(clmean = ~mean(.), cwc = ~.-mean(.)))) # %>% bind_cols(Sdumm)
    
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
            S = NA_character_,
            id = id
        ),
        weights = weights,
        a0 = a0,
        a1 = a1
    )
    # Multiple mediators
    # Create independent M -----
    # this estimates mediation for the first mediator M1
    num_M <- length(med@vars@M)
    m_each_Mjo <- list(
        IIE_Mjo_m = c(az = "", am = "data_0", ay = "data_1"),
        IIE_Mjo_p = c(az = "", am = "data_1", ay = "data_1"),
        IDE_m = c(az = "", am = "data_0", ay = "data_0"),
        IDE_p = c(az = "", am = "data_0", ay = "data_1")
    ) 
    med1 <- med
    # med <- add_indepM.within(med, control, which.M = which.M, if_indepM = TRUE)
    m_each_IIE_M1 <- list(
        IIE_M1_m = c(az = "data_0_indepM", am = "data_0_indepM", ay = "data_1_indepM") 
        ,
        IIE_M1_p = c(az = "data_0_indepM", am = "data_1_indepM", ay = "data_1_indepM"))
    
    
    med2 <- med
    m_each_IIE_M2 <- list(
        IIE_M2_m = c(az = "data_1_indepM", am = "data_0_indepM", ay = "data_1_indepM") 
        ,
        IIE_M2_p = c(az = "data_1_indepM", am = "data_1_indepM", ay = "data_1_indepM"))
    
    
    folds <- make.fold_K(data, Snames = "S", control$crossfit_folds)
    #folds <- origami::make_folds(med@data, V = control$crossfit_folds)
    if (control$crossfit_folds == 1) {
        folds[[1]]$training_set <- folds[[1]]$validation_set
    }
    
    # Fit models ------------------------
    ## IIE_Mjo and IDE -----------------------
    
    #  mus <- estimate_mu_cl(med, m_each, folds, learners, control, which.M, covariates_cl)
    mus_Mjo <- estimate_mu_cl(med, m_each_Mjo, folds, learners_outcome, control, which.M = NULL, covariates_cl = F)
    # mus_Mjo3 <- estimate_mu_cl(med, m_each_Mjo, folds, learners_outcome, control, which.M = 2)
    
    # mean(mus_Mjo[[2]]$b_am_jo - mus_Mjo[[1]]$b_am_jo)
    # mean(mus_Mjo[[4]]$b_am_jo - mus_Mjo[[3]]$b_am_jo)
    # mean(mus_Mjo3[[2]]$b1 - mus_Mjo3[[1]]$b1)
    # mean(mus_Mjo3[[4]]$b1 - mus_Mjo3[[3]]$b1)
    
    # Estimate Riesz representer alphas
    # a_c <- a.c(med, varnames, cluster_opt = "noncluster.mlr", 
    #            # cluster_opt = "cwc", 
    #            folds, glue("SL.{learners}"), bounded = TRUE)
    #learners <- c("glm", "mean", "lightgbm")
    a_c <- a.c(med, folds, learners_treatment, control, covariates_cl = F)
    
    # learners <- c("glm", "xgboost")
    # a_mc <- a.mc(med, # varnames, cluster_opt = "noncluster.mlr", bounded = TRUE,
    #              # cluster_opt = "cwc", 
    #              folds, glue("SL.{learners}"), )
    a_mc <- a.mc(med, folds, learners_treatment, control, covariates_cl = F)
    
    
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
    mus_IIE_M1 <- estimate_mu_cl(med1, m_each_IIE_M1, folds, learners_outcome, control, which.M = 2, covariates_cl = F)
    
    # mean(mus_IIE_M1[[2]]$b1 - mus_IIE_M1[[1]]$b1)
    
    # nn_module <- sequential_module()
    # alphas_each <- estimate_alpha_cl(med, m_each, folds, nn_module, control, which.M, covariates_cl)
    # alphas_each <- estimate_alpha_each(med, m_each[c(5:6)], folds, nn_module, control, which.M)
    
    # aa=lrn("riesz.nn"); aa$param_set; mlr_learners
    control$learners_riesz <- list("nn" = list(.key = "nn", verbose = FALSE, interactions = 1, learning_rate = 0.001, hidden = 30, layers = 1))
    
    alphas_each_IIE_M1 <- riesz_crossfit(med1, m_each = m_each_IIE_M1, folds, which.M = 2, learners_riesz = control$learners_riesz, control, riesz_alpha1 = F, riesz_alpha3 = TRUE, alpha_ac, a_c, a_mc, covariates_cl = F)
    
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
    mus_IIE_M2 <- estimate_mu_cl(med2, m_each_IIE_M2, folds, learners_outcome, control, which.M = 1, covariates_cl = F)
    # mean(mus_IIE_M2[[2]]$b1 - mus_IIE_M2[[1]]$b1)
    
    control$learners_riesz <- list("nn" = list(.key = "nn", verbose = FALSE, interactions = 1, learning_rate = 0.001, hidden = 30, layers = 1))
    
    alphas_each_IIE_M2 <- riesz_crossfit(med2, m_each = m_each_IIE_M2, folds, which.M = 1, learners_riesz = control$learners_riesz, control, riesz_alpha1 = F, riesz_alpha3 = TRUE, alpha_ac, a_c, a_mc, covariates_cl = F)
    
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
    
    estimates <- inference.cl(eifs, data$S, average = "individual", conf.level = 0.95)
    # estimates_C <- inference.cl(eifs_each, data_in[["S"]], average = "cluster")
    estimates <- estimates %>% rownames_to_column(var = "estimand") %>% 
        arrange(estimand)
    
    estimates
    
}
