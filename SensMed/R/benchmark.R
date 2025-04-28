
benchmark.effect <- function(med, 
                           benchmark_covariates = list(a=covariates[c(11)]), #c("X_BYSES1", "X_BYPARED"),
                           param = "IIE_M2", mus, alphas_each, m_each, estimates, folds, control, learners_outcome, only_YM_confounded = TRUE, mediatedY = c(1)){
    
    
    # with the benchmark covariates ----------------
    ## nu2 for alpha -------------
    # alphas_each = alphas_each_IIE_M2_r
    # alphas_each = alphas_each_IIE_Mjo_r
    # alphas_each = alphas_each_IIE_M1_r[2]
    # cal_nu2 <- cal.nu2(alphas_each)
    # param = "IIE_Mjo"
    if (m_each[[1]][["az"]] == "") {
        if (param == "IIE_Mjo") {
            nu2_list <- cal.nu2.IIE_Mjo(alphas_each[ mediatedY ])
            nu2 <- nu2_list$nu2
            nu2plugin <- nu2_list$nu2plugin
        }
        if (param == "IDE") {
            nu2_list <- cal.nu2.IDE(alphas_each[ mediatedY ])
            nu2 <- nu2_list$nu2
            nu2plugin <- nu2_list$nu2plugin
        }
        which.M <- NULL
    } else {
        # nu2 <- map(1:2, \(i) cal.nu2.IIE_M1(alphas_each[i])$nu2 %>% unlist() )
        nu2 <- cal.nu2.IIE_M1(alphas_each[ mediatedY ])$nu2
        which.M <- 2
    }
    
    # colMeans((alphas_each[[2]] - alphas_each[[1]])^2)

    ## R2.Y -----------------
    # mus <- mus_IIE_M2; m_each <- m_each_IIE_M2
    # mus <- mus_IIE_M1; m_each <- m_each_IIE_M1
    # mus <- mus_Mjo[c(1:2)]; m_each <- m_each_Mjo[1:2]
    # med <- med2_cwc
    Y <- med@data[[med@vars@Y]]
    # resY <- cal.resY(Y, m_each, mus) # map(1:2, \(i) cal.resY(Y, m_each[i], mus[i]  ))
    # R2_Y <- cal.R2Y(Y, m_each, mus) # map(1:2, \(i) cal.R2Y(Y, m_each[i], mus[i]) %>% unlist())
    
    
    # omit the benckmark covariates ------------
    benchmarks <- list()
    i <- 1
    cli::cli_progress_message("Analyzing benchmarks using covariate: {benchmark_covariates[[i]]}")
    
    # med <- med1_cwc
    for (i in seq_along(benchmark_covariates)) {
        if (length(grep("_cwc$", med@vars@C)) > 0) {
            covar <- glue("{benchmark_covariates[[i]]}_cwc")
        } else {
            covar <- benchmark_covariates[[i]]
        }
        
        # covar <- "X_BYPARED"
        # med <- med1_cwc
        # if (YM_confound_only) {
        #     
        # }
        # if (!YM_confound_only) {
        #     
        # }
        # med_oc <- med2_cwc
        # m_each <- m_each_IIE_M2
        
        # only_YM_confounded = TRUE
        if (only_YM_confounded) {
            med_oc <- med
            med_oc@vars@C <- setdiff(med_oc@vars@C, covar)
            mus_long_short <- mus
        } 
        if (!only_YM_confounded) {
            
            med_oc <- med 
            med_oc@vars@C <- setdiff(med_oc@vars@C, covar)
            
            ## Permute with omitted benchmark covariate
            cluster_id <- med@vars@cluster_id
            if (!m_each[[1]][["az"]] == "") {
                set.seed(1234)
                if (is.na(cluster_id)) {
                    # create M_2^\pi | A, X^\short_-j
                    med_oc <- add_indepM(med_oc, control, which.M = which.M, cwc = FALSE, Sname_dummies = NULL)
                }
                if (!is.na(cluster_id)) {
                    med_oc <- add_indepM.within(med_oc, control, which.M = which.M, if_indepM = FALSE)
                }
                
                # theta_{2, theta_3^\ss}
                # med_tmp <- med_oc
                
            }
            
            med_tmp <- med_oc
            med_tmp@vars@C <- med@vars@C
            
            mus_long_short <- estimate_mu.long.short(covar, med_tmp, m_each, folds, learners_outcome, control, which.M = which.M, covariates_cl = FALSE)
        }
        
        
        ## R2_Y_oc ----
        set.seed(1234)
        mus_oc <- estimate_mu_cl(med_oc, m_each, folds, learners_outcome, control, which.M = which.M, covariates_cl = FALSE, benchmark_covar = covar)
        
        resY_oc <- cal.resY(Y, m_each, mus_oc[ mediatedY ]) # map(1:2, \(i) cal.resY(Y, m_each[i], mus_oc[i]  ))
        R2_Y_oc <- cal.R2Y(Y, m_each, mus_oc[ mediatedY ]) # map(1:2, \(i) cal.R2Y(Y, m_each[i], mus_oc[i]) %>% unlist())
        
        # nu2_oc ---------------
        nn_module <- sequential_module(dropout = 0.05)
        # control$device="cpu"
        set.seed(1234)
        alphas_oc <- estimate_alpha_each(med_oc, m_each, folds, nn_module, control, which.M = which.M, benchmark_covar = covar)
        
        if (m_each[[1]][["az"]] == "") {
            # mediatedY <- grep("am0_ay1", names(alphas_each))
            if (param == "IIE_Mjo") {
                nu2_oc_list <- cal.nu2.IIE_Mjo(alphas_oc[ mediatedY ])
                nu2_oc <- nu2_oc_list$nu2
                nu2plugin_oc <- nu2_oc_list$nu2plugin
            }
            if (param == "IDE") {
                nu2_oc_list <- cal.nu2.IDE(alphas_oc[ mediatedY ]) 
                nu2_oc <- nu2_oc_list$nu2
                nu2plugin_oc <- nu2_oc_list$nu2plugin
            }
            # nu2_oc <- map(1:2, \(i) cal.nu2.IIE_Mjo(alphas_oc[i])$nu2 %>% unlist() )
            which.M <- NULL
        } else {
            nu2_oc <- cal.nu2.IIE_M1(alphas_oc[ mediatedY ])$nu2
            # nu2_oc <- map(1:2, \(i) cal.nu2.IIE_M1(alphas_oc[i])$nu2 %>% unlist() )
            which.M <- 2
        }

       # mean(( alphas_oc[[2]]$alpha1 - alphas_oc[[1]]$alpha1 )* data[[covar]])
        
        eifs_oc <- eif_each_cal(med_oc, mus_oc, alphas_oc, m_each)
        
        # param <- "IIE_M1"
        eifs_oc[[ param ]] <- eifs_oc[[2]] - eifs_oc[[1]]
        eifs_oc <- eifs_oc %>% do.call(cbind, .) %>% data.frame()
        
        # cluster_id <- med_oc@vars@cluster_id
        if (is.na(cluster_id)) {
            S_cluster <- NA
        } else {
            S_cluster <- med@data[["S"]]
        }
        estimates_oc <- inference.cl(eifs_oc, S = S_cluster, average = "individual", conf.level = 0.95)
        
        
        # estimates_C <- inference.cl(eifs_each, data_in[["S"]], average = "cluster")
        estimates_oc <- estimates_oc %>% 
            rownames_to_column(var = "estimand") %>%
            arrange(estimand)
        
        ## C_alpha 
        
        # alpha_l_list_diff <- alphas_each[[2]] - alphas_each[[1]]
        # alpha_s_list_diff <- alphas_oc[[2]] - alphas_oc[[1]]
        # C2_alpha <- map(1:length(nu2_oc), \(x) mean((alpha_l_list_diff[[x]] - alpha_s_list_diff[[x]])^2) / mean(alpha_s_list_diff[[x]]^2) )
        # C2_alpha <- map(1:length(nu2_oc), \(x) (mean(alpha_l_list_diff[[x]]^2) - mean(alpha_s_list_diff[[x]]^2)) / mean(alpha_s_list_diff[[x]]^2) )
        # use the non-debiased, plug-in. follow Hong et al 2018. not do effect. do potential outcomes
        # Csq_alpha_m <- map(1:length(nu2), \(x) (mean(alphas_each[[1]][[x]]^2) - mean(alphas_oc[[1]][[x]]^2)) / mean(alphas_oc[[1]][[x]]^2) )
        # Csq_alpha_p <- map(1:length(nu2), \(x) (mean(alphas_each[[2]][[x]]^2) - mean(alphas_oc[[2]][[x]]^2)) / mean(alphas_oc[[2]][[x]]^2) )  %>% unlist()
        
        
        # Csq_alpha <- map(1:length(nu2), \(x) (nu2[[x]] - nu2_oc[[x]]) / nu2_oc[[x]] ) %>% unlist()
        # # cf_m is 1 - R^2_{alpha ~ alpha.s}
        # cf_m_bench <- Csq_alpha / (1 + Csq_alpha)
        
        # only Y-M confounding
        # Csq_alpha <- (nu2[["nu3"]] - nu2_oc[["nu3"]]) / nu2_oc[["nu3"]]
        # cf_m_bench <- Csq_alpha / (1 + Csq_alpha)
        
        # names(cf_m_bench) <- names(nu2)
        
        # map(1:3, \(x) (mean((alphas_each[[2]][[x]] - alphas_oc[[2]][[x]])^2)) / mean(alphas_oc[[2]][[x]]^2) )
        
        
        ## C_theta 
        
        # Ypseudo_list <- data.frame(Ypseudo4 = Y, 
        #                            Ypseudo3 = mus_long_short[[1]][["b3"]], 
        #                            Ypseudo2 = mus_long_short[[1]][["b2"]])
        
        # mediatedY <- grep("am0_ay1", names(mus_long_short))
        # Ypseudo_list <- data.frame(Ypseudo3 = Y, 
        #                            Ypseudo2 = mus_long_short[[mediatedY]][["b_ay_jo"]])
        # # cf_y
        # mus_mediatedY <- mus_long_short[[mediatedY]] %>% 
        #     rename("mu_fit2" = mu_fit_ay_jo,
        #            "mu_fit1" = mu_fit_am_jo)
        # 
        # mediatedY <- grep("am0_ay1", names(mus_oc))
        # mus_oc_mediatedY <- mus_oc[[mediatedY]] %>% 
        #     rename("mu_fit2" = mu_fit_ay_jo,
        #            "mu_fit1" = mu_fit_am_jo)
        # 
        # Csq_theta <- map(1:2, \(x=1) {
        #     resid_x_num <- mus_mediatedY[[glue("mu_fit{x}")]] - mus_oc_mediatedY[[glue("mu_fit{x}")]]
        #     resid_x_den <- Ypseudo_list[[glue("Ypseudo{x+1}")]] - mus_oc_mediatedY[[glue("mu_fit{x}")]]
        #     mean(resid_x_num^2) / mean(resid_x_den^2)
        # })
        # cf_y_bench <- Csq_theta %>% unlist()
        # names(cf_y_bench) <- glue("theta{1:2}")
        
        # only Y-M confounding
        # mus[[2]]$mu_fit3 - mus[[2]]$mu_fit3
        
        # rho  
        # not debiased
        # alphas_mediatedY <- alphas_each[[ grep("am0_ay1", names(alphas_each)) ]] %>% 
        #     select(alpha1, alpha3) %>% 
        #     rename("alpha2" = alpha3)
        # alphas_oc_mediatedY <- alphas_oc[[ grep("am0_ay1", names(alphas_oc)) ]] %>% 
        #     select(alpha1, alpha3) %>% 
        #     rename("alpha2" = alpha3)
        # 
        # rho <- map(1:2, \(x) {
        #     resid_theta <- mus_mediatedY[[glue("mu_fit{x}")]] - mus_oc_mediatedY[[glue("mu_fit{x}")]]
        #     resid_alpha <- alphas_mediatedY[[x]] - alphas_oc_mediatedY[[x]]
        #     cor(resid_theta, resid_alpha)
        # })
        # rho <- unlist(rho)
        # names(rho) <- glue("rho{1:2}")
        
        # Bias ----
        estimates_ocwc <- left_join(estimates_oc, estimates, by = "estimand", suffix = c("_oc", "_wc"))
        Bias <- estimates_ocwc %>% dplyr::filter(str_detect(estimand, "I")) %>%
        mutate(Bias = estimate_oc - estimate_wc) %>% pull(Bias)
        
        
        # from dml.sensmaker
        
        # Benchmark gain metrics -----
        resY <- cal.resY(Y, m_each, mus_long_short[ mediatedY ])
        
        V_g <- apply(resY_oc, 2, var) - apply(resY, 2, var)
        # map(1:length(resY), \(y) apply(resY_oc[[y]], 2, var) - apply(resY[[y]], 2, var))
        V_a <- map(1:length(nu2), \(y) 
                   {dd <- nu2[[y]] - nu2_oc[[y]]
                       if(dd < 0) { dd <- nu2plugin[[y]] - nu2plugin_oc[[y]] }
                       dd
                       }
                   )

        Cor <- map(1:length(V_g), \(y) {
            valid <- V_g[[y]] > 0 & V_a[[y]] > 0
            Cor <- 0
            Cor[valid] <- (abs(Bias[valid])/sqrt(V_g[[y]][valid]*V_a[[y]][valid]))*sign(Bias[valid])
        }) #%>% unlist()
        names(Cor) <- glue("rho{1:length(Cor)}")
        
        R2_Y <- cal.R2Y(Y, m_each, mus_long_short[ mediatedY ])

        Gain_mu <- map(1:length(R2_Y), \(x) pmax(0, (R2_Y[[x]] - R2_Y_oc[[x]]) / (1 - R2_Y[[x]]))) %>% unlist()
        names(Gain_mu) <- names(R2_Y)
        # Gain_mu3 <- Gain_mu$mu3
        
        
        Gain_alpha <- map(1:length(nu2), \(x) {
            dd <- (nu2[[x]] - nu2_oc[[x]]) / nu2_oc[[x]]
            if (dd < 0) { 
                dd <- (nu2plugin[[x]] - nu2plugin_oc[[x]]) / nu2plugin_oc[[x]]
            }
            dd
        }) %>% unlist()
        names(Gain_alpha) <- names(nu2)
        # Gain_alpha3 <- Gain_alpha$nu3
        
        # if (Gain_alpha3 == 0) {
        #     nusq_3 <- mean((alphas_each[[2]]$alpha3 - alphas_each[[1]]$alpha3)^2)
        #     nusq_oc_3 <- mean((alphas_oc[[2]]$alpha3 - alphas_oc[[1]]$alpha3)^2)
        #     Gain_alpha3 <- (nusq_3 - nusq_oc_3) / nusq_oc_3
        # }
        benchmarks[[paste(covar, collapse = ".")]] <- list(Gain_mu = Gain_mu,
                      Gain_alpha = Gain_alpha,
                      rho = Cor,
                      estimates_ocwc  = estimates_ocwc ,
                      Bias = Bias)
        
        # benchmarks[[paste(covar, collapse = ".")]] <- list(C2_alpha = Gain_alpha3,
        #                             C2_theta = Gain_mu3,
        #                             # rho = rho,
        #                             estimates_ocwc  = estimates_ocwc ,
        #                             Bias = Bias)
        
    }
    
    benchmarks
}
