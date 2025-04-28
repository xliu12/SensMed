#' @export
#'
Sens.Mediate <- function(
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
    benchmark_covariates = NULL,
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
          M = unlist(mediators),
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
    
    # Create independent M -----
    # this estimates mediation for the first mediator M1
    num_M <- length(mediators)

    if ((control$indepM == FALSE)) {
        m_each_Mjo <- list(
            IIE_Mjo_m = c(az = "", am = "data_0", ay = "data_1"),
            IIE_Mjo_p = c(az = "", am = "data_1", ay = "data_1"),
            IDE_m = c(az = "", am = "data_0", ay = "data_0"),
            IDE_p = c(az = "", am = "data_0", ay = "data_1")
        )

        # IIE_M1, independently draw M2 from its marginal distribution
        set.seed(1234)
        if (is.na(cluster_id)) {
          med1 <- add_indepM(med, control, which.M = 2, cwc = FALSE, Sname_dummies = NULL)
        }
        if (!is.na(cluster_id)) {
          med1 <- add_indepM.within(med, control, which.M = 2, if_indepM = FALSE)
        }
        
        # cor(med1@data_indepM[mediators])
        # cor(med1@data[mediators])

        m_each_IIE_M1 <- list(
            IIE_M1_m = c(az = "data_0_indepM", am = "data_0_indepM", ay = "data_1_indepM"),
            IIE_M1_p = c(az = "data_0_indepM", am = "data_1_indepM", ay = "data_1_indepM")
        )
        # m_each_IIE_M1 <- list(
        #     IIE_M1_pm = c(az = "data_0_indepM", am.p = "data_1_indepM", am.m = "data_0_indepM", ay = "data_1_indepM")
        # )

        # IIE_M2, independently draw M1 from its marginal distribution
        med2 <- med
        med2@vars@M <- rev(med@vars@M)
        set.seed(1234)
        if (is.na(cluster_id)) {
          med2 <- add_indepM(med2, control, which.M = c(1:num_M)[-1], cwc = FALSE, Sname_dummies = NULL)
        }
        if (!is.na(cluster_id)) {
          med2 <- add_indepM.within(med2, control, which.M = 2, if_indepM = FALSE)
        }

        m_each_IIE_M2 <- list(
            IIE_M2_m = c(az = "data_1_indepM", am = "data_0_indepM", ay = "data_1_indepM"),
            IIE_M2_p = c(az = "data_1_indepM", am = "data_1_indepM", ay = "data_1_indepM")
        )
        # med3 <- add_indepM(med, control, which.M = c(1:num_M)[-3], cwc = FALSE)
        # medmu <- add_indepM(med, control, which.M = c(1:num_M), cwc = FALSE)
        # m_each_IIE_mu <- list(
        #   IIE_M1_m = c(az = "data_0_indepM", am = "data_0_indepM", ay = "data_1_indepM")
        #   ,
        #   IIE_M1_p = c(az = "data_1_indepM", am = "data_1_indepM", ay = "data_1_indepM"))
    }

    # Assume conditionally independent M -----------
    # run single-mediator analysis for each
    if (num_M >= 2 & (control$indepM == TRUE)) {
        m_each_Mjo <- list(
            IIE_Mjo_m = c(az = "", am = "data_0", ay = "data_1"),
            IIE_Mjo_p = c(az = "", am = "data_1", ay = "data_1"),
            IDE_m = c(az = "", am = "data_0", ay = "data_0"),
            IDE_p = c(az = "", am = "data_0", ay = "data_1")
        )

        # IIE_M1
        med1 <- med
        med1@vars@M <- med@vars@M[1] 
        # med <- add_indepM.within(med, control, which.M = which.M, if_indepM = TRUE)
        m_each_IIE_M1 <- list(
            IIE_M1_m = c(az = "", am = "data_0", ay = "data_1"),
            IIE_M1_p = c(az = "", am = "data_1", ay = "data_1")
        )
        # m_each_IIE_M1 <- list(
        #   IIE_M1_pm = c(az = "data_0", am.p = "data_1", am.m = "data_0", ay = "data_1")
        # )

        # IIE_M2, independently draw M1 from its marginal distribution
        med2 <- med
        med2@vars@M <- med@vars@M[2]
        # med2 <- add_indepM(med, control, which.M = c(1:num_M)[-2], cwc = FALSE, Sname_dummies = NULL)

        m_each_IIE_M2 <- list(
            IIE_M2_m = c(az = "", am = "data_0", ay = "data_1"),
            IIE_M2_p = c(az = "", am = "data_1", ay = "data_1")
        )
        # m_each_IIE_M2 <- list(
        #     IIE_M1_pm = c(az = "data_1", am.p = "data_1", am.m = "data_0", ay = "data_1")
        # )
    }
    
    # only one mediator
    if (num_M == 1) {
      m_each_Mjo <- list(
        IIE_Mjo_m = c(az = "", am = "data_0", ay = "data_1"),
        IIE_Mjo_p = c(az = "", am = "data_1", ay = "data_1"),
        IDE_m = c(az = "", am = "data_0", ay = "data_0"),
        IDE_p = c(az = "", am = "data_0", ay = "data_1")
      )
      
    }

    # folds <- make.fold_K(data, Snames = "S", control$crossfit_folds)
    set.seed(1234)
    if (is.na(cluster_id)) {
      folds <- origami::make_folds(med@data, V = control$crossfit_folds)
    }
    if (!is.na(cluster_id)) {
      folds <- origami::make_folds(med@data, cluster_ids = med@data$S, V = control$crossfit_folds)
      
    }
    if (control$crossfit_folds == 1) {
        folds[[1]]$training_set <- folds[[1]]$validation_set
    }

    # devtools::load_all("/Users/xl9663/Library/CloudStorage/Box-Box/Labs/Mediated/R_M2more/crumble-main")
    # ?crumble::crumble
    # dat_crumble <- data %>% ungroup() %>% mutate(id = 1:nrow(data))
    # crumble_res <- crumble::crumble(
    #     data = dat_crumble,
    #     trt = treatment,
    #     outcome = outcome, mediators = mediators[1], moc = NULL,
    #     covar = covariates,
    #     d0 = \(data, trt) rep(0, nrow(data)),
    #     d1 = \(data, trt) rep(1, nrow(data)),
    #     effect = c("N"),
    #     learners = learners_outcome,
    #     nn_module = crumble::sequential_module(dropout = 0.05),
    #     id = "id",
    #     control = crumble::crumble_control(crossfit_folds = 5, mlr3superlearner_folds = 5, device = "cpu", epochs = 10)
    # )

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

    # IIE_Mjo and IDE -----------------------
    # Estimate mu
    if (is.null(learners_outcome)) {
        learners_outcome <- c("glm", "nnet")
    }
    
    med_cwc <- med
    if (!is.na(cluster_id)) {
      med_cwc@vars@A <- glue("{med@vars@A}_cwc")
      med_cwc@vars@M <- glue("{med@vars@M}_cwc")
      med_cwc@vars@Y <- glue("{med@vars@Y}_cwc")
      med_cwc@vars@C <- glue("{med@vars@C}_cwc")
    }
    
    set.seed(1234)
    # mus_Mjo <- estimate_mu_cl(med, m_each_Mjo, folds, learners_outcome, control, which.M = NULL, covariates_cl = FALSE)
    
    mus_Mjo <- estimate_mu_cl(med_cwc, m_each_Mjo, folds, learners_outcome, control, which.M = NULL, covariates_cl = FALSE)

    # mean(mus_Mjo[[2]]$b_am_jo - mus_Mjo[[1]]$b_am_jo)
    # mean(mus_Mjo[[4]]$b_am_jo - mus_Mjo[[3]]$b_am_jo)
    # mean(mus_Mjo3[[2]]$b1 - mus_Mjo3[[1]]$b1)
    # mean(mus_Mjo3[[4]]$b1 - mus_Mjo3[[3]]$b1)

    # Estimate Riesz representer alphas
    nn_module <- sequential_module(dropout = 0.05)
    # control$epochs <- 12
    set.seed(1234)
    # alphas_each_IIE_Mjo_r <- estimate_alpha_each(med, m_each_Mjo[c(1:2)], folds, nn_module, control, which.M = NULL)
    
    alphas_each_IIE_Mjo_r <- estimate_alpha_each(med_cwc, m_each_Mjo[c(1:2)], folds, nn_module, control, which.M = NULL)
    set.seed(1234)
    alphas_each_IDE_r <- estimate_alpha_each(med_cwc, m_each_Mjo[c(3:4)], folds, nn_module, control, which.M = NULL)

    # colMeans((alphas_each_IIE_Mjo_r[[2]] - alphas_each_IIE_Mjo_r[[1]]) * med_cwc@data[[med_cwc@vars@C[1]]])
    # learners_treatment <- c("glm", "nnet", "mean")
    # a_c <- a.c(med, folds, learners_treatment, control, covariates_cl = FALSE)
    # a_mc <- a.mc(med, M = med@vars@M, folds, learners_treatment, control, covariates_cl = FALSE)
    # alpha_ac <- alpha.ac(A = med@data[[med@vars@A]], a_c)
    # mean((alpha_ac$`A=0`) * data$Y)

    # alphas_each_IIE_Mjo <- riesz_crossfit(med, m_each = m_each_Mjo[c(1:2)], folds, which.M = NULL, learners_riesz = control$learners_riesz, control, riesz_alpha1 = FALSE, riesz_alpha2 = FALSE, riesz_alpha3 = FALSE, alpha_ac, a_c, a_mc, a_m1c, a_m2c, r_zmac, covariates_cl = FALSE)
    # alphas_each_IDE <- riesz_crossfit(med, m_each = m_each_Mjo[c(3:4)], folds, which.M = NULL, learners_riesz = control$learners_riesz, control, riesz_alpha1 = FALSE, riesz_alpha2 = FALSE, riesz_alpha3 = FALSE, alpha_ac, a_c, a_mc, a_m1c, a_m2c, r_zmac, covariates_cl = FALSE)

    # mean(alphas_each_IDE[[2]]$alpha3 * data[[med@vars@Y]]) - mean(alphas_each_IDE[[1]]$alpha3 * data[[med@vars@Y]])
    # mean(alphas_each_IIE_Mjo[[2]]$alpha3 * data[[med@vars@Y]]) - mean(alphas_each_IIE_Mjo[[1]]$alpha3 * data[[med@vars@Y]])

    # control$learners_riesz <- list("nn" = list(.key = "nn", verbose = FALSE, interactions = 1, learning_rate = 0.001, hidden = 20, layers = 1))

    # alphas_each_IIE_Mjo_r <- riesz_crossfit(med, m_each = m_each_Mjo[c(1:2)], folds, which.M = NULL, learners_riesz = control$learners_riesz, control, riesz_alpha1 = TRUE, riesz_alpha2 = TRUE, riesz_alpha3 = TRUE, alpha_ac, a_c, a_mc, a_m1c, a_m2c, r_zmac, covariates_cl = FALSE)

    # alphas_each_IDE_r <- riesz_crossfit(med, m_each = m_each_Mjo[c(3:4)], folds, which.M = NULL, learners_riesz = control$learners_riesz, control, riesz_alpha1 = TRUE, riesz_alpha2 = TRUE, riesz_alpha3 = TRUE, alpha_ac, a_c, a_mc, a_m1c, a_m2c, r_zmac, covariates_cl = FALSE)

    # mean(alphas_each_IDE_r[[2]]$alpha3 * data[[med@vars@Y]]) - mean(alphas_each_IDE_r[[1]]$alpha3 * data[[med@vars@Y]])
    # mean(alphas_each_IIE_Mjo_r[[2]]$alpha3 * data[[med@vars@Y]]) - mean(alphas_each_IIE_Mjo_r[[1]]$alpha3 * data[[med@vars@Y]])
    
    
    ## EIF -----
    eifs_IDE <- eif_each_cal(med_cwc, mus_Mjo[c(3:4)], alphas_each_IDE_r, m_each_Mjo[3:4])
    eifs_IDE[["IDE"]] <- eifs_IDE[[2]] - eifs_IDE[[1]]
    
    # eifs_IIE_Mjo <- eif_each_cal(med, mus_Mjo[c(1:2)], alphas_each_IIE_Mjo_r, m_each_Mjo[1:2])
    eifs_IIE_Mjo <- eif_each_cal(med_cwc, mus_Mjo[c(1:2)], alphas_each_IIE_Mjo_r, m_each_Mjo[1:2])
    eifs_IIE_Mjo[["IIE_Mjo"]] <- eifs_IIE_Mjo[[2]] - eifs_IIE_Mjo[[1]]
    

    if (num_M > 1) {
      # IIE_M1 --------------
      med1_cwc <- med1
      if (!is.na(cluster_id)) {
        med1_cwc@vars@A <- glue("{med1@vars@A}_cwc")
        med1_cwc@vars@M <- glue("{med1@vars@M}_cwc")
        med1_cwc@vars@Y <- glue("{med1@vars@Y}_cwc")
        med1_cwc@vars@C <- glue("{med1@vars@C}_cwc")
      }
      
      # med1_cwcSdumm <- med1_cwc
      # med1_cwcSdumm@vars@C <- c(glue("{covariates}_cwc"), Sname_dummies)
      if(length(med1_cwc@vars@M) == 1) {
        which.M <- NULL
      } else {
        which.M <- 2
      }
      # learners_outcome <- c("glm", "nnet")
      set.seed(1234)
      mus_IIE_M1 <- estimate_mu_cl(med1_cwc, m_each_IIE_M1, folds, learners_outcome, control, which.M = which.M, covariates_cl = FALSE)
      
      # mean(mus_IIE_M1$az0_am.p1_am.m0_ay1$b1)
      # mean(mus_IIE_M1[[2]]$b1)
      # mean(mus_IIE_M1[[1]]$b1)
      # mean(mus_IIE_M1[[2]]$b1 - mus_IIE_M1[[1]]$b1)
      # a_m1c <- a.mc(med, M = med@vars@M[1], folds, learners_treatment, control, covariates_cl = F)
      # a_m2c <- a.mc(med, M = med@vars@M[2], folds, learners_treatment, control, covariates_cl = F)
      
      # for IIE_M1
      # r_zm1ac <- r.zmac(med = med1, permuted = med@vars@M[1], folds, learners_treatment, control, covariates_cl = F)
      
      # control$epochs <- 12
      # control$learning_rate <- 0.01
      nn_module <- sequential_module(dropout = 0.05)
      
      set.seed(1234)
      alphas_each_IIE_M1_r <- estimate_alpha_each(med1_cwc, m_each_IIE_M1, folds, nn_module, control, which.M = which.M)
      # alphas_each_IIE_M1_r <- estimate_alpha_each(med1, m_each_IIE_M1, folds, nn_module, control, which.M = which.M)
      
      # mean(alphas_each[[2]]$alpha3 * data[[med@vars@Y]]) - mean( alphas_each[[1]]$alpha3 * data[[med@vars@Y]])
      # mean(alphas_each_IIE_M1_r[[2]]$alpha1 * data[[covariates[1]]] ) - mean( alphas_each_IIE_M1_r[[1]]$alpha1 * data[[covariates[1]]] )
      
      # library(SuperRiesz); aa=lrn("riesz.nn"); aa$param_set; mlr_learners
      # control$learners_riesz <- list("nn" = list(.key = "nn", verbose = FALSE, learning_rate = 0.01, hidden = 20, layers = 1, epochs = 10L, dropout = 0.1))
      
      # alphas_each_IIE_M1 <- riesz_crossfit(med1, m_each = m_each_IIE_M1, folds, which.M = 2, learners_riesz = control$learners_riesz, control, riesz_alpha1 = FALSE, riesz_alpha2 = FALSE, riesz_alpha3 = FALSE, alpha_ac, a_c, a_mc, a_m1c = a_m1c, a_m2c = a_m2c, r_zm1ac, covariates_cl = F)
      # m_each_IIE_M1 <- list(
      #   IIE_M1_pm = c(az = "data_0", am.p = "data_1", am.m = "data_0", ay = "data_1")
      # )
      
      # alphas_each_IIE_M1_r1 <- riesz_crossfit(med1, m_each = m_each_IIE_M1, folds, which.M = 2, learners_riesz = control$learners_riesz, control, riesz_alpha1 = TRUE, riesz_alpha2 = TRUE, riesz_alpha3 = TRUE, alpha_ac, a_c, a_mc, a_m1c = a_m1c, a_m2c = a_m2c, r_zm1ac, covariates_cl = FALSE)
      
      # mean(alphas_each_IIE_M1_r$az0_am.p1_am.m0_ay1$alpha1 * mus_IIE_M1$az0_am.p1_am.m0_ay1$mu_fit1)
      # mean(alphas_each_IIE_M1_r1$az0_am.p1_am.m0_ay1$alpha1 * data$X1)
      
      # mean(alphas_each_IIE_M1_r[[1]]$alpha3 * data[[med@vars@Y]])
      # mean(alphas_each_IIE_M1_r[[2]]$alpha3 * data[[med@vars@Y]]) - mean( alphas_each_IIE_M1_r[[1]]$alpha3 * data[[med@vars@Y]])
      
      
      # IIE_M2 (switch) ---------------
      med2@vars@M
      med2_cwc <- med2
      if (!is.na(cluster_id)) {
        med2_cwc@vars@A <- glue("{med2@vars@A}_cwc")
        med2_cwc@vars@M <- glue("{med2@vars@M}_cwc")
        med2_cwc@vars@Y <- glue("{med2@vars@Y}_cwc")
        med2_cwc@vars@C <- glue("{med2@vars@C}_cwc")
      }
      
      if(length(med2_cwc@vars@M) == 1) {
        which.M <- NULL
      } else {
        which.M <- 2
      }
      set.seed(1234)
      mus_IIE_M2 <- estimate_mu_cl(med2_cwc, m_each_IIE_M2, folds, learners_outcome, control, which.M = which.M, covariates_cl = FALSE)
      
      # mean(mus_IIE_M2[[1]]$b1)
      # mean(mus_IIE_M2[[2]]$b1)
      
      # nn_module <- sequential_module(dropout = 0.05)
      
      # control$epochs = 10
      set.seed(1234)
      alphas_each_IIE_M2_r <- estimate_alpha_each(med2_cwc, m_each_IIE_M2, folds, nn_module, control, which.M = which.M)
      
      # r_zm2ac <- r.zmac(med = med2, permuted = med@vars@M[1], folds, learners_treatment, control, covariates_cl = F)
      # alphas_each_IIE_M2 <- riesz_crossfit(med2, m_each = m_each_IIE_M2, folds, which.M = 2, learners_riesz = control$learners_riesz, control, riesz_alpha1 = FALSE, riesz_alpha2 = FALSE, riesz_alpha3 = FALSE, alpha_ac, a_c, a_mc, a_m1c = a_m2c, a_m2c = a_m1c, r_zm2ac, covariates_cl = F)
      
      # control$learners_riesz <- list("nn" = list(.key = "nn", verbose = FALSE, interactions = 1, learning_rate = 0.001, hidden = 50, layers = 1))
      
      # alphas_each_IIE_M2_r <- riesz_crossfit(med2, m_each = m_each_IIE_M2, folds, which.M = 2, learners_riesz = control$learners_riesz, control, riesz_alpha1 = TRUE, riesz_alpha2 = TRUE, riesz_alpha3 = TRUE, alpha_ac, a_c, a_mc, a_m1c = a_m2c, a_m2c = a_m1c, r_zm2ac, covariates_cl = FALSE)
      # mean(alphas_each_IIE_M2_r[[2]]$alpha3 * data$Y) - mean( alphas_each_IIE_M2_r[[1]]$alpha3 * data$Y)
      
      # IIE_M2
      # med2_cwc <- med2
      # med2_cwc@vars <- med_vars(A = "A_cwc",
      #                           Y = "Y_cwc",
      #                           M = glue("M{1:length(mediators)}_cwc"),
      #                           C = glue("{covariates}_cwc"),
      #                           cluster_id = "S", id = id)
      # med2_cwcSdumm <- med2_cwc
      # med2_cwcSdumm@vars@C <- c(glue("{covariates}_cwc"), Sname_dummies)
      
      #  mus <- estimate_mu_cl(med, m_each, folds, learners, control, which.M, covariates_cl)
      
      # mus_IIE_M2 <- estimate_mu_cl(med2, m_each_IIE_M2, folds, learners_outcome, control, which.M = 1, covariates_cl = F)
      # mean(mus_IIE_M2[[2]]$b1 - mus_IIE_M2[[1]]$b1)
      
      # control$learners_riesz <- list("nn" = list(.key = "nn", verbose = FALSE, interactions = 1, learning_rate = 0.001, hidden = 50, layers = 1))
      # r_zm2ac <- r.zmac(med = med2, permuted = med@vars@M[2], folds, learners_treatment, control)
      # r_zm2ac <- r_zm1ac
      
      # alphas_each_IIE_M2 <- riesz_crossfit(med2, m_each = m_each_IIE_M2, folds, which.M = 1, learners_riesz = control$learners_riesz, control, riesz_alpha1 = FALSE, riesz_alpha2 = FALSE, riesz_alpha3 = FALSE, alpha_ac, a_c, a_mc, a_m1c, a_m2c, r_zm2ac, covariates_cl = FALSE)
      
      # mean(alphas_each_IIE_M2[[2]]$alpha3 * data$Y) - mean( alphas_each_IIE_M2[[1]]$alpha3 * data$Y)
      
      ## EIF -----
      
      eifs_IIE_M1 <- eif_each_cal(med1_cwc, mus_IIE_M1, alphas_each_IIE_M1_r, m_each_IIE_M1)
      eifs_IIE_M1[["IIE_M1"]] <- eifs_IIE_M1[[2]] - eifs_IIE_M1[[1]]
      
      eifs_IIE_M2 <- eif_each_cal(med2_cwc, mus_IIE_M2, alphas_each_IIE_M2_r, m_each_IIE_M2)
      eifs_IIE_M2[["IIE_M2"]] <- eifs_IIE_M2[[2]] - eifs_IIE_M2[[1]]
      
      eifs_IIE_Mjo[["IIE_mu"]] <- eifs_IIE_Mjo[["IIE_Mjo"]] - eifs_IIE_M1[["IIE_M1"]] - eifs_IIE_M2[["IIE_M2"]]
      
      eifs <- c(eifs_IIE_M1, eifs_IIE_M2) %>% do.call(cbind, .) %>% data.frame()
      
      eifs <- c(eifs_IIE_M1, eifs_IIE_M2, eifs_IDE, eifs_IIE_Mjo) %>%
        do.call(cbind, .) %>%
        data.frame()
    }
    
    eifs <- c(eifs_IDE, eifs_IIE_Mjo) %>%
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

    # crumble_res$outcome_reg$theta_n$natural$fit3_natural[1,]
    
    ## RV values ----
    sens_fit_res <- list()
    rv_res <- list()
    
    if (num_M > 1) {
      # IIE_M1
      sens_fit_res[["IIE_M1"]] <- sens_IIE_Mjo(med1_cwc, mus_IIE_M1, alphas_each_IIE_M1_r, m_each_IIE_M1, params = "IIE_M1", cf.y, cf.m, rho2, only_YM_confounded = FALSE, mediatedY = 1:2)
      
      # rv_res[["IIE_M1"]] <- rv.value(med_cwc, mus_IIE_M1, alphas_each_IIE_M1_r, m_each_IIE_M1, params = "IIE_M1", theta.null = 0, cf = seq(0, 0.99, by = 0.01), only_YM_confounded = FALSE, mediatedY = 1:2)
      
      rv_res[["IIE_M1"]] <- rv.value(med1_cwc, mus_IIE_M1, alphas_each_IIE_M1_r, m_each_IIE_M1, params = "IIE_M1", theta.null = 0, cf = seq(0, 0.99, by = 0.01), only_YM_confounded = FALSE, mediatedY = 1:2)
      
      # IIE_M2
      sens_fit_res[["IIE_M2"]] <- sens_IIE_Mjo(med2_cwc, mus_IIE_M2, alphas_each_IIE_M2_r, m_each_IIE_M2, params = "IIE_M2", cf.y, cf.m, rho2, only_YM_confounded = FALSE, mediatedY = 1:2)
      
      rv_res[["IIE_M2"]] <- rv.value(med2_cwc, mus_IIE_M2, alphas_each_IIE_M2_r, m_each_IIE_M2, params = "IIE_M2", theta.null = 0, cf = seq(0, 0.99, by = 0.01), only_YM_confounded = FALSE, mediatedY = 1:2)
      
    }
    
    # IDE
    sens_fit_res[["IDE"]] <- sens_IIE_Mjo(med_cwc, mus_Mjo[3:4], alphas_each_IDE_r, m_each_Mjo[3:4], params = "IDE", cf.y, cf.m, rho2, only_YM_confounded = FALSE, mediatedY = 1:2)
    
    rv_res[["IDE"]] <- rv.value(med_cwc, mus_Mjo[3:4], alphas_each_IDE_r, m_each_Mjo[3:4], params = "IDE", theta.null = 0, cf = seq(0, 0.99, by = 0.01), only_YM_confounded = FALSE, mediatedY = 1:2)
    
    # IIE_Mjo
    
    sens_fit_res[["IIE_Mjo"]] <- sens_IIE_Mjo(med_cwc, mus_Mjo[1:2], alphas_each_IIE_Mjo_r, m_each_Mjo[1:2], params = "IIE_Mjo", cf.y, cf.m, rho2, only_YM_confounded = FALSE, mediatedY = 1:2)
    
    rv_res[["IIE_Mjo"]] <- rv.value(med_cwc, mus_Mjo[1:2], alphas_each_IIE_Mjo_r, m_each_Mjo[1:2], params = "IIE_Mjo", theta.null = 0, cf = seq(0, 0.99, by = 0.01), only_YM_confounded = FALSE, mediatedY = 1:2)

    ## Minimal sensitivity reporting ----
    rv_values <- map_dbl(rv_res, \(x) x$RV_null0_sig0.05)
    tab_sensitivity_reporting <- estimates %>%
        dplyr::filter(str_detect(estimand, "I")) %>%
      mutate(CI = glue("[{round(`95% CI_low`, 2)}, {round(`95% CI_high`, 2)}]")) %>% 
        mutate(RV_null0_sig0.05 = rv_values[estimand])

    # tab_sensitivity_reporting %>% write_csv(file = "Example/Results/tab_sensitivity_reporting.csv")
    # tab_sensitivity_reporting %>% write_csv(file = "Example/Results_NoCluster/tab_sensitivity_reporting.csv")
    
    
    
    # Benchmark ------------------------
    # benchmark_covariates <- list(a=c("X_BYSES1","X_BYPARED","X_BYINCOME", "X_BYRACE_Hispanic", "X_BYRACE_White", "X_BYRACE_Other", "X_BYGNSTAT", "X_BYP49", "X_BYP08", "X_BYGNSTAT", "X_BYHOMLIT", "X_BYSEX"))
    benchmark_covariates <- list(a=c("X_BYSEX")) #"X_BYPARED","X_BYINCOME", "X_BYP44C", 
    Benchmark_res <- list()
    set.seed(12)
    Benchmark_res[["IIE_Mjo"]] <- benchmark.effect(
      med_cwc,
      benchmark_covariates, param = "IIE_Mjo", mus_Mjo[1:2], alphas_each_IIE_Mjo_r, m_each_Mjo[1:2], estimates, folds, control, learners_outcome, only_YM_confounded = FALSE, mediatedY = 1:2
    )
    
    Benchmark_res[["IIE_Mjo"]][[1]][1:2]
    # sum(Benchmark_res_IIE_M1_tmp[[1]]$Gain_mu * Benchmark_res_IIE_M1_tmp[[1]]$Gain_alpha * unlist(Benchmark_res_IIE_M1_tmp[[1]]$rho))
    set.seed(12)
    Benchmark_res[["IDE"]] <- benchmark.effect(
      med_cwc,
      benchmark_covariates, param = "IDE", mus_Mjo[3:4], alphas_each_IDE_r, m_each_Mjo[3:4], estimates, folds, control, learners_outcome, only_YM_confounded = FALSE, mediatedY = 1:2
    )
    
    Benchmark_res[["IDE"]][[1]][1:2]
    
    if (num_M > 1) {
      # benchmark_covariates = list(a=c("X_BYSES1","X_BYPARED","X_BYINCOME"))
      
      Benchmark_res_IIE_M1_tmp <- benchmark.effect(
        med1_cwc,
        benchmark_covariates, param = "IIE_M1", mus_IIE_M1, alphas_each_IIE_M1_r, m_each_IIE_M1, estimates, folds, control, learners_outcome, only_YM_confounded = FALSE, mediatedY = 1:2
      )
      sum(Benchmark_res_IIE_M1_tmp[[1]]$Gain_mu * Benchmark_res_IIE_M1_tmp[[1]]$Gain_alpha * unlist(Benchmark_res_IIE_M1_tmp[[1]]$rho))
      
      Benchmark_res_IIE_M2_tmp <- benchmark.effect(
        med2_cwc,
        benchmark_covariates, param = "IIE_M2", mus_IIE_M2, alphas_each_IIE_M2_r, m_each_IIE_M2, estimates, folds, control, learners_outcome, only_YM_confounded = FALSE, mediatedY = 1:2
      )
    }
    
    
    # write.table(Benchmark_res_IIE_M1_tmp, file = "Results/Benchmark_res_IIE_M1.txt", append = TRUE)
    
    ## Contour plots -----
    sens_grid <- list()

    cf_grid <- expand.grid(cf.y = seq(0.01, 0.99, by = 0.01), cf.m = seq(0.01, 0.99, by = 0.01))
    label.unadjusted <- "Original"
    
    if (num_M > 1) {
      ### IIE_M1 contour ----
      sens_grid[["IIE_M1"]] <- map2_df(
        .x = cf_grid$cf.m, .y = cf_grid$cf.y,
        .f = ~ sens_IIE_Mjo(med1_cwc, mus_IIE_M1, alphas_each_IIE_M1_r, m_each_IIE_M1, params = "IIE_M1", 
                            cf.y = .y, cf.m = .x, 
                            rho2 = 1, only_YM_confounded = FALSE, mediatedY = 1:2)
      )
      
      
      z_grid <- sens_grid[["IIE_M1"]] %>% 
        select(cf.y, cf.m, `theta_m.90%CI_low`) %>% 
        pivot_wider(names_from = cf.m, values_from = `theta_m.90%CI_low`) %>% 
        as.matrix(.)
      z_grid <- z_grid[, -1]
      
      bench.cf.y <- Benchmark_res_IIE_M1_tmp[[1]]$Gain_mu 
      bench.cf.m <- Benchmark_res_IIE_M1_tmp[[1]]$Gain_alpha
      add_benchmark_max <- sens_IIE_Mjo(med2_cwc, mus_IIE_M1, alphas_each_IIE_M1_r, m_each_IIE_M1, params = "IIE_M1", 
                                        cf.y = max(bench.cf.y), cf.m = max(bench.cf.m), 
                                        rho2 = 1, only_YM_confounded = FALSE, mediatedY = 1:2) #[["theta_m.90%CI_low"]]
      # add_benchmark_max$`theta_m.90%CI_low`
      add_benchmark_min <- sens_IIE_Mjo(med2_cwc, mus_IIE_M1, alphas_each_IIE_M1_r, m_each_IIE_M1, params = "IIE_M1", 
                                        cf.y = min(bench.cf.y), cf.m = min(bench.cf.m), 
                                        rho2 = 1, only_YM_confounded = FALSE, mediatedY = 1:2) #[["theta_m.90%CI_low"]]
      
      plot_estimate <- rv_res[["IIE_M1"]]$RVsens_res %>% 
        filter(cf.y == 0, cf.m == 0) %>% 
        pull(conf_bound.l)
      
      pdf(file = "Example/Results/contour_IIE_M1.pdf", width = 9, height = 9)
      
      contour_plot(grid_values.x = unique(cf_grid$cf.m),
                   grid_values.y = unique(cf_grid$cf.y),
                   z_axis = z_grid,
                   levels = c(-6, -3.5, -2, -1.5, -1, -0.8, -0.6, -0.4, -0.2,  0, 0.1),
                   labels = NULL,
                   threshold = 0,
                   col.contour = "blue4")
      
      points(0, 0, pch = 17, col = "black", cex = 0.8)
      label.bump.x <- 0.03
      label.bump.y <- 0.05
      text(0 + label.bump.x, 0 + label.bump.y,
           paste0(label.unadjusted,
                  "\n(", round(plot_estimate, 3), ")"),
           cex = 0.7)
      
      sensemakr::add_bound_to_contour(r2dz.x = add_benchmark_max$cf.m,
                                      r2yz.dx = add_benchmark_max$cf.y,
                                      bound_value = add_benchmark_max$`theta_m.90%CI_low`,
                                      bound_label = "Maximum \nBenchmark",
                                      label.text = TRUE,
                                      label.bump.x = 0.03,
                                      label.bump.y = 0.04,
                                      cex.label.text = 0.7,
                                      round = 3)
      
      sensemakr::add_bound_to_contour(r2dz.x = add_benchmark_min$cf.m,
                                      r2yz.dx = add_benchmark_min$cf.y,
                                      bound_value = add_benchmark_min$`theta_m.90%CI_low`,
                                      bound_label = "Minimum \nBenchmark",
                                      label.text = TRUE,
                                      label.bump.x = 0.07,
                                      label.bump.y = 0.01,
                                      cex.label.text = 0.7,
                                      round = 3)
      
      dev.off()
      
      ### IIE_M2 contour ----
      sens_grid[["IIE_M2"]] <- map2_df(
        .x = cf_grid$cf.m, .y = cf_grid$cf.y,
        .f = ~ sens_IIE_Mjo(med2_cwc, mus_IIE_M2, alphas_each_IIE_M2_r, m_each_IIE_M2, params = "IIE_M2", 
                            cf.y = .y, cf.m = .x, 
                            rho2 = 1, only_YM_confounded = FALSE, mediatedY = 1:2)
      )
      
      
      z_grid <- sens_grid[["IIE_M2"]] %>% 
        select(cf.y, cf.m, `theta_m.90%CI_low`) %>% 
        pivot_wider(names_from = cf.m, values_from = `theta_m.90%CI_low`) %>% 
        as.matrix(.)
      z_grid <- z_grid[, -1]
      
      bench.cf.y <- Benchmark_res_IIE_M2_tmp[[1]]$Gain_mu 
      bench.cf.m <- Benchmark_res_IIE_M2_tmp[[1]]$Gain_alpha
      add_benchmark_max <- sens_IIE_Mjo(med2_cwc, mus_IIE_M2, alphas_each_IIE_M2_r, m_each_IIE_M2, params = "IIE_M2", 
                                        cf.y = max(bench.cf.y), cf.m = max(bench.cf.m), 
                                        rho2 = 1, only_YM_confounded = FALSE, mediatedY = 1:2) #[["theta_m.90%CI_low"]]
      # add_benchmark_max$`theta_m.90%CI_low`
      add_benchmark_min <- sens_IIE_Mjo(med2_cwc, mus_IIE_M2, alphas_each_IIE_M2_r, m_each_IIE_M2, params = "IIE_M2", 
                                        cf.y = min(bench.cf.y), cf.m = min(bench.cf.m), 
                                        rho2 = 1, only_YM_confounded = FALSE, mediatedY = 1:2) #[["theta_m.90%CI_low"]]
      
      plot_estimate <- rv_res[["IIE_M2"]]$RVsens_res %>% 
        filter(cf.y == 0, cf.m == 0) %>% 
        pull(conf_bound.l)
      
      pdf(file = "Example/Results/contour_IIE_M2.pdf", width = 9, height = 9)
      
      contour_plot(grid_values.x = unique(cf_grid$cf.m),
                   grid_values.y = unique(cf_grid$cf.y),
                   z_axis = z_grid,
                   levels = c( -3.5, -2, -1.5, -1,   -0.75, -0.5, -0.25, -0.1, 0),
                   labels = NULL,
                   threshold = 0,
                   col.contour = "blue4")
      
      points(0, 0, pch = 17, col = "black", cex = 0.8)
      label.bump.x <- 0.02
      label.bump.y <- 0.05
      text(0 + label.bump.x, 0 + label.bump.y,
           paste0(label.unadjusted,
                  "\n(", round(plot_estimate, 3), ")"),
           cex = 0.7)
      
      sensemakr::add_bound_to_contour(r2dz.x = add_benchmark_max$cf.m,
                                      r2yz.dx = add_benchmark_max$cf.y,
                                      bound_value = add_benchmark_max$`theta_m.90%CI_low`,
                                      bound_label = "Maximum \nBenchmark",
                                      label.text = TRUE,
                                      label.bump.x = 0.03,
                                      label.bump.y = 0.04,
                                      cex.label.text = 0.7,
                                      round = 3)
      
      sensemakr::add_bound_to_contour(r2dz.x = add_benchmark_min$cf.m,
                                      r2yz.dx = add_benchmark_min$cf.y,
                                      bound_value = add_benchmark_min$`theta_m.90%CI_low`,
                                      bound_label = "Minimum \nBenchmark",
                                      label.text = TRUE,
                                      label.bump.x = 0.07,
                                      label.bump.y = 0.01,
                                      cex.label.text = 0.7,
                                      round = 3)
      
      dev.off()
      
    }
    
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
    
    bench.cf.y <- Benchmark_res[["IDE"]][[1]]$Gain_mu 
    bench.cf.m <- Benchmark_res[["IDE"]][[1]]$Gain_alpha
    add_benchmark_max <- sens_IIE_Mjo(med, mus_Mjo[3:4], alphas_each_IDE_r, m_each_Mjo[3:4], params = "IDE", 
                                      cf.y = max(bench.cf.y), cf.m = max(bench.cf.m), 
                                      rho2 = 1, only_YM_confounded = FALSE, mediatedY = 1:2) #[["theta_m.90%CI_low"]]
    # add_benchmark_max$`theta_m.90%CI_low`
    add_benchmark_min <- sens_IIE_Mjo(med, mus_Mjo[3:4], alphas_each_IDE_r, m_each_Mjo[3:4], params = "IDE", 
                                      cf.y = min(bench.cf.y), cf.m = min(bench.cf.m), 
                                      rho2 = 1, only_YM_confounded = FALSE, mediatedY = 1:2) #[["theta_m.90%CI_low"]]
    
    plot_estimate <- rv_res[["IDE"]]$RVsens_res %>% 
      filter(cf.y == 0, cf.m == 0) %>% 
      pull(conf_bound.l)
    
    quantile(z_grid)
    # pdf(file = "Example/Results/contour_IDE.pdf", width = 5, height = 5)
    pdf(file = "Example/Results_NoCluster/contour_IDE_nocluster.pdf", width = 10, height = 10)
    
    contour_plot(grid_values.x = unique(cf_grid$cf.m),
                 grid_values.y = unique(cf_grid$cf.y),
                 z_axis = z_grid,
                 levels = c(-10, -7, -5, -4, -3,  -2,  -1,  -0.5, 0, 0.6),
                 labels = NULL,
                 threshold = 0,
                 col.contour = "blue4")
    
    
    points(0, 0, pch = 17, col = "black", cex = 0.8)
    label.bump.x <- 0.02
    label.bump.y <- 0.03
    text(0 + label.bump.x, 0 + label.bump.y,
         paste0(label.unadjusted,
                "\n(", round(plot_estimate, 3), ")"),
         cex = 0.6)
    
    sensemakr::add_bound_to_contour(r2dz.x = add_benchmark_max$cf.m,
                                    r2yz.dx = add_benchmark_max$cf.y,
                                    bound_value = add_benchmark_max$`theta_m.90%CI_low`,
                                    bound_label = "Maximum \nBenchmark",
                                    label.text = TRUE,
                                    label.bump.x = 0.07,
                                    label.bump.y = 0.04,
                                    cex.label.text = 0.7,
                                    round = 3)
    
    sensemakr::add_bound_to_contour(r2dz.x = add_benchmark_min$cf.m,
                                    r2yz.dx = add_benchmark_min$cf.y,
                                    bound_value = add_benchmark_min$`theta_m.90%CI_low`,
                                    bound_label = "Minimum \nBenchmark",
                                    label.text = TRUE,
                                    label.bump.x = 0.06,
                                    label.bump.y = 0,
                                    cex.label.text = 0.7,
                                    col = "green4",
                                    round = 3)
    title(sub = "Benchmark covariate: Sex")
    
    dev.off()
    
    ### IIE_Mjo contour -----
    sens_grid[["IIE_Mjo"]] <- map2_df(
      .x = cf_grid$cf.m, .y = cf_grid$cf.y,
      .f = ~ sens_IIE_Mjo(med, mus_Mjo[1:2], alphas_each_IIE_Mjo_r, m_each_Mjo[1:2],
                          params = "IIE_Mjo",
                          cf.y = .y, cf.m = .x, rho2 = 1
      )
    )
    
    
    z_grid <- sens_grid[["IIE_Mjo"]] %>% 
      select(cf.y, cf.m, `theta_m.90%CI_low`) %>% 
      pivot_wider(names_from = cf.m, values_from = `theta_m.90%CI_low`) %>% 
      as.matrix(.)
    z_grid <- z_grid[, -1]
    
    bench.cf.y <- Benchmark_res[["IIE_Mjo"]][[1]]$Gain_mu 
    bench.cf.m <- Benchmark_res[["IIE_Mjo"]][[1]]$Gain_alpha
    add_benchmark_max <- sens_IIE_Mjo(med, mus_Mjo[1:2], alphas_each_IIE_Mjo_r, m_each_Mjo[1:2], params = "IIE_Mjo", 
                                      cf.y = max(bench.cf.y), cf.m = max(bench.cf.m), 
                                      rho2 = 1, only_YM_confounded = FALSE, mediatedY = 1:2) #[["theta_m.90%CI_low"]]
    # add_benchmark_max$`theta_m.90%CI_low`
    add_benchmark_min <- sens_IIE_Mjo(med, mus_Mjo[1:2], alphas_each_IIE_Mjo_r, m_each_Mjo[1:2], params = "IIE_Mjo", 
                                      cf.y = min(bench.cf.y), cf.m = min(bench.cf.m), 
                                      rho2 = 1, only_YM_confounded = FALSE, mediatedY = 1:2) #[["theta_m.90%CI_low"]]
    
    plot_estimate <- rv_res[["IIE_Mjo"]]$RVsens_res %>% 
      filter(cf.y == 0, cf.m == 0) %>% 
      pull(conf_bound.l)
    
    quantile(z_grid)
    # pdf(file = "Example/Results/contour_IIE_Mjo.pdf", width = 5, height = 5)
    pdf(file = "Example/Results_NoCluster/contour_IIE_Mjo_nocluster.pdf", width = 10, height = 10)
    
    contour_plot(grid_values.x = unique(cf_grid$cf.m),
                 grid_values.y = unique(cf_grid$cf.y),
                 z_axis = z_grid,
                 levels = c(-10, -5, -3, -2,  -1,  -0.5, 0, 0.2),
                 labels = NULL,
                 threshold = 0,
                 col.contour = "blue4")
    
    points(0, 0, pch = 17, col = "black", cex = 0.8)
    label.bump.x <- 0.01
    label.bump.y <- 0.05
    text(0 + label.bump.x, 0 + label.bump.y,
         paste0(label.unadjusted,
                "\n(", round(plot_estimate, 3), ")"),
         cex = 0.6)
    
    sensemakr::add_bound_to_contour(r2dz.x = add_benchmark_max$cf.m,
                                    r2yz.dx = add_benchmark_max$cf.y,
                                    bound_value = add_benchmark_max$`theta_m.90%CI_low`,
                                    bound_label = "Maximum \nBenchmark",
                                    label.text = TRUE,
                                    label.bump.x = 0.05,
                                    label.bump.y = 0.04,
                                    cex.label.text = 0.7,
                                    round = 3)
    
    sensemakr::add_bound_to_contour(r2dz.x = add_benchmark_min$cf.m,
                                    r2yz.dx = add_benchmark_min$cf.y,
                                    bound_value = add_benchmark_min$`theta_m.90%CI_low`,
                                    bound_label = "Minimum \nBenchmark",
                                    label.text = TRUE,
                                    label.bump.x = 0.07,
                                    label.bump.y = 0,
                                    cex.label.text = 0.7,
                                    round = 3)
    title(sub = "Benchmark covariate: Sex")
    
    dev.off()
    
    
    # output ----------------
    med_sens <- list(estimates, rv_res, sens_fit_res)
    names(med_sens) <- c("estimates_short", "robustness_value", "sens_fit_res")
    out <- mget(ls(envir = environment()))

    out
    # med_sens
}
