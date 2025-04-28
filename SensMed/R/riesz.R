
riesz_crossfit <- function(med, m_each, folds, which.M, learners_riesz, control, riesz_alpha1 = FALSE, riesz_alpha2 = FALSE, riesz_alpha3 = FALSE, alpha_ac, a_c, a_mc, a_m1c, a_m2c, r_zmac, covariates_cl = FALSE) {
    cf_list <- vector("list", control$crossfit_folds)
    i <- 1
    cli::cli_progress_step("Computing alphas ... {i}/{control$crossfit_folds} folds")
    for (i in seq_along(folds)) {
        # train <- training(med, folds, i)
        # valid <- validation(med, folds, i)
        # m_each1 <- m_each[[1]]
        cf_list[[i]] <- lapply(
            m_each,
            \(m_each1) riesz_estimate(med, folds, i_fold = i, m_each1, which.M, learners_riesz, control, riesz_alpha1, riesz_alpha2, riesz_alpha3, alpha_ac, a_c, a_mc, a_m1c, a_m2c, r_zmac, covariates_cl)
        )
        
        names(cf_list[[i]]) <- map_chr(m_each, ~paste(glue( '{names(.)}{gsub("data_", "", .)}' ), collapse = "_") )
        #unlist(lapply(m_each, \(x) paste0(gsub("data_", "", x), collapse = "")))
        cli::cli_progress_update()
    }
    cli::cli_progress_done()
    
    # cf_list[[1]]$az0_indepM_am0_indepM_ay1_indepM
    alpha_allfolds <- combine_folds_alpha(cf_list, folds)
    alpha_allfolds
}

riesz_estimate <- function(med, folds, i_fold, m_each1, which.M = NULL, learners_riesz, control, riesz_alpha1 = FALSE, riesz_alpha2 = FALSE, riesz_alpha3 = FALSE, alpha_ac, a_c, a_mc, a_m1c, a_m2c, r_zmac, covariates_cl = FALSE) {
    
    train <- training(med, folds, i_fold)
    valid <- validation(med, folds, i_fold)
    
    vars <- med@vars
    # for alpha1
    # varnames_1 <- c(glue("{c(vars@A, covariates_cl$Xij)}"), glue("{c(vars@A, vars@M[-which.M])}_clmean")
    #                # covariates_cl$Wj, covariates_cl$Sdumm, glue("{vars@M[-which.M]}_clmean")
    # )
    # varnames_2 <- c(glue("{c(vars@A, covariates_cl$Xij, vars@M[-which.M])}"), glue("{c(vars@A, vars@M)}_clmean")
    #                 #covariates_cl$Wj, covariates_cl$Sdumm, glue("{vars@M[which.M]}_clmean")
    # )
    # varnames_y <- c(glue("{c(vars@A, covariates_cl$Xij, vars@M)}"), glue("{c(vars@A, vars@M, vars@Y)}_clmean")
    #                 #covariates_cl$Wj, covariates_cl$Sdumm, glue("{vars@Y}_clmean")
    # )
    
    varnames_1 <- c(vars@A, vars@C)
    if (!is.null(which.M)) {
        varnames_2 <- c(vars@A, vars@C, vars@M[-which.M])
    }
    varnames_y <- c(vars@A, vars@C, vars@M)
    
    if (covariates_cl) {
        varnames_1 <- na.omit(c(glue("{c(vars@M, vars@Y)}_clmean"), glue("{c(vars@A, vars@C)}_cwc"), vars@Wj))
        if ( !is.null(which.M) ) {
            varnames_2 <- na.omit(c(c(glue("{c(vars@M[which.M], vars@Y)}_clmean"), glue("{c(vars@A, vars@C, vars@M[-which.M])}_cwc"), vars@Wj)))
        }
        varnames_y <- na.omit(c(c(glue("{c(vars@Y)}_clmean"), glue("{c(vars@A, vars@C, vars@M)}_cwc"), vars@Wj)))
    }
    
    # varnames_1 <- c(vars@A, gsub("A", "A0", vars@A), vars@C)
    # varnames_2 <- c(vars@A, gsub("A", "A0", vars@A), vars@C, vars@M[-which.M])
    # varnames_y <- c(vars@A, gsub("A", "A0", vars@A), vars@C, vars@M)
    
    # if (str_detect(m_each1["am"], "_0")) {
    #     varnames_1 <- gsub("A", "A0", varnames_1)
    # }
    # if (str_detect(m_each1["az"], "_0")) {
    #     varnames_2 <- gsub("A", "A0", varnames_2)
    # }
    # if (str_detect(m_each1["ay"], "_0")) {
    #     varnames_y <- gsub("A", "A0", varnames_y)
    # }
    
    # alpha1 ----
    if (riesz_alpha1 == TRUE) {
        alpha1 <- riesz_fit(
            natural_train = train[[ "data" ]][, varnames_1],
            natural_valid = valid[[ "data" ]][, varnames_1],
            shifted_train = list("shifted" = train[[ m_each1["am"] ]][, varnames_1]),
            shifted_valid = list("shifted" = valid[[ m_each1["am"] ]][, varnames_1]),
            prev_riesz_train = data.frame(alpha_prev_train = rep(1, nrow(train[[ "data" ]]))),
            learners_riesz = learners_riesz,
            riesz_folds = control$riesz_folds
        )
        
        # # NOT WORK
        # alpha1 <- riesz_fit(
        #     natural_train = train[[ "data" ]][, varnames_1],
        #     natural_valid = valid[[ "data" ]][, varnames_1],
        # 
        #     shifted_train = list(
        #         "a0" = train[[ m_each1["am.m"] ]][, varnames_1],
        #         "a1" = train[[ m_each1["am.p"] ]][, varnames_1]
        #     ),
        #     shifted_valid = list(
        #         "a0" = valid[[ m_each1["am.m"] ]][, varnames_1],
        #         "a1" = valid[[ m_each1["am.p"] ]][, varnames_1]
        #     ),
        #     prev_riesz_train = data.frame(alpha_prev_train = rep(1, nrow(train[[ "data" ]]))),
        #     learners_riesz = learners_riesz,
        #     riesz_folds = control$riesz_folds
        # )
        
    } 
    
    
    if (riesz_alpha1 == FALSE) {
        a_val <- substr(m_each1["am"], 6, 6)
        alpha1 <- list()
        alpha1$predictions_train <- alpha_ac[[ glue("A={a_val}") ]][folds[[i_fold]]$training_set]
        
        alpha1$predictions_valid <- alpha_ac[[ glue("A={a_val}") ]][folds[[i_fold]]$validation_set]
        
        alpha1$predictions_shifted_valid <- alpha.ac(
            A = valid[[ m_each1["am"] ]][[med@vars@A]], 
            a_c = a_c[folds[[i_fold]]$validation_set, ])[[ glue("A={a_val}") ]]
        
    }
    # mean(alpha1$predictions_train * train[[ "data" ]]$Y)
    
    # alpha2 ----
    if (m_each1[["az"]] != "") {
        if (riesz_alpha2 == TRUE) {
            alpha2 <- riesz_fit(
                natural_train = train[[ "data" ]][, varnames_2],
                natural_valid = valid[[ "data" ]][, varnames_2],
                shifted_train = list("shifted" = train[[ m_each1["az"] ]][, varnames_2]),
                shifted_valid = list("shifted" = valid[[ m_each1["az"] ]][, varnames_2]),
                prev_riesz_train = data.frame(alpha_prev_train = alpha1$predictions_train),
                learners_riesz = learners_riesz, 
                riesz_folds = control$riesz_folds
            )
        }
        
        if (riesz_alpha2 == FALSE) {
            az <- substr(m_each1["az"], 6, 6)
            am <- substr(m_each1["am"], 6, 6)
            
            alpha2 <- list()
            alpha2$predictions_train <- alpha.amc(
                A = train[[ "data" ]][[med@vars@A]], 
                ay = az, am = am, 
                a_c = a_c[folds[[i_fold]]$training_set, ],
                a_mc = a_m1c[folds[[i_fold]]$training_set, ])
            
            alpha2$predictions_valid <- alpha.amc(
                A = valid[[ "data" ]][[med@vars@A]], 
                ay = az, am = am, 
                a_c = a_c[folds[[i_fold]]$validation_set, ],
                a_mc = a_m1c[folds[[i_fold]]$validation_set, ])
            
            
            alpha2$predictions_shifted_valid <- alpha.amc(
                A = valid[[ m_each1["az"] ]][[med@vars@A]], 
                ay = az, am = am, 
                a_c = a_c[folds[[i_fold]]$validation_set, ],
                a_mc = a_m1c[folds[[i_fold]]$validation_set, ])
        }
    }
    if (m_each1[["az"]] == "") {
        alpha2 <- alpha1
    }
    
    # alpha3 ---------
    if (riesz_alpha3 == TRUE) {
        alpha3 <- riesz_fit(
            natural_train = train[[ "data" ]][, varnames_y],
            natural_valid = valid[[ "data" ]][, varnames_y],
            shifted_train = list("shifted" = train[[ m_each1["ay"] ]][, varnames_y]),
            shifted_valid = list("shifted" = valid[[ m_each1["ay"] ]][, varnames_y]),
            prev_riesz_train = data.frame(alpha_prev_train = alpha2$predictions_train),
            learners_riesz = learners_riesz, 
            riesz_folds = control$riesz_folds
        )
    }
    if (riesz_alpha3 == FALSE) {
        if (m_each1[["az"]] == "") {
            ay <- substr(m_each1["ay"], 6, 6)
            am <- substr(m_each1["am"], 6, 6)
            alpha3 <- list()
            alpha3$predictions_train <- alpha.amc(
                A = train[[ "data" ]][[med@vars@A]], 
                ay = ay, am = am, 
                a_c = a_c[folds[[i_fold]]$training_set, ],
                a_mc = a_mc[folds[[i_fold]]$training_set, ])
            
            alpha3$predictions_valid <- alpha.amc(
                A = valid[[ "data" ]][[med@vars@A]], 
                ay = ay, am = am, 
                a_c = a_c[folds[[i_fold]]$validation_set, ],
                a_mc = a_mc[folds[[i_fold]]$validation_set, ])
            
            
            alpha3$predictions_shifted_valid <- alpha.amc(
                A = valid[[ m_each1["ay"] ]][[med@vars@A]], 
                ay = ay, am = am, 
                a_c = a_c[folds[[i_fold]]$validation_set, ],
                a_mc = a_mc[folds[[i_fold]]$validation_set, ])
        }
        if (m_each1[["az"]] != "") {
            ay <- substr(m_each1["ay"], 6, 6)
            am <- substr(m_each1["am"], 6, 6)
            az <- substr(m_each1["az"], 6, 6)
            alpha3 <- list()
            alpha3$predictions_train <- alpha.azmc(
                A = train[[ "data" ]][[med@vars@A]], 
                ay = ay, am = am, az = az,
                a_c = a_c[folds[[i_fold]]$training_set, ],
                a_mc = a_m1c[folds[[i_fold]]$training_set, ],
                a_zc = a_m2c[folds[[i_fold]]$training_set, ],
                r_zmac = r_zmac[folds[[i_fold]]$training_set, ])
            
            alpha3$predictions_valid <- alpha.azmc(
                A = valid[[ "data" ]][[med@vars@A]], 
                ay = ay, am = am, az = az,
                a_c = a_c[folds[[i_fold]]$validation_set, ],
                a_mc = a_m1c[folds[[i_fold]]$validation_set, ],
                a_zc = a_m2c[folds[[i_fold]]$validation_set, ],
                r_zmac = r_zmac[folds[[i_fold]]$validation_set, ])
            
            
            alpha3$predictions_shifted_valid <- alpha.azmc(
                A = valid[[ m_each1["ay"] ]][[med@vars@A]], 
                ay = ay, am = am, az = az,
                a_c = a_c[folds[[i_fold]]$validation_set, ],
                a_mc = a_m1c[folds[[i_fold]]$validation_set, ],
                a_zc = a_m2c[folds[[i_fold]]$validation_set, ],
                r_zmac = r_zmac[folds[[i_fold]]$validation_set, ])
        }
    }
    
    if (m_each1[["az"]] != "") {
        alpha_list_valid <- list(#zmy = paste(substr(m_each1, 6, 6), collapse = ""),
            alpha1 = alpha1$predictions_valid,
            alpha2 = alpha2$predictions_valid,
            alpha3 = alpha3$predictions_valid,
            alpha1_target = alpha1$predictions_shifted_valid,
            alpha2_target = alpha2$predictions_shifted_valid,
            alpha3_target = alpha3$predictions_shifted_valid)
        
    }
    if (m_each1[["az"]] == "") {
        alpha_list_valid <- list(#zmy = paste(substr(m_each1, 6, 6), collapse = ""),
            alpha1 = alpha1$predictions_valid,
            # alpha2 = alpha2$predictions_valid,
            alpha2 = alpha3$predictions_valid,
            alpha1_target = alpha1$predictions_shifted_valid,
            # alpha2_target = alpha2$predictions_shifted_valid,
            alpha2_target = alpha3$predictions_shifted_valid)
        
    }
    
    
    alpha_list_valid
} 

# remotes::install_github("herbps10/SuperRiesz")
# from lmtp R package  https://github.com/nt-williams/lmtp/tree/riesz 
#' @importFrom SuperRiesz super_riesz

riesz_fit <- function(natural_train, shifted_train, natural_valid, shifted_valid, prev_riesz_train, learners_riesz = c("nn", "linear"), riesz_folds = 5) {
    
    natural_train <- one_hot_encode(natural_train)
    natural_valid <- one_hot_encode(natural_valid)
    if (length(shifted_train) == 1) {
        shifted_train_list <- map(shifted_train, one_hot_encode)
        shifted_valid_list <- map(shifted_valid, one_hot_encode)
        names(shifted_train_list) <- "shifted"
        names(shifted_valid_list) <- "shifted"
        
        m <- \(alpha, data) {
            alpha(data("shifted")) * data("weight")[,1]
        }
    }
    if (length(shifted_train) > 1) {
        shifted_train_list <- map(shifted_train, one_hot_encode)
        # names(shifted_train)
        names(shifted_train_list) <- c("a0", "a1") # be the same as names(shifted_train)
        shifted_valid_list <- map(shifted_valid, one_hot_encode)
        names(shifted_valid_list) <- c("a0", "a1") 
        
        m <- \(alpha, data) {
            (alpha(data("a1")) - alpha(data("a0"))) #* data("weight")[,1]
        }
    }
    
    sl <- SuperRiesz::super_riesz(
        data = natural_train,
        library = learners_riesz,
        alternatives = shifted_train_list,
        extra = list(weight = data.frame(weight = prev_riesz_train[, 1])),
        m = m,
        folds = riesz_folds
    )
    
    predictions_valid <- predict(sl, natural_valid)
    predictions_train <- predict(sl, natural_train)
    if (length(shifted_train) == 1) {
        predictions_shifted_valid <- predict(sl, shifted_valid_list[["shifted"]]) 
        predictions_shifted_train <- predict(sl, shifted_train_list[["shifted"]]) 
    }
    
    if (length(shifted_train) == 2) {
        predictions_shifted_valid <- predict(sl, shifted_valid_list[["a1"]]) - predict(sl, shifted_valid_list[["a0"]])
        predictions_shifted_train <- predict(sl, shifted_train_list[["a1"]]) - predict(sl, shifted_train_list[["a0"]])
    }
    
    list(
        predictions_valid = predictions_valid,
        predictions_train = predictions_train,
        predictions_shifted_valid = predictions_shifted_valid,
        predictions_shifted_train = predictions_shifted_train,
        fits = sl,
        coef = sl$weights,
        risk = sl$risk
    )

}


