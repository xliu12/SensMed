
estimate_mu.long.short <- function(benchmark_covariates, 
                                   med, m_each, folds, learners_regressions, control, which.M, covariates_cl = FALSE) {
    mu_folds <- vector("list", control$crossfit_folds)
    i <- 1
    cli::cli_progress_step("Fitting outcome regressions... {i}/{control$crossfit_folds} folds")
    for (i in seq_along(mu_folds)) {
        train <- training(med, folds, i)
        valid <- validation(med, folds, i)
        mu_folds[[i]] <- mu_eachM.long.short(benchmark_covariates, 
                                  train, valid, med@vars, m_each,  learners_regressions, control, which.M, covariates_cl)
        cli::cli_progress_update()
    }
    cli::cli_progress_done()
    mu_allfolds <- combine_folds_theta(mu_folds, folds)
    mu_allfolds
}

mu_eachM.long.short <- function(benchmark_covariates, 
                                train, valid, vars, m_each, learners, control, which.M, covariates_cl = FALSE) {
    
    varnames_1 <- c(vars@A, vars@C)
    varnames_1_short <- setdiff(varnames_1, benchmark_covariates)
    
    if ( !is.null(which.M) ) {
        varnames_2 <- c(vars@A, vars@C, vars@M[-which.M])
        varnames_2_short <- setdiff(varnames_2, benchmark_covariates)
    }
    varnames_y <- c(vars@A, vars@C, vars@M)
    varnames_y_short <- setdiff(varnames_y, benchmark_covariates)
    
    y_type <- ifelse(if.binary(train$data[[vars@Y]]), "binomial", "continuous")
    
    mu_y_short <- mlr3superlearner::mlr3superlearner(
        data = train$data[, na.omit(c(varnames_y_short, vars@Y))],
        target = vars@Y,
        library = learners,
        outcome_type = y_type,
        folds = control$mlr3superlearner_folds,
        newdata = valid,
        group = NULL,
        discrete = TRUE
    )
    mu_y_long <- mlr3superlearner::mlr3superlearner(
        data = train$data[, na.omit(c(varnames_y, vars@Y))],
        target = vars@Y,
        library = learners,
        outcome_type = y_type,
        folds = control$mlr3superlearner_folds,
        newdata = valid,
        group = NULL,
        discrete = TRUE
    )
    
    mu_list <- vector("list", length = length(m_each))
    
    i <- 1
    for (i in 1:length(m_each)) {
        # m_each[[i]][L+1]
        m_each1 <- m_each[[i]]
        
        if (m_each1[["az"]] != "") {
            
            b_list_train <- b_list_valid <- mu_fit <- vector("list", length = 3) 
            names(b_list_valid) <- glue("b{1:length(b_list_valid)}")
            names(mu_fit) <- glue("mu_fit{1:length(mu_fit)}")
            
            mu_fit[[3]] <- mu_y_long$preds$data
            
            b_list_train[[ 3 ]] <- predict(mu_y_short, train[[ m_each[[i]]["ay"]  ]])
            b_list_valid[[ 3 ]] <- mu_y_short$preds[[ m_each[[i]]["ay"]  ]]
            
            
            mu2_long_short <- mlr3superlearner::mlr3superlearner(
                data = data.frame(train$data[, na.omit(varnames_2)], tmp_pseudo_y = b_list_train[[3]] ),
                target = "tmp_pseudo_y",
                library = learners,
                outcome_type = "continuous",
                folds = control$mlr3superlearner_folds,
                newdata = valid,
                group = NULL,
                discrete = TRUE
            )
            
            b_list_train[[2]] <- predict(mu2_long_short, train[[ m_each[[i]]["az"] ]])
            b_list_valid[[2]] <- mu2_long_short$preds[[ m_each[[i]]["az"] ]]
            
            mu_fit[[2]] <- mu2_long_short$preds$data
            
            mu2_short <- mlr3superlearner::mlr3superlearner(
                data = data.frame(train$data[, na.omit(varnames_2_short)], tmp_pseudo_y = b_list_train[[3]] ),
                target = "tmp_pseudo_y",
                library = learners,
                outcome_type = "continuous",
                folds = control$mlr3superlearner_folds,
                newdata = valid,
                group = NULL,
                discrete = TRUE
            )
            mu1_long_short <- mlr3superlearner::mlr3superlearner(
                data = data.frame(train$data[, na.omit(varnames_1)], 
                                  tmp_pseudo_y = predict(mu2_short, 
                                                         train[[ m_each[[i]]["az"] ]]) ),
                target = "tmp_pseudo_y",
                library = learners,
                outcome_type = "continuous",
                folds = control$mlr3superlearner_folds,
                newdata = valid,
                group = NULL,
                discrete = TRUE
            )
            
            if ( length(m_each[[i]]) == 3 ) {
                b_list_train[[1]] <- predict(mu1_long_short, train[[ m_each[[i]]["am"] ]])
                b_list_valid[[1]] <- mu1_long_short$preds[[ m_each[[i]]["am"] ]]
                
            }
            if ( length(m_each[[i]]) == 4 ) {
                b_list_train[[1]] <- predict(mu1_long_short, train[[ m_each[[i]]["am.p"] ]]) - predict(mu1_long_short, train[[ m_each[[i]]["am.m"] ]])
                b_list_valid[[1]] <- mu1_long_short$preds[[ m_each[[i]]["am.p"] ]] - mu1_long_short$preds[[ m_each[[i]]["am.m"] ]]
            }
            
            mu_fit[[1]] <- mu1_long_short$preds$data
            
            mu_list[[i]] <- data.frame(do.call(cbind, mu_fit), do.call(cbind, b_list_valid))
            
            
            # if (m_each1[["az"]] == m_each1[["am"]]) {
            #   b_list_train <- b_list_valid <- mu_fit <- vector("list", length = 2) 
            #   names(b_list_valid) <- c("b_am_indepM", "b_ay_indepM")
            #   names(mu_fit) <- glue("mu_fit_am_indepM", "mu_fit_ay_indepM")
            #   
            #   b_list_train[[ "b_ay_indepM" ]] <- predict(mu_y, train[[ m_each[[i]]["ay"]  ]])
            #   b_list_valid[[ "b_ay_indepM" ]] <- mu_y$preds[[ m_each[[i]]["ay"]  ]]
            #   mu_fit[[ "mu_fit_ay_indepM" ]] <- mu_y$preds$data
            #   
            #   
            #   mu2 <- mlr3superlearner::mlr3superlearner(
            #     data = data.frame(train$data[, na.omit(varnames_1)], tmp_pseudo_y = b_list_train[[ "b_ay_indepM" ]] ),
            #     target = "tmp_pseudo_y",
            #     library = learners,
            #     outcome_type = "continuous",
            #     folds = control$mlr3superlearner_folds,
            #     newdata = valid,
            #     group = NULL,
            #     discrete = TRUE
            #   )
            #   
            #   b_list_train[[ "b_am_indepM" ]] <- predict(mu2, train[[ m_each[[i]]["am"] ]])
            #   b_list_valid[[ "b_am_indepM" ]] <- mu2$preds[[ m_each[[i]]["am"] ]]
            #   mu_fit[[ "mu_fit_am_indepM" ]] <- mu2$preds$data
            #   
            #   mu_list[[i]] <- data.frame(do.call(cbind, mu_fit), do.call(cbind, b_list_valid))
            # }
        }
        
        if (m_each1[["az"]] == "") {
            b_list_train <- b_list_valid <- mu_fit <- vector("list", length = 2) 
            names(b_list_valid) <- c("b_am_jo", "b_ay_jo")
            names(mu_fit) <- glue("mu_fit_am_jo", "mu_fit_ay_jo")
            
            mu_fit[[ "mu_fit_ay_jo" ]] <- mu_y_long$preds$data
            
            b_list_train[[ "b_ay_jo" ]] <- predict(mu_y_short, train[[ m_each[[i]]["ay"]  ]])
            b_list_valid[[ "b_ay_jo" ]] <- mu_y_short$preds[[ m_each[[i]]["ay"]  ]]
            
            
            mu2_long_short <- mlr3superlearner::mlr3superlearner(
                data = data.frame(train$data[, na.omit(varnames_1)], tmp_pseudo_y = b_list_train[[ "b_ay_jo" ]] ),
                target = "tmp_pseudo_y",
                library = learners,
                outcome_type = "continuous",
                folds = control$mlr3superlearner_folds,
                newdata = valid,
                group = NULL,
                discrete = TRUE
            )
            
            b_list_train[[ "b_am_jo" ]] <- predict(mu2_long_short, train[[ m_each[[i]]["am"] ]])
            b_list_valid[[ "b_am_jo" ]] <- mu2_long_short$preds[[ m_each[[i]]["am"] ]]
            
            mu_fit[[ "mu_fit_am_jo" ]] <- mu2_long_short$preds$data
            
            mu_list[[i]] <- data.frame(do.call(cbind, mu_fit), do.call(cbind, b_list_valid))
            
        }  
        
    }
    names(mu_list) <- map_chr(m_each, ~paste(glue( '{names(.)}{gsub("data_", "", .)}' ), collapse = "_") )
    mu_list
}




