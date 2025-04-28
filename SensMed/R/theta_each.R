
estimate_theta <- function(med, m_each, folds, learners_regressions, control, which.M) {
    b_theta_folds <- vector("list", control$crossfit_folds)
    i <- 1
    cli::cli_progress_step("Fitting outcome regressions... {i}/{control$crossfit_folds} folds")
    for (i in seq_along(b_theta_folds)) {
        train <- training(med, folds, i)
        valid <- validation(med, folds, i)
        b_theta_folds[[i]] <- theta_indepM(train, valid, med@vars, m_each,  learners_regressions, control, which.M)
        cli::cli_progress_update()
    }
    cli::cli_progress_done()
    b_theta_allfolds <- combine_folds_theta(b_theta_folds, folds)
    b_theta_allfolds
}

theta_indepM <- function(train, valid, vars, m_each, learners, control, which.M) {
    
    y_type <- ifelse(if.binary(train$data[[vars@Y]]), "binomial", "continuous")
    theta_y <- mlr3superlearner::mlr3superlearner(
        data = train$data[, na.omit(c(vars@A, vars@C, vars@M, vars@Y))],
        target = vars@Y,
        library = learners,
        outcome_type = y_type,
        folds = control$mlr3superlearner_folds,
        newdata = valid,
        group = NULL
    )
    
    b_thetas <- vector("list", length = length(m_each))
    i <- 1
    for (i in 1:length(m_each)) {
        # m_each[[i]][L+1]
        b_list_train <- b_list_valid <- theta_fit <- vector("list", length = 3) 
        names(b_list_valid) <- glue("b{1:length(b_list_valid)}")
        names(theta_fit) <- glue("theta_fit{1:length(theta_fit)}")
        
        theta_m <- vector("list", length = 2) 
        names(theta_m) <- glue("theta{1:2}") 
        
        b_list_train[[ 3 ]] <- predict(theta_y, train[[ m_each[[i]]["ay"]  ]])
        b_list_valid[[ 3 ]] <- theta_y$preds[[ m_each[[i]]["ay"]  ]]
        theta_fit[[3]] <- theta_y$preds$data
        
        theta_m[[2]] <- mlr3superlearner::mlr3superlearner(
            data = data.frame(train$data[, na.omit(c(vars@A, vars@C, vars@M[-which.M])) ], tmp_pseudo_y = b_list_train[[3]] ),
            target = "tmp_pseudo_y",
            library = learners,
            outcome_type = "continuous",
            folds = control$mlr3superlearner_folds,
            newdata = valid,
            group = NULL
        )
        
        b_list_train[[2]] <- predict(theta_m[[2]], train[[ m_each[[i]]["az"] ]])
        b_list_valid[[2]] <- theta_m[[2]]$preds[[ m_each[[i]]["az"] ]]
        theta_fit[[2]] <- theta_m[[2]]$preds$data
        
        theta_m[[1]] <- mlr3superlearner::mlr3superlearner(
            data = data.frame(train$data[, na.omit(c(vars@A, vars@C)) ], tmp_pseudo_y = b_list_train[[2]] ),
            target = "tmp_pseudo_y",
            library = learners,
            outcome_type = "continuous",
            folds = control$mlr3superlearner_folds,
            newdata = valid,
            group = NULL
        )
        
        b_list_train[[1]] <- predict(theta_m[[1]], train[[ m_each[[i]]["am"] ]])
        b_list_valid[[1]] <- theta_m[[1]]$preds[[ m_each[[i]]["am"] ]]
        theta_fit[[1]] <- theta_m[[1]]$preds$data
        
        b_thetas[[i]] <- data.frame(do.call(cbind, theta_fit), do.call(cbind, b_list_valid))
         
    }
    names(b_thetas) <- map_chr(m_each, ~paste(glue( '{names(.)}{gsub("data_", "", .)}' ), collapse = "_") )
    b_thetas
}



