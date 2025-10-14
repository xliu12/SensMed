
estimate_mu.long.short <- function(benchmark_covariates,
                                   med, m_each, folds, learners_regressions, control, covariates_cl = FALSE) {
    mu_folds <- vector("list", control$crossfit_folds)
    i <- 1
    # cli::cli_progress_step("Fitting outcome regressions... {i}/{control$crossfit_folds} folds")
    for (i in seq_along(mu_folds)) {
        train <- training(med, folds, i)
        valid <- validation(med, folds, i)
        mu_folds[[i]] <- mu_eachM.long.short(benchmark_covariates,
                                  train, valid, med@vars, m_each,  learners_regressions, control, covariates_cl)
        # cli::cli_progress_update()
    }
    # cli::cli_progress_done()
    mu_allfolds <- combine_folds_theta(mu_folds, folds)
    mu_allfolds
}

mu_eachM.long.short <- function(benchmark_covariates,
                                train, valid, vars, m_each, learners, control, covariates_cl = FALSE) {

    varnames_1 <- c(vars@A, vars@C)
    varnames_1_short <- setdiff(varnames_1, benchmark_covariates)

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
        m_each1 <- m_each[[i]]
        b_list_train <- b_list_valid <- mu_fit <- vector("list", length = 2)
        names(b_list_valid) <- c("b_am_jo", "b_ay_jo")
        names(mu_fit) <- c("mu_fit_am_jo", "mu_fit_ay_jo")

        mu_fit[[ "mu_fit_ay_jo" ]] <- mu_y_long$preds$data

        # set the pseudo outcome to the short regression theta2s's prediction
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
    names(mu_list) <- names(m_each)
    mu_list
}




