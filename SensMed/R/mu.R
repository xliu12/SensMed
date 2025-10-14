
estimate_mu_cl <- function(med, m_each, folds, learners_regressions, control, covariates_cl = FALSE, benchmark_covar = NULL) {
    mu_folds <- vector("list", control$crossfit_folds)
    i <- 1
    # cli::cli_progress_step("Fitting outcome regressions... {i}/{control$crossfit_folds} folds")
    for (i in seq_along(mu_folds)) {
        train <- training(med, folds, i)
        valid <- validation(med, folds, i)
        mu_folds[[i]] <- mu_eachM(train, valid, med@vars, m_each,  learners_regressions, control, covariates_cl, benchmark_covar)
        # cli::cli_progress_update()
    }
    # cli::cli_progress_done()
    mu_allfolds <- combine_folds_theta(mu_folds, folds)
    mu_allfolds
}

mu_eachM <- function(train, valid, vars, m_each, learners, control, covariates_cl = FALSE, benchmark_covar = NULL) {
  # vars <- med@vars
  if (length(benchmark_covar) > 0) {
    shortC <- setdiff(vars@C, benchmark_covar)
  } else {
    shortC <- vars@C
  }

  varnames_1 <- c(vars@A, shortC)
  varnames_y <- c(vars@A, shortC, vars@M)

  # if (covariates_cl) {
  #   varnames_1 <- na.omit(c(glue("{c(vars@A, shortC)}_cwc"), glue("{c(vars@Y, vars@M)}_clmean"), vars@Wj))
  #   varnames_y <- na.omit(c(glue("{c(vars@A, shortC, vars@M)}_cwc"), glue("{c(vars@Y)}_clmean"), vars@Wj))
  # }

    y_type <- ifelse(if.binary(train$data[[vars@Y]]), "binomial", "continuous")

    mu_y <- mlr3superlearner::mlr3superlearner(
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

      b_list_train <- b_list_valid <- mu_fit <- vector("list", length = 2)
      names(b_list_valid) <- c("b_am_jo", "b_ay_jo")
      names(mu_fit) <- c("mu_fit_am_jo", "mu_fit_ay_jo")

      b_list_train[[ "b_ay_jo" ]] <- predict(mu_y, train[[ m_each[[i]]["ay"]  ]])
      b_list_valid[[ "b_ay_jo" ]] <- mu_y$preds[[ m_each[[i]]["ay"]  ]]
      mu_fit[[ "mu_fit_ay_jo" ]] <- mu_y$preds$data

      set.seed(12)
      mu2 <- mlr3superlearner::mlr3superlearner(
        data = data.frame(train$data[, na.omit(varnames_1)], tmp_pseudo_y = b_list_train[[ "b_ay_jo" ]] ),
        target = "tmp_pseudo_y",
        library = learners,
        outcome_type = "continuous",
        folds = control$mlr3superlearner_folds,
        newdata = valid,
        group = NULL,
        discrete = TRUE
      )

      b_list_train[[ "b_am_jo" ]] <- predict(mu2, train[[ m_each[[i]]["am"] ]])
      b_list_valid[[ "b_am_jo" ]] <- mu2$preds[[ m_each[[i]]["am"] ]]
      mu_fit[[ "mu_fit_am_jo" ]] <- mu2$preds$data

      mu_list[[i]] <- data.frame(do.call(cbind, mu_fit), do.call(cbind, b_list_valid))

    }
    # names(mu_list) <- map_chr(m_each, ~paste(glue( '{names(.)}{gsub("data_", "", .)}' ), collapse = "_") )
    names(mu_list) <- names(m_each)

    mu_list
}




