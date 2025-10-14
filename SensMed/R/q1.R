
estimate_q1.IIE_Mjo <- function(med, alphas_IIE, folds, learners_regressions, control, covariates_cl = FALSE, benchmark_covar = NULL) {

    med@data$b2_alpha2 <- alphas_IIE[["alpha2_target"]] # add to data as pseudo outcome, so will appear in the training data

    q1_folds <- vector("list", control$crossfit_folds)
    i <- 1
    # cli::cli_progress_step("Fitting outcome regressions... {i}/{control$crossfit_folds} folds")
    for (i in seq_along(q1_folds)) {
        train <- training(med, folds, i)
        valid <- validation(med, folds, i)
        q1_folds[[i]] <- q1.IIE_Mjo(train, valid, med@vars,  learners_regressions, control, covariates_cl, benchmark_covar)
        # cli::cli_progress_update()
    }
    # cli::cli_progress_done()
    q1_allfolds <- combine_folds_theta(q1_folds, folds)
    q1_allfolds
}

estimate_q1.IDE <- function(med, alphas_IDE, folds, learners_regressions, control, covariates_cl = FALSE, benchmark_covar = NULL) {

   med@data$b2_alpha2 <- alphas_IDE[["alpha2_target"]] # add to data as pseudo outcome, so will appear in the training data

  q1_folds <- vector("list", control$crossfit_folds)
  i <- 1
  # cli::cli_progress_step("Fitting outcome regressions... {i}/{control$crossfit_folds} folds")
  for (i in seq_along(q1_folds)) {
    train <- training(med, folds, i)
    valid <- validation(med, folds, i)
    q1_folds[[i]] <- q1.IDE(train, valid, med@vars,  learners_regressions, control, covariates_cl, benchmark_covar)
    # cli::cli_progress_update()
  }
  # cli::cli_progress_done()
  q1_allfolds <- combine_folds_theta(q1_folds, folds)
  q1_allfolds
}

q1.IIE_Mjo <- function(train, valid, vars, learners, control, covariates_cl = FALSE, benchmark_covar = NULL) {
  # vars <- med@vars
  # b2_alpha2 <- alphas[[2]][["alpha2_target"]] - alphas[[1]][["alpha2_target"]]
    # b2_alpha2_train <- alphas_IIE[["alpha2_target"]] # have added as pseudo-outcome in the estimate. function above

  if (length(benchmark_covar) > 0) {
    shortC <- setdiff(vars@C, benchmark_covar)
  } else {
    shortC <- vars@C
  }

  q1 <- mlr3superlearner::mlr3superlearner(
    data = train$data[, na.omit(c(vars@A, shortC, "b2_alpha2"))],
    target = "b2_alpha2",
    library = learners,
    outcome_type = "continuous",
    folds = control$mlr3superlearner_folds,
    newdata = valid, # validation data
    group = NULL,
    discrete = TRUE
  )

  m_each <- list(IIE_Mjo_m = c(az = "", am = "data_0", ay = "data_1"),
                 IIE_Mjo_p = c(az = "", am = "data_1", ay = "data_1"))

  b1_q1 <- q1$preds[[ m_each[[2]]["am"] ]] - q1$preds[[ m_each[[1]]["am"] ]]
  q1_fit <- q1$preds[[ "data" ]]
  b2_alpha2 <- valid$data[["b2_alpha2"]]

  q1_list <- data.frame(b2_alpha2 = b2_alpha2, q1_fit = q1_fit, b1_q1 = b1_q1)

  # names(q1_list) <- map_chr(~paste(glue( '{names(.)}{gsub("data_", "", .)}' ), collapse = "_") )
  q1_list
}




q1.IDE <- function(train, valid, vars, learners, control, covariates_cl = FALSE, benchmark_covar = NULL) {
  # vars <- med@vars

  # b2_alpha2 <- alphas_IDE[["alpha2_target"]] # have added as pseudo-outcome in the estimate. function above

  if (length(benchmark_covar) > 0) {
    shortC <- setdiff(vars@C, benchmark_covar)
  } else {
    shortC <- vars@C
  }

  q1 <- mlr3superlearner::mlr3superlearner(
    data = train$data[, na.omit(c(vars@A, shortC, "b2_alpha2"))],
    target = "b2_alpha2",
    library = learners,
    outcome_type = "continuous",
    folds = control$mlr3superlearner_folds,
    newdata = valid,
    group = NULL,
    discrete = TRUE
  )

  m_each <- list(IDE_m = c(az = "", am = "data_0", ay = "data_0"),
                 IDE_p = c(az = "", am = "data_0", ay = "data_1"))

  # for IDE, am is the same
  b1_q1 <- q1$preds[[ m_each[[2]]["am"] ]]  # m_each[[2]]["am"] == m_each[[1]]["am"]
  q1_fit <- q1$preds[[ "data" ]]
  b2_alpha2 <- valid$data[["b2_alpha2"]] # the pseudo-outcome in the validation data

  q1_list <- data.frame(b2_alpha2 = b2_alpha2, q1_fit = q1_fit, b1_q1 = b1_q1)

  # names(q1_list) <- map_chr(~paste(glue( '{names(.)}{gsub("data_", "", .)}' ), collapse = "_") )
  q1_list
}



