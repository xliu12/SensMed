# --------------
estimate_h1.IIE <- function(med, a_mc, mus_IIE, folds, learners_regressions, control, covariates_cl = FALSE, benchmark_covar = NULL) {
    # add the pseudo-outcome to the data
    Y <- med@data[[med@vars@Y]]
    A <- med@data[[med@vars@A]]
    med@data[["wYresid"]] <- 1*(A==1)/bound(a_mc[, "a(1|m,c)"]) * (Y - mus_IIE$mu_fit_ay_jo)

    h1_folds <- vector("list", control$crossfit_folds)
    i <- 1
    for (i in seq_along(h1_folds)) {
        train <- training(med, folds, i)
        valid <- validation(med, folds, i)
        h1_folds[[i]] <- h1.IIE(train, valid, med@vars, learners_regressions, control, covariates_cl, benchmark_covar)
    }

    h1_IIE <- combine_folds_theta(h1_folds, folds)
    names(h1_IIE) <- NULL
    h1_IIE
}

estimate_h1.IDE <- function(med, a_mc, mus_IDE, folds, learners_regressions, control, covariates_cl = FALSE, benchmark_covar = NULL) {
    # add the pseudo-outcome to the data
    Y <- med@data[[med@vars@Y]]
    A <- med@data[[med@vars@A]]
    med@data[["wYresid"]] <- (1*(A==1)/bound(a_mc[, "a(1|m,c)"]) - 1*(A==0)/bound(a_mc[, "a(0|m,c)"])) * (Y - mus_IDE$mu_fit_ay_jo)

    h1_folds <- vector("list", control$crossfit_folds)
    i <- 1
    for (i in seq_along(h1_folds)) {
        train <- training(med, folds, i)
        valid <- validation(med, folds, i)
        h1_folds[[i]] <- h1.IDE(train, valid, med@vars,  learners_regressions, control, covariates_cl, benchmark_covar)
    }
    h1_IDE <- combine_folds_theta(h1_folds, folds)
    names(h1_IDE) <- NULL
    h1_IDE
}

h1.IIE <- function(train, valid, vars, learners, control, covariates_cl = FALSE, benchmark_covar = NULL) {

    if (length(benchmark_covar) > 0) {
        shortC <- setdiff(vars@C, benchmark_covar)
    } else {
        shortC <- vars@C
    }
    # for IIE b2 = theta2(1,m,x)

    h1 <- mlr3superlearner::mlr3superlearner(
        data = train$data[, na.omit(c(vars@A, shortC, "wYresid"))],
        # data = train$data[train$data[[vars@A]]==1, na.omit(c(shortC, "wYresid"))],
        target = "wYresid",
        library = learners,
        outcome_type = "continuous",
        folds = control$mlr3superlearner_folds,
        newdata = valid,
        group = NULL,
        discrete = TRUE
    )

    h1_fit <- h1$preds[[ "data" ]]
    # h1_fit[valid$data[[vars@A]] == 0] <- 0

    h1_df <- data.frame(h1_fit = h1_fit)
}




h1.IDE <- function(train, valid, vars, learners, control, covariates_cl = FALSE, benchmark_covar = NULL) {
    # vars <- med@vars

    if (length(benchmark_covar) > 0) {
        shortC <- setdiff(vars@C, benchmark_covar)
    } else {
        shortC <- vars@C
    }

    h1 <- mlr3superlearner::mlr3superlearner(
        data = train$data[, na.omit(c(vars@A, shortC, "wYresid"))],
        # data = train$data[train$data[[vars@A]]==1, na.omit(c(shortC, "wYresid"))],
        target = "wYresid",
        library = learners,
        outcome_type = "continuous",
        folds = control$mlr3superlearner_folds,
        newdata = valid,
        group = NULL,
        discrete = TRUE
    )

    h1_fit <- h1$preds[[ "data" ]]
    # h1_fit[valid$data[[vars@A]] == 0] <- 0

    h1_df <- data.frame(h1_fit = h1_fit)
}


# Old -----------
# estimate_h1.IIE <- function(med, a_c, a_mc, mus_IIE, folds, learners_regressions, control, covariates_cl = FALSE, benchmark_covar = NULL) {
#     # add the pseudo-outcome to the data
#     Y <- med@data[[med@vars@Y]]
#     A <- med@data[[med@vars@A]]
#     med@data[["wYresid"]] <- 1*(A==1)*a_c[, "a(1|c)"]/bound(a_mc[, "a(1|m,c)"]) * (Y - mus_IIE$mu_fit_ay_jo)
#
#     h1_folds <- vector("list", control$crossfit_folds)
#     i <- 1
#     for (i in seq_along(h1_folds)) {
#         train <- training(med, folds, i)
#         valid <- validation(med, folds, i)
#         h1_folds[[i]] <- h1.IIE(train, valid, med@vars, learners_regressions, control, covariates_cl, benchmark_covar)
#     }
#
#     h1_IIE <- combine_folds_theta(h1_folds, folds)
#     names(h1_IIE) <- NULL
#     h1_IIE
# }
#
# estimate_h1.IDE <- function(med, a_c, a_mc, mus_IDE, folds, learners_regressions, control, covariates_cl = FALSE, benchmark_covar = NULL) {
#     # add the pseudo-outcome to the data
#     Y <- med@data[[med@vars@Y]]
#     A <- med@data[[med@vars@A]]
#     med@data[["wYresid"]] <- (1*(A==1)*a_c[, "a(1|c)"]/bound(a_mc[, "a(1|m,c)"]) - 1*(A==0)*a_c[, "a(0|c)"]/bound(a_mc[, "a(0|m,c)"])) * (Y - mus_IDE$mu_fit_ay_jo)
#
#     h1_folds <- vector("list", control$crossfit_folds)
#     i <- 1
#     for (i in seq_along(h1_folds)) {
#         train <- training(med, folds, i)
#         valid <- validation(med, folds, i)
#         h1_folds[[i]] <- h1.IDE(train, valid, med@vars,  learners_regressions, control, covariates_cl, benchmark_covar)
#     }
#     h1_IDE <- combine_folds_theta(h1_folds, folds)
#     names(h1_IDE) <- NULL
#     h1_IDE
# }
#
# h1.IIE <- function(train, valid, vars, learners, control, covariates_cl = FALSE, benchmark_covar = NULL) {
#
#     if (length(benchmark_covar) > 0) {
#         shortC <- setdiff(vars@C, benchmark_covar)
#     } else {
#         shortC <- vars@C
#     }
#     # for IIE b2 = theta2(1,m,x)
#
#     h1 <- mlr3superlearner::mlr3superlearner(
#         data = train$data[, na.omit(c(vars@A, shortC, "wYresid"))],
#         # data = train$data[train$data[[vars@A]]==1, na.omit(c(shortC, "wYresid"))],
#         target = "wYresid",
#         library = learners,
#         outcome_type = "continuous",
#         folds = control$mlr3superlearner_folds,
#         newdata = valid,
#         group = NULL,
#         discrete = TRUE
#     )
#
#     h1_fit <- h1$preds[[ "data" ]]
#     # h1_fit[valid$data[[vars@A]] == 0] <- 0
#
#     h1_df <- data.frame(h1_fit = h1_fit)
# }
#
#
#
#
# h1.IDE <- function(train, valid, vars, learners, control, covariates_cl = FALSE, benchmark_covar = NULL) {
#     # vars <- med@vars
#
#     if (length(benchmark_covar) > 0) {
#         shortC <- setdiff(vars@C, benchmark_covar)
#     } else {
#         shortC <- vars@C
#     }
#
#     h1 <- mlr3superlearner::mlr3superlearner(
#         data = train$data[, na.omit(c(vars@A, shortC, "wYresid"))],
#         # data = train$data[train$data[[vars@A]]==1, na.omit(c(shortC, "wYresid"))],
#         target = "wYresid",
#         library = learners,
#         outcome_type = "continuous",
#         folds = control$mlr3superlearner_folds,
#         newdata = valid,
#         group = NULL,
#         discrete = TRUE
#     )
#
#     h1_fit <- h1$preds[[ "data" ]]
#     # h1_fit[valid$data[[vars@A]] == 0] <- 0
#
#     h1_df <- data.frame(h1_fit = h1_fit)
# }





# estimate_h1.IIE <- function(med, mus_IIE, folds, learners_regressions, control, covariates_cl = FALSE, benchmark_covar = NULL) {
#     # add the pseudo-outcome to the data
#     Y <- med@data[[med@vars@Y]]
#     # A <- med@data[[med@vars@A]]
#     med@data[["theta1"]] <- mus_IIE$mu_fit_am_jo
#
#     h1_folds <- vector("list", control$crossfit_folds)
#     i <- 1
#     for (i in seq_along(h1_folds)) {
#         train <- training(med, folds, i)
#         valid <- validation(med, folds, i)
#         h1_folds[[i]] <- h1.IIE(train, valid, med@vars, learners_regressions, control, covariates_cl, benchmark_covar)
#     }
#
#     h1_IIE <- combine_folds_theta(h1_folds, folds)
#     names(h1_IIE) <- NULL
#     h1_IIE
# }
#
# estimate_h1.IDE <- function(med, mus_IDE, folds, learners_regressions, control, covariates_cl = FALSE, benchmark_covar = NULL) {
#     # add the pseudo-outcome to the data
#     Y <- med@data[[med@vars@Y]]
#     # A <- med@data[[med@vars@A]]
#     med@data[["theta1"]] <- mus_IDE$mu_fit_am_jo
#
#     h1_folds <- vector("list", control$crossfit_folds)
#     i <- 1
#     for (i in seq_along(h1_folds)) {
#         train <- training(med, folds, i)
#         valid <- validation(med, folds, i)
#         h1_folds[[i]] <- h1.IDE(train, valid, med@vars,  learners_regressions, control, covariates_cl, benchmark_covar)
#     }
#     h1_IDE <- combine_folds_theta(h1_folds, folds)
#     names(h1_IDE) <- NULL
#     h1_IDE
# }
#
# h1.IIE <- function(train, valid, vars, learners, control, covariates_cl = FALSE, benchmark_covar = NULL) {
#
#     if (length(benchmark_covar) > 0) {
#         shortC <- setdiff(vars@C, benchmark_covar)
#     } else {
#         shortC <- vars@C
#     }
#     # for IIE b2 = theta2(1,m,x)
#
#     h1 <- mlr3superlearner::mlr3superlearner(
#         data = train$data[, na.omit(c(vars@M, shortC, "theta1"))],
#         # data = train$data[train$data[[vars@A]]==1, na.omit(c(shortC, "wYresid"))],
#         target = "theta1",
#         library = learners,
#         outcome_type = "continuous",
#         folds = control$mlr3superlearner_folds,
#         newdata = valid,
#         group = NULL,
#         discrete = TRUE
#     )
#
#     h1_fit <- h1$preds[[ "data" ]]
#     # h1_fit[valid$data[[vars@A]] == 0] <- 0
#
#     h1_df <- data.frame(h1_fit = h1_fit)
# }
#
#
#
#
# h1.IDE <- function(train, valid, vars, learners, control, covariates_cl = FALSE, benchmark_covar = NULL) {
#     # vars <- med@vars
#
#     if (length(benchmark_covar) > 0) {
#         shortC <- setdiff(vars@C, benchmark_covar)
#     } else {
#         shortC <- vars@C
#     }
#
#     h1 <- mlr3superlearner::mlr3superlearner(
#         data = train$data[, na.omit(c(vars@M, shortC, "theta1"))],
#         # data = train$data[train$data[[vars@A]]==1, na.omit(c(shortC, "wYresid"))],
#         target = "theta1",
#         library = learners,
#         outcome_type = "continuous",
#         folds = control$mlr3superlearner_folds,
#         newdata = valid,
#         group = NULL,
#         discrete = TRUE
#     )
#
#     h1_fit <- h1$preds[[ "data" ]]
#     # h1_fit[valid$data[[vars@A]] == 0] <- 0
#
#     h1_df <- data.frame(h1_fit = h1_fit)
# }



