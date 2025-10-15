
fit_alpha_each <- function(med, a_c, a_mc, m_ac, m_each, folds, learners_weights, control, covariates_cl = FALSE, benchmark_covariates = NULL) {

    alpha_each_list <- lapply(m_each, \(m_each1) {
        am <- as.numeric(substr(m_each1[["am"]], 6, 6))
        ay <- as.numeric(substr(m_each1[["ay"]], 6, 6))
            alpha_each <- data.frame(
                alpha1 = alpha.ac(A=med@data[[med@vars@A]], am, a_c),
                alpha2 = alpha.amc(A=med@data[[med@vars@A]], ay, am, a_c, a_mc, m_ac),
                alpha1_target = alpha.ac(A=am, am, a_c),
                alpha2_target = alpha.amc(A=ay, ay, am, a_c, a_mc, m_ac))
        }
    )

    alpha_each_list
}



a.c <- function(med, folds, learners, control, covariates_cl = FALSE) {
    # data_in <- med@data
    a_c <- matrix(nrow = nrow(med@data), ncol = 2)
    colnames(a_c) <- c("a(0|c)", "a(1|c)")

    v <- 1
    for (v in seq_along(folds)) {
        train <- training(med, folds, v)
        valid <- validation(med, folds, v)

        vars <- med@vars
        varnames_1 <- c(vars@A, vars@C)

        if (covariates_cl) {
            varnames_1 <- na.omit(c(vars@A, glue("{c(vars@A)}_clmean"), glue("{c(vars@C)}_cwc"), vars@Wj))
        }

        fit <- mlr3superlearner::mlr3superlearner(
            data = train$data[, varnames_1],
            target = med@vars@A,
            library = learners,
            outcome_type = "binomial",
            folds = control$mlr3superlearner_folds,
            newdata = valid,
            group = NULL,
            discrete = TRUE
        )

        a_c[folds[[v]]$validation_set, "a(0|c)"] <- 1 - fit$preds$data
        a_c[folds[[v]]$validation_set, "a(1|c)"] <- fit$preds$data
    }
    a_c
}


a.mc <- function(med, folds, learners, control, covariates_cl = FALSE) {
  
    a_mc <- matrix(nrow = nrow(med@data), ncol = 2)
    colnames(a_mc) <- c("a(0|m,c)", "a(1|m,c)")

    vars <- med@vars
    varnames_y <- na.omit(c(vars@A, vars@C, vars@M))

    if (covariates_cl) {
        varnames_y <- na.omit(c(c(vars@A), glue("{c(vars@A, vars@M)}_clmean"), glue("{c(vars@C, vars@M)}_cwc"), vars@Wj))
    }

    v <- 1
    for (v in seq_along(folds)) {
        train <- training(med, folds, v)
        valid <- validation(med, folds, v)
        fit <- mlr3superlearner::mlr3superlearner(
            data = train$data[, varnames_y],
            target = med@vars@A,
            library = learners,
            outcome_type = "binomial",
            folds = control$mlr3superlearner_folds,
            newdata = valid,
            group = NULL,
            discrete = TRUE
        )

        a_mc[folds[[v]]$validation_set, "a(0|m,c)"] <- 1 - fit$preds$data
        a_mc[folds[[v]]$validation_set, "a(1|m,c)"] <- fit$preds$data
    }
    a_mc
}


m.ac <- function(med, folds, learners, control, covariates_cl = FALSE) {
    # data_in <- med@data
    m_ac <- matrix(nrow = nrow(med@data), ncol = 2)
    colnames(m_ac) <- c("m(m|1,c)", "m(m|0,c)")

    vars <- med@vars
    varnames_y <- na.omit(c(vars@A, vars@C, vars@M))

    if (covariates_cl) {
        varnames_y <- na.omit(c(c(vars@A), glue("{c(vars@A, vars@M)}_clmean"), glue("{c(vars@C, vars@M)}_cwc"), vars@Wj))
    }

    v <- 1
    for (v in seq_along(folds)) {
        train <- training(med, folds, v)
        valid <- validation(med, folds, v)

        fit <- mlr3superlearner::mlr3superlearner(
            data = train$data[, varnames_y],
            target = med@vars@M,
            library = learners,
            outcome_type = "binomial",
            folds = control$mlr3superlearner_folds,
            newdata = valid,
            group = NULL,
            discrete = TRUE
        )

        m_ac[folds[[v]]$validation_set, "m(m|1,c)"] <- (valid$data[[med@vars@M]]==1)*fit$preds$data_1 + (valid$data[[med@vars@M]]==0)*(1 - fit$preds$data_1)

        m_ac[folds[[v]]$validation_set, "m(m|0,c)"] <- (valid$data[[med@vars@M]]==1)*fit$preds$data_0 + (valid$data[[med@vars@M]]==0)*(1 - fit$preds$data_0)
    }
    m_ac
}




# weights 

alpha.ac <- function(A, am, a_c) {
    ipw_a <- 1*(A==am) / pmax(a_c[, glue("a({am}|c)")], 1e-2) #tol = 1e-2
    # normalize(ipw_a) # no weight normalization to compute nu2s
    ipw_a
}

alpha.amc <- function(A, ay, am, a_c, a_mc, m_ac) {

    if (!is.null(m_ac)) {
        h_am <- 1*(A==ay) * m_ac[, glue("m(m|{am},c)")] / pmax(a_c[, glue("a({ay}|c)")] * m_ac[, glue("m(m|{ay},c)")], 1e-2)
    } else {
        h_am <- ( 1*(A==ay) * a_mc[, glue("a({am}|m,c)")] ) / 
            pmax(a_c[, glue("a({am}|c)")] * a_mc[, glue("a({ay}|m,c)")], 1e-2)
    }


    # normalize(h_am) # no weight normalization here to compute nu2s
    h_am
}






# NOT USED

# r.zmac <- function(med, permuted, #varnames, cluster_opt = "cwc", bounded = TRUE,
#                     folds, learners, control, covariates_cl = FALSE) {
#
#     r_zmac <- matrix(NA, nrow = nrow(med@data), ncol = 1)
#     colnames(r_zmac) <- "r_zmac"
#     # two_sample[["tmp_id"]] <- rep(1:nrow(med@data), times = 2)
#
#     v <- 1
#     for (v in seq_along(folds)) {
#         # train <- origami::training(data_in, folds[[v]])
#         # valid <- origami::validation(data_in, folds[[v]])
#         train <- training(med, folds, v)
#         valid <- validation(med, folds, v)
#
#         pseudo <- train$data
#         # # P <- linear_permutation(pseudo[, c(med@vars@A, med@vars@C)])
#         for (j in permuted) {
#             pseudo[, j] <- sample_M(train$data[, j, drop = TRUE])
#             # pseudo[, j] <- as.vector(P %*% as.matrix(pseudo[, j]))
#         }
#
#         # pseudo <- train$data_indepM
#
#         stacked_train <- rbind(train$data, pseudo)
#         stacked_train[["Lambda"]] <- rep(c(0, 1), each = nrow(train$data))
#
#         varnames_azmc <- c("Lambda", med@vars@A, med@vars@M, med@vars@C)
#         if (covariates_cl) {
#             varnames_azmc <- c("Lambda", med@vars@A, med@vars@M, med@vars@C)
#         }
#
#         fit_azmc <- mlr3superlearner::mlr3superlearner(
#             data = stacked_train[, varnames_azmc],
#             target = "Lambda",
#             library = learners,
#             outcome_type = "binomial",
#             folds = control$mlr3superlearner_folds,
#             newdata = valid,
#             group = NULL,
#             discrete = TRUE
#         )
#         # r_zmac[folds[[v]]$validation_set, "r_zmac"] <- fit_azmc$preds$data / pmax((1 - fit_azmc$preds$data), 0.01)
#
#         varnames_amc <- c("Lambda", med@vars@A, permuted, med@vars@C)
#         fit_amc <- mlr3superlearner::mlr3superlearner(
#             data = stacked_train[, varnames_amc],
#             target = "Lambda",
#             library = learners,
#             outcome_type = "binomial",
#             folds = control$mlr3superlearner_folds,
#             newdata = valid,
#             group = NULL,
#             discrete = TRUE
#         )
#         r_zmac[folds[[v]]$validation_set, "r_zmac"] <- fit_azmc$preds$data * (1 - fit_amc$preds$data) / pmax((1 - fit_azmc$preds$data) * fit_amc$preds$data, 0.001)
#
#
#     }
#     r_zmac
# }

# alpha.azmc <- function(A, ay, am, az, a_c, a_mc, a_zc, r_zmac) {
#
#     iorw_m <- a_mc[, glue("a({am}|m,c)")] * a_c[, glue("a({ay}|c)")] / pmax(a_mc[, glue("a({ay}|m,c)")] * a_c[, glue("a({am}|c)")], 1e-2)
#     iorw_z <- a_zc[, glue("a({az}|m,c)")] * a_c[, glue("a({ay}|c)")] / pmax(a_zc[, glue("a({ay}|m,c)")] * a_c[, glue("a({az}|c)")], 1e-2)
#     ipw_a <- 1 * (A == ay) / pmax(a_c[, glue("a({ay}|c)")], 1e-2)
#
#     h_azm <- ipw_a * iorw_m * iorw_z * r_zmac
#     # h_azm <- r_zmac * ( 1*(A==ay) * a_mc[, glue("a({am}|m,c)")] * a_zc[, glue("a({az}|m,c)")] )  /
#     #     pmax(a_c[, glue("a({ay}|c)")] * a_mc[, glue("a({ay}|m,c)")] * a_zc[, glue("a({ay}|m,c)")], 0.01)
#
#     # h_azm
#     normalize(h_azm)
# }
