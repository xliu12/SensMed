
cal.sigma1sq.IIE <- function(med, a_mc, mus_IIE, mus_Mjo) {
    # h1_IIE <- estimate_h1.IIE(med, a_mc,
    #                           mus_IIE, folds, learners_regressions, control, covariates_cl, benchmark_covar)

    Y <- med@data[[med@vars@Y]]
    A <- med@data[[med@vars@A]]
    w1_mc <- 1*(A==1)/pmax(a_mc[, "a(1|m,c)"], 1e-2)
    # \hat e(M,X) = \E[\hat\theta_1(A,X_\ss) \mid M,X_\ss] = \sum_{a=0,1}\hat\theta_1(a,X_\ss)\p(A=a\mid M,X_\ss)
    e_mc <- mus_Mjo$y01$b_am_jo * a_mc[, "a(0|m,c)"] + mus_Mjo$y11$b_am_jo * a_mc[, "a(1|m,c)"]

    eif_sigma1sq <- 2 * w1_mc * (Y - mus_IIE$mu_fit_ay_jo) * (mus_IIE$b_ay_jo - e_mc) + (mus_IIE$b_ay_jo - mus_IIE$mu_fit_am_jo)^2

    # mean(eif_sigma1sq)
    # mean((mus_IIE$b_ay_jo - mus_IIE$mu_fit_am_jo)^2)

    # b2 <- mus_IIE$b_ay_jo # mus[[2]]$b_ay_jo = mus[[1]]$b_ay_jo for IIE
    # theta1 <- mus_IIE$mu_fit_am_jo
    #
    # eif_sigma1sq <- 2*(wYresid - h1_IIE) * (b2 - theta1) + (b2 - theta1)^2

    sigma1sq <- mean(eif_sigma1sq)
    list(sigma1sq = sigma1sq, eif_sigma1sq = eif_sigma1sq)
}

cal.sigma1sq.IDE <- function(med, a_mc, mus_IDE, mus_Mjo) {
    # h1_IDE <- estimate_h1.IDE(med, a_mc,
    #                           mus_IDE, folds, learners_regressions, control, covariates_cl, benchmark_covar)

    Y <- med@data[[med@vars@Y]]
    A <- med@data[[med@vars@A]]
    w_mc_ay1 <- 1*(A==1)/pmax(a_mc[, "a(1|m,c)"], 1e-2)
    w_mc_ay0 <- 1*(A==0)/pmax(a_mc[, "a(0|m,c)"], 1e-2)
    # \hat e(M,X) = \E[\hat\theta_1(A,X_\ss) \mid M,X_\ss] = \sum_{a=0,1}\hat\theta_1(a,X_\ss)\p(A=a\mid M,X_\ss)
    e_mc <- (mus_Mjo$y01$b_am_jo - mus_Mjo$y00$b_am_jo) * a_mc[, "a(0|m,c)"] + (mus_Mjo$y11$b_am_jo - mus_Mjo$y10$b_am_jo) * a_mc[, "a(1|m,c)"]

    eif_sigma1sq <- 2 * (w_mc_ay1 * (Y - mus_IDE$mu_fit_ay_jo) - w_mc_ay0 * (Y - mus_IDE$mu_fit_ay_jo)) * (mus_IDE$b_ay_jo - e_mc) + (mus_IDE$b_ay_jo - mus_IDE$mu_fit_am_jo)^2

    mean(eif_sigma1sq)
    # mean((mus_IDE$b_ay_jo - mus_IDE$mu_fit_am_jo)^2)

    # wYresid <- (1*(A==1)/bound(a_mc[, "a(1|m,c)"], tol = 1e-2) - 1*(A==0)/bound(a_mc[, "a(0|m,c)"], tol = 1e-2)) * (Y - mus_IDE$mu_fit_ay_jo)

    # b2 <- mus_IDE$b_ay_jo # mus[[2]]$b_ay_jo = mus[[1]]$b_ay_jo for IDE
    # theta1 <- mus_IDE$mu_fit_am_jo
    #
    # eif_sigma1sq <- 2*(wYresid - h1_IDE) * (b2 - theta1) + (b2 - theta1)^2

    sigma1sq <- mean(eif_sigma1sq)
    list(sigma1sq = sigma1sq, eif_sigma1sq = eif_sigma1sq)

}

cal.sigma2s.IIE_Mjo <- function(Y, mus_IIE) {
    # for IIE, ay is the same for both potential outcomes, am differs
    resY <- R2_Y <- list("sigma1" = NA, "sigma2" = NA)

    resY[["sigma2"]] <- (Y - mus_IIE$mu_fit_ay_jo) # is the same as mus[[1]]$mu_fit_ay_jo
    R2_Y[["sigma2"]] <- max(1 - var(resY[["sigma2"]]) / var(Y), 0)

    b2 <- mus_IIE$b_ay_jo # mus[[2]]$b_ay_jo = mus[[1]]$b_ay_jo for IIE
    theta1 <- mus_IIE$mu_fit_am_jo
    resY[["sigma1"]] <- (b2 - theta1)
    R2_Y[["sigma1"]] <- max(1 - var(resY[["sigma1"]]) / var(b2), 0)

    # verify
    # mean((gen_data$dat$Y -gen_data$mus_IIE$mu_fit_ay_jo)^2 )
    # mean((gen_data$mus_IIE$b_ay_jo -gen_data$mus_IIE$mu_fit_am_jo)^2 )

    eif_sigma2s <- map(resY, \(x) x^2)
    sigma2s <- map(eif_sigma2s, \(x) mean(x))

    list(resY = resY, R2_Y = R2_Y, sigma2s = sigma2s, eif_sigma2s = eif_sigma2s)
}

cal.sigma2s.IDE <- function(Y, mus_IDE) {
    # for IDE, am is the same for both potential outcomes,  ay differs
    resY <- R2_Y <- list("sigma1" = NA, "sigma2" = NA)

    resY[["sigma2"]] <- (Y - mus_IDE$mu_fit_ay_jo) # is the same as mus[[1]]$mu_fit_ay_jo
    R2_Y[["sigma2"]] <- max(1 - var(resY[["sigma2"]]) / var(Y), 0)

    b2 <- mus_IDE$b_ay_jo
    theta1 <- mus_IDE$mu_fit_am_jo
    resY[["sigma1"]] <- (b2 - theta1)
    R2_Y[["sigma1"]] <- max(1 - var(resY[["sigma1"]]) / var(b2), 0)

    # verify
    # mean((gen_data$dat$Y -gen_data$mus_IDE$mu_fit_ay_jo)^2 )
    # mean((gen_data$mus_IDE$b_ay_jo -gen_data$mus_IDE$mu_fit_am_jo)^2 )

    eif_sigma2s <- map(resY, \(x) x^2)
    sigma2s <- map(eif_sigma2s, \(x) mean(x))

    list(resY = resY, R2_Y = R2_Y, sigma2s = sigma2s, eif_sigma2s = eif_sigma2s)
}


# NOT USED --------
# cal.sigma1sq.IIE <- function(med, a_mc, mus_IIE, folds, learners_regressions, control, covariates_cl = FALSE, benchmark_covar = NULL) {
#     h1_IIE <- estimate_h1.IIE(med, a_mc,
#                               mus_IIE, folds, learners_regressions, control, covariates_cl, benchmark_covar)
#
#     Y <- med@data[[med@vars@Y]]
#     A <- med@data[[med@vars@A]]
#     wYresid <- 1*(A==1)/bound(a_mc[, "a(1|m,c)"], tol = 1e-2) * (Y - mus_IIE$mu_fit_ay_jo)
#     # wYresid <- 1*(A==1)*a_c[, "a(1|c)"]/bound(a_mc[, "a(1|m,c)"]) * (Y - mus_IIE$mu_fit_ay_jo)
#
#     b2 <- mus_IIE$b_ay_jo # mus[[2]]$b_ay_jo = mus[[1]]$b_ay_jo for IIE
#     theta1 <- mus_IIE$mu_fit_am_jo
#
#     eif_sigma1sq <- 2*(wYresid - h1_IIE) * (b2 - theta1) + (b2 - theta1)^2
#     # wYresid2 <- wYresid*(Y - mus_IIE$mu_fit_ay_jo)
#     # eif_sigma1sq <- wYresid2 - h1_IIE^2 + (b2 - theta1)^2
#     # eif_sigma1sq <- 2*wYresid*(b2 - h1_IIE) + (b2 - theta1)^2
#     sigma1sq <- mean(eif_sigma1sq)
#     list(sigma1sq = sigma1sq, eif_sigma1sq = eif_sigma1sq)
# }
#
# cal.sigma1sq.IDE <- function(med, a_mc, mus_IDE, folds, learners_regressions, control, covariates_cl = FALSE, benchmark_covar = NULL) {
#     h1_IDE <- estimate_h1.IDE(med, a_mc,
#                               mus_IDE, folds, learners_regressions, control, covariates_cl, benchmark_covar)
#
#     Y <- med@data[[med@vars@Y]]
#     A <- med@data[[med@vars@A]]
#     wYresid <- (1*(A==1)/bound(a_mc[, "a(1|m,c)"], tol = 1e-2) - 1*(A==0)/bound(a_mc[, "a(0|m,c)"], tol = 1e-2)) * (Y - mus_IDE$mu_fit_ay_jo)
#     # wYresid <- (1*(A==1)*a_c[, "a(1|c)"]/bound(a_mc[, "a(1|m,c)"]) - 1*(A==0)*a_c[, "a(0|c)"]/bound(a_mc[, "a(0|m,c)"])) * (Y - mus_IDE$mu_fit_ay_jo)
#
#     b2 <- mus_IDE$b_ay_jo # mus[[2]]$b_ay_jo = mus[[1]]$b_ay_jo for IDE
#     theta1 <- mus_IDE$mu_fit_am_jo
#
#     eif_sigma1sq <- 2*(wYresid - h1_IDE) * (b2 - theta1) + (b2 - theta1)^2
#     # eif_sigma1sq <- 2*(wYresid - h1_IDE + b2 - theta1) * (b2 - theta1) + (b2 - theta1)^2
#     # eif_sigma1sq <- 2*(wYresid - h1_IDE) * b2 + (b2 - theta1)^2
#     # eif_sigma1sq <- 2*wYresid*(b2 - h1_IDE) + (b2 - theta1)^2
#     sigma1sq <- mean(eif_sigma1sq)
#     list(sigma1sq = sigma1sq, eif_sigma1sq = eif_sigma1sq)
#
# }

# cal.resY <- function(Y, m_each, mus) {
#
#
#     if (m_each[[1]][["az"]] != "") {
#         resY <- vector("list", length = 3)
#         names(resY) <- glue("mu{1:3}")
#
#         resY[["mu1"]] <- sapply(1, \(x) {
#             resid_y <- (mus[[x]]$b2 - mus[[x]]$mu_fit1)
#             resid_y
#         })
#         resY[["mu2"]] <- sapply(1, \(x) {
#             resid_y <- (mus[[x]]$b3 - mus[[x]]$mu_fit2)
#             resid_y
#         })
#         resY[["mu3"]] <- sapply(1, \(x) {
#             resid_y <- (Y - mus[[x]]$mu_fit3)
#             resid_y
#         })
#
#     }
#
#     if (m_each[[1]][["az"]] == "") {
#         resY <- vector("list", length = 2)
#         names(resY) <- glue("mu{2:3}")
#
#         resY[["mu2"]] <- sapply(1, \(x) {
#             resid_y <- (mus[[x]]$b_ay_jo - mus[[x]]$mu_fit_am_jo)
#             resid_y
#         })
#         resY[["mu3"]] <- sapply(1, \(x) {
#             resid_y <- (Y - mus[[x]]$mu_fit_ay_jo)
#             resid_y
#         })
#
#
#     }
#
#     reduce(resY, cbind)
# }
#
#
# # cal.R2Y <- function(Y, m_each, mus) {
#
#
#     if (m_each[[1]][["az"]] != "") {
#         R2_Y <- vector("list", length = 3)
#         names(R2_Y) <- glue("mu{1:3}")
#
#         R2_Y[["mu1"]] <- sapply(1, \(x) {
#             resid_y <- (mus[[x]]$b2 - mus[[x]]$mu_fit1)
#             r2_y <- max(1 - var(resid_y) / var(mus[[x]]$b2), 0)
#             r2_y
#         })
#
#         R2_Y[["mu2"]] <- sapply(1, \(x) {
#             resid_y <- (mus[[x]]$b3 - mus[[x]]$mu_fit2)
#             r2_y <- max(1 - var(resid_y) / var(mus[[x]]$b3), 0)
#             r2_y
#         })
#         R2_Y[["mu3"]] <- sapply(1, \(x) {
#             resid_y <- (Y - mus[[x]]$mu_fit3)
#             r2_y <- max(1 - var(resid_y) / var(Y), 0)
#             r2_y
#         })
#     }
#     if (m_each[[1]][["az"]] == "") {
#         R2_Y <- vector("list", length = 2)
#         names(R2_Y) <- glue("mu{2:3}")
#
#         R2_Y[["mu2"]] <- sapply(1, \(x) {
#             resid_y <- (mus[[x]]$b_ay_jo - mus[[x]]$mu_fit_am_jo)
#             r2_y <- max(1 - var(resid_y) / var(mus[[x]]$b_ay_jo), 0)
#             r2_y
#         })
#
#         R2_Y[["mu3"]] <- sapply(1, \(x) {
#             resid_y <- (Y - mus[[x]]$mu_fit_ay_jo)
#             r2_y <- max(1 - var(resid_y) / var(Y), 0)
#             r2_y
#         })
#
#
#     }
#
#     R2_Y
# #}

