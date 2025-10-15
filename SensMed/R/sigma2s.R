
cal.sigma1sq.IIE <- function(med, a_mc, mus_IIE, mus_Mjo) {
    
    Y <- med@data[[med@vars@Y]]
    A <- med@data[[med@vars@A]]
    w1_mc <- 1*(A==1)/pmax(a_mc[, "a(1|m,c)"], 1e-2)
    e_mc <- mus_Mjo$y01$b_am_jo * a_mc[, "a(0|m,c)"] + mus_Mjo$y11$b_am_jo * a_mc[, "a(1|m,c)"]

    eif_sigma1sq <- 2 * w1_mc * (Y - mus_IIE$mu_fit_ay_jo) * (mus_IIE$b_ay_jo - e_mc) + (mus_IIE$b_ay_jo - mus_IIE$mu_fit_am_jo)^2

    sigma1sq <- mean(eif_sigma1sq)
    list(sigma1sq = sigma1sq, eif_sigma1sq = eif_sigma1sq)
}

cal.sigma1sq.IDE <- function(med, a_mc, mus_IDE, mus_Mjo) {
    
    Y <- med@data[[med@vars@Y]]
    A <- med@data[[med@vars@A]]
    w_mc_ay1 <- 1*(A==1)/pmax(a_mc[, "a(1|m,c)"], 1e-2)
    w_mc_ay0 <- 1*(A==0)/pmax(a_mc[, "a(0|m,c)"], 1e-2)
    e_mc <- (mus_Mjo$y01$b_am_jo - mus_Mjo$y00$b_am_jo) * a_mc[, "a(0|m,c)"] + (mus_Mjo$y11$b_am_jo - mus_Mjo$y10$b_am_jo) * a_mc[, "a(1|m,c)"]

    eif_sigma1sq <- 2 * (w_mc_ay1 * (Y - mus_IDE$mu_fit_ay_jo) - w_mc_ay0 * (Y - mus_IDE$mu_fit_ay_jo)) * (mus_IDE$b_ay_jo - e_mc) + (mus_IDE$b_ay_jo - mus_IDE$mu_fit_am_jo)^2

    mean(eif_sigma1sq)
    
    sigma1sq <- mean(eif_sigma1sq)
    list(sigma1sq = sigma1sq, eif_sigma1sq = eif_sigma1sq)

}

cal.sigma2s.IIE_Mjo <- function(Y, mus_IIE) {
    resY <- R2_Y <- list("sigma1" = NA, "sigma2" = NA)

    resY[["sigma2"]] <- (Y - mus_IIE$mu_fit_ay_jo) # is the same as mus[[1]]$mu_fit_ay_jo
    R2_Y[["sigma2"]] <- max(1 - var(resY[["sigma2"]]) / var(Y), 0)

    b2 <- mus_IIE$b_ay_jo # mus[[2]]$b_ay_jo = mus[[1]]$b_ay_jo for IIE
    theta1 <- mus_IIE$mu_fit_am_jo
    resY[["sigma1"]] <- (b2 - theta1)
    R2_Y[["sigma1"]] <- max(1 - var(resY[["sigma1"]]) / var(b2), 0)

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

    eif_sigma2s <- map(resY, \(x) x^2)
    sigma2s <- map(eif_sigma2s, \(x) mean(x))

    list(resY = resY, R2_Y = R2_Y, sigma2s = sigma2s, eif_sigma2s = eif_sigma2s)
}

