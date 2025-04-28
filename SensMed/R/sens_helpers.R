
cal.nu2 <- function(alphas) {
    eif_nu2 <- list()
    # nu2
    eif_nu2[["alpha3"]] <- sapply(alphas[1:2], 
                                     \(x) 2 * x$alpha3_target * x$alpha2 - x$alpha3^2)
    eif_nu2[["alpha2"]] <- sapply(alphas[1:2], 
                                     \(x) 2 * x$alpha2_target * x$alpha1 - x$alpha2^2)
    eif_nu2[["alpha1"]] <- sapply(alphas[1:2], 
                                     \(x) 2 * x$alpha1_target * 1 - x$alpha1^2)
    
    ## currently, focusing on the individual-average effects
    # individual average for the nuisance parameters
    nu2 <- map(eif_nu2, \(x) abs(apply(x, 2, mean)))
    
    list(nu2 = nu2, eif_nu2 = eif_nu2)
}

cal.nu2.IIE_M1 <- function(alphas) {
    # alphas=alphas_each_IIE_M1_r
    natural <- list()
    if (length(alphas) == 2) {
        natural <- map(1:3, \(k) alphas[[2]][[glue("alpha{k}")]] - alphas[[1]][[glue("alpha{k}")]] ) 
        # shifted <- map(1:3, \(k) alphas[[2]][[glue("alpha{k}_target")]] - alphas[[1]][[glue("alpha{k}_target")]] ) 
        
        shifted <- list()
        shifted[["alpha1_target"]] <- alphas[[2]][["alpha1_target"]] - alphas[[1]][["alpha1_target"]] 
        shifted[["alpha2_target"]] <- alphas[[2]][["alpha2_target"]] * alphas[[2]][["alpha1"]]  - alphas[[1]][["alpha2_target"]] * alphas[[1]][["alpha1"]]
        shifted[["alpha3_target"]] <- alphas[[2]][["alpha3_target"]] * alphas[[2]][["alpha2"]]  - alphas[[1]][["alpha3_target"]] * alphas[[1]][["alpha2"]]
    }
    if (length(alphas) == 1) {
        natural <- map(1:3, \(k) alphas[[1]][[glue("alpha{k}")]] ) 
        # shifted <- map(1:3, \(k) alphas[[1]][[glue("alpha{k}_target")]] ) 
        shifted <- list()
        shifted[["alpha1_target"]] <- alphas[[1]][["alpha1_target"]] 
        shifted[["alpha2_target"]] <- alphas[[1]][["alpha2_target"]] * alphas[[1]][["alpha1"]]
        shifted[["alpha3_target"]] <- alphas[[1]][["alpha3_target"]] * alphas[[1]][["alpha2"]]
    }
    
    names(natural) <- glue("alpha{1:3}")
    names(shifted) <- glue("alpha{1:3}_target")
    
    
    eif_nu2 <- vector("list", length = 3)
    names(eif_nu2) <- glue("nu{1:3}") # make sure all things are ordered as 1,2,3
    # nu2
    # eif_nu2[["nu1"]] <- 2 * shifted$alpha1_target * 1 - natural$alpha1^2
    # eif_nu2[["nu2"]] <- 2 * shifted$alpha2_target * natural$alpha1 - natural$alpha2^2
    # eif_nu2[["nu3"]] <- 2 * shifted$alpha3_target * natural$alpha2 - natural$alpha3^2
    
    eif_nu2[["nu1"]] <- 2 * shifted$alpha1_target - natural$alpha1^2
    eif_nu2[["nu2"]] <- 2 * shifted$alpha2_target - natural$alpha2^2
    eif_nu2[["nu3"]] <- 2 * shifted$alpha3_target - natural$alpha3^2
    
    ## currently, focusing on the individual-average effects
    # individual average for the nuisance parameters
    
    nu2 <- map(eif_nu2, \(x) abs(mean(x)))
    # nu2 <- vector("list", length = 3)
    names(nu2) <- glue("nu{1:3}")
    # nu2[["nu1"]] <- mean(natural$alpha1^2)
    # nu2[["nu2"]] <- mean(natural$alpha2^2)
    # nu2[["nu3"]] <- mean(natural$alpha3^2)
    
    
    list(nu2 = nu2, eif_nu2 = eif_nu2)
}

# alphas <- alphas_each_IIE_Mjo[1:2]

cal.nu2.IIE_Mjo <- function(alphas) {
    # alphas=alphas_each_IIE_M1_r
    natural <- list()
    if (length(alphas) == 2) {
        natural <- map(1:3, \(k) alphas[[2]][[glue("alpha{k}")]] - alphas[[1]][[glue("alpha{k}")]] ) 
        # shifted <- map(1:3, \(k) alphas[[2]][[glue("alpha{k}_target")]] - alphas[[1]][[glue("alpha{k}_target")]] ) 
        shifted <- list()
        # shifted[["alpha1_target"]] <- alphas[[2]][["alpha1_target"]] - alphas[[1]][["alpha1_target"]] 
        shifted[["alpha2_target"]] <- alphas[[2]][["alpha2_target"]] - alphas[[1]][["alpha2_target"]] 
        # shifted[["alpha3_target"]] <- alphas[[2]][["alpha3_target"]] * alphas[[2]][["alpha2"]]  - alphas[[1]][["alpha3_target"]] * alphas[[1]][["alpha2"]]
        shifted[["alpha3_target"]] <- (alphas[[2]][["alpha3_target"]] - alphas[[1]][["alpha3_target"]]) * (alphas[[2]][["alpha2"]] - alphas[[1]][["alpha2"]])
    }
    if (length(alphas) == 1) {
        natural <- map(1:3, \(k) alphas[[1]][[glue("alpha{k}")]] ) 
        # shifted <- map(1:3, \(k) alphas[[1]][[glue("alpha{k}_target")]] ) 
        shifted <- list()
        
        shifted[["alpha2_target"]] <- alphas[[1]][["alpha2_target"]] 
        shifted[["alpha3_target"]] <- alphas[[1]][["alpha3_target"]] * alphas[[1]][["alpha2"]]
    }
    
    names(natural) <- glue("alpha{1:3}")
    # names(shifted) <- glue("alpha{1:3}_target")
    
    
    eif_nu2 <- vector("list", length = 2)
    names(eif_nu2) <- glue("nu{2:3}")
    # nu2 
    eif_nu2[["nu2"]] <- 2 * shifted$alpha2_target - natural$alpha2^2 # alpha1 = alpha2 for Mjo
    eif_nu2[["nu3"]] <- 2 * shifted$alpha3_target - natural$alpha3^2
    
    ## currently, focusing on the individual-average effects
    # individual average for the nuisance parameters
    nu2 <- map(eif_nu2, \(x) abs(mean(x)))
    # nu2 <- map(eif_nu2, \(x) (mean(x)))
    nu2plugin <- map(natural, \(x) (mean(x^2)))
    
    list(nu2 = nu2, eif_nu2 = eif_nu2, nu2plugin = nu2plugin)
}

cal.nu2.IDE <- function(alphas) {
    # alphas=alphas_each_IIE_M1_r
    natural <- list()
    if (length(alphas) == 2) {
        natural <- list()
        natural[["alpha2"]] <- alphas[[2]][["alpha2"]] # for IDE, alpha_2 = alpha_1 = 1(A=0)/p(A|X)
        natural[["alpha3"]] <- alphas[[2]][["alpha3"]] - alphas[[1]][["alpha3"]]
        
        # natural <- map(1:3, \(k) alphas[[2]][[glue("alpha{k}")]] - alphas[[1]][[glue("alpha{k}")]] ) 
        # shifted <- map(1:3, \(k) alphas[[2]][[glue("alpha{k}_target")]] - alphas[[1]][[glue("alpha{k}_target")]] ) 
        shifted <- list()
        shifted[["alpha2_target"]] <- alphas[[2]][["alpha2_target"]] # - alphas[[1]][["alpha2_target"]] 
        shifted[["alpha3_target"]] <- alphas[[2]][["alpha2"]] * alphas[[2]][["alpha3_target"]] - alphas[[2]][["alpha2"]] *  alphas[[1]][["alpha3_target"]]
        # shifted[["alpha3_target"]] <- (alphas[[2]][["alpha3_target"]] - alphas[[1]][["alpha3_target"]]) * (alphas[[2]][["alpha2"]] - alphas[[1]][["alpha2"]])
    }
    if (length(alphas) == 1) {
        natural <- map(1:3, \(k) alphas[[1]][[glue("alpha{k}")]] ) 
        names(natural) <- glue("alpha{1:3}")
        # shifted <- map(1:3, \(k) alphas[[1]][[glue("alpha{k}_target")]] ) 
        shifted <- list()
        
        shifted[["alpha2_target"]] <- alphas[[1]][["alpha2_target"]] 
        shifted[["alpha3_target"]] <- alphas[[1]][["alpha3_target"]] * alphas[[1]][["alpha2"]]
    }
    
    # 
    # names(shifted) <- glue("alpha{1:3}_target")
    
    
    eif_nu2 <- vector("list", length = 2)
    names(eif_nu2) <- glue("nu{2:3}")
    # nu2 
    eif_nu2[["nu2"]] <- 2 * shifted$alpha2_target - natural$alpha2^2 # alpha1 = alpha2 for Mjo
    eif_nu2[["nu3"]] <- 2 * shifted$alpha3_target - natural$alpha3^2
    
    ## currently, focusing on the individual-average effects
    # individual average for the nuisance parameters
    nu2 <- map(eif_nu2, \(x) abs(mean(x)))
    nu2plugin <- map(natural, \(x) (mean(x^2)))
    
    list(nu2 = nu2, eif_nu2 = eif_nu2, nu2plugin = nu2plugin)
    
}

cal.resY <- function(Y, m_each, mus) {
    
    
    if (m_each[[1]][["az"]] != "") {
        resY <- vector("list", length = 3)
        names(resY) <- glue("mu{1:3}")
        
        resY[["mu1"]] <- sapply(1, \(x) {
            resid_y <- (mus[[x]]$b2 - mus[[x]]$mu_fit1)
            resid_y
        }) 
        resY[["mu2"]] <- sapply(1, \(x) {
            resid_y <- (mus[[x]]$b3 - mus[[x]]$mu_fit2)
            resid_y
        })
        resY[["mu3"]] <- sapply(1, \(x) {
            resid_y <- (Y - mus[[x]]$mu_fit3)
            resid_y
        })

    }
    
    if (m_each[[1]][["az"]] == "") {
        resY <- vector("list", length = 2)
        names(resY) <- glue("mu{2:3}")
        
        resY[["mu2"]] <- sapply(1, \(x) {
            resid_y <- (mus[[x]]$b_ay_jo - mus[[x]]$mu_fit_am_jo)
            resid_y
        }) 
        resY[["mu3"]] <- sapply(1, \(x) {
            resid_y <- (Y - mus[[x]]$mu_fit_ay_jo)
            resid_y
        })
        
        
    }
    
    reduce(resY, cbind)
}


cal.R2Y <- function(Y, m_each, mus) {
    
    
    if (m_each[[1]][["az"]] != "") {
        R2_Y <- vector("list", length = 3)
        names(R2_Y) <- glue("mu{1:3}")
        
        R2_Y[["mu1"]] <- sapply(1, \(x) {
            resid_y <- (mus[[x]]$b2 - mus[[x]]$mu_fit1)
            r2_y <- max(1 - var(resid_y) / var(mus[[x]]$b2), 0)
            r2_y
        })
        
        R2_Y[["mu2"]] <- sapply(1, \(x) {
            resid_y <- (mus[[x]]$b3 - mus[[x]]$mu_fit2)
            r2_y <- max(1 - var(resid_y) / var(mus[[x]]$b3), 0)
            r2_y
        })
        R2_Y[["mu3"]] <- sapply(1, \(x) {
            resid_y <- (Y - mus[[x]]$mu_fit3)
            r2_y <- max(1 - var(resid_y) / var(Y), 0)
            r2_y
        })
    }
    if (m_each[[1]][["az"]] == "") {
        R2_Y <- vector("list", length = 2)
        names(R2_Y) <- glue("mu{2:3}")
        
        R2_Y[["mu2"]] <- sapply(1, \(x) {
            resid_y <- (mus[[x]]$b_ay_jo - mus[[x]]$mu_fit_am_jo)
            r2_y <- max(1 - var(resid_y) / var(mus[[x]]$b_ay_jo), 0)
            r2_y
        })
        
        R2_Y[["mu3"]] <- sapply(1, \(x) {
            resid_y <- (Y - mus[[x]]$mu_fit_ay_jo)
            r2_y <- max(1 - var(resid_y) / var(Y), 0)
            r2_y
        })
        
        
    }
    
    R2_Y
}


