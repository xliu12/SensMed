sens_each <- function(med, b_thetas, alphas_each, m_each, cf.y = 0.04, cf.m = 0.03, rho2 = 1) {
    # nu2
    # eif_nu2 <- 2 * alpha_each1$alpha3 * alpha_each1$alpha2 - alpha_each1$alpha3^2 
    # eif_nu2 <- list()
    # eif_nu2[[2]] <- (2 *alphas_each[[2]]$alpha3* alphas_each[[2]]$alpha2 - alphas_each[[2]]$alpha3^2) 
    # eif_nu2[[1]] <- (2 *alphas_each[[1]]$alpha3* alphas_each[[1]]$alpha2 - alphas_each[[1]]$alpha3^2)
    eif_nu2 <- (2 *alphas_each[[2]]$alpha3_target * alphas_each[[2]]$alpha2 - 2 *alphas_each[[1]]$alpha3_target * alphas_each[[1]]$alpha2) - (alphas_each[[2]]$alpha3 - alphas_each[[1]]$alpha3)^2
    # nu2 <- mean(eif_nu2)
    nu2 <- mean(eif_nu2)
    if (nu2 < 0) {
        eif_nu2 <- -eif_nu2
        nu2 <- -nu2 # (alphas_each[[2]]$alpha3 - alphas_each[[1]]$alpha3)^2
    }
    psi_nu2 <- eif_nu2 - nu2
    
    # sigma2
    Y <- med@data[[med@vars@Y]]
    eif_sigma2 <- (Y - b_thetas[[1]]$theta_fit3)^2 # b_thetas[[1]]$theta_fit3 = b_thetas[[2]]$theta_fit3 when m_each[[1]]$ay = m_each[[2]]$ay
    sigma2 <- mean(eif_sigma2)
    psi_sigma2 <- eif_sigma2 - sigma2
    
    # S2
    S2 <- sigma2 * nu2
    psi_S2 <- sigma2 * psi_nu2 + nu2 * psi_sigma2
    
    # eif
    # eif_short <- eif_each(med, b_theta_each1, alpha_each1, m_each1)
    eif_short_each <- eif_each_cal(med, b_thetas, alphas_each, m_each)
    names(eif_short_each)
    eif_short <- eif_short_each[[2]] - eif_short_each[[1]]
    theta_short <- mean(eif_short)
    psi_theta_short <- eif_short - theta_short
    
    # bias bound
    bf <- bias.factor(cf.y = cf.y, cf.m = cf.m, rho2 = rho2)
    S  <- sqrt(S2)
    bias_bound <- S*bf
    
    # theta_short \pm bias.bound (plug-in, without debias for S2)
    plugin_long_p <- theta_short + bias_bound
    plugin_long_m <- theta_short - bias_bound
    
    
    # bounds IFs
    psi_bias_bound <- (bf/2) * (1/S) * psi_S2
    psi_long_p    <- psi_theta_short + psi_bias_bound
    psi_long_m    <- psi_theta_short - psi_bias_bound
    
    sens_res <- data.frame(
        theta_short = theta_short,
        plugin_long_p = plugin_long_p,  
        plugin_long_m = plugin_long_m,
        eif_short = eif_short, 
        psi_long_p = psi_long_p,
        psi_long_m = psi_long_m,
        psi_bias_bound = psi_bias_bound, 
        bias_bound = bias_bound
    )
    fit_short <- data.frame(nu2 = nu2, yfit_short =  b_thetas[[1]]$theta_fit3)
    
    list(sens_res = sens_res, fit_short = fit_short)
}
