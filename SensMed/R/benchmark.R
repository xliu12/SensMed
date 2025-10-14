#' Obtain sensitivity parameter values based on observed covariates
#'
#' Compute gain metrics for an observed covariate (or a set of observed covariates) used as a benchmark.
#'
#' @param short_list An object named "short_list" returned by the function \code{run.sensmed}.
#' @param benchmark_covariates A \code{character} vector containing the column names of the observed covariates used as the benchmark covariate.
#'
#' @export
#'
run.benchmark <- function(short_list,
                             benchmark_covariates){

    short_list$benchmark_covariates <- benchmark_covariates

    list2env(short_list, envir = environment())

    # get short components without Xb; long-short regression \theta_{1,\uu,\ss,-b} (A,X_\ss) on the benchmark_covariates
    short_XBomitted <- estimate.short(data,
                                      treatment,
                                      outcome,
                                      mediators,
                                      covariates = base::setdiff(covariates, benchmark_covariates),
                                      cluster_id,
                                      learners_outcome,
                                      learners_weights,
                                      learners_else,
                                      # for benchmark, the long regressions/weights are those from short
                                      long_mus_IIE = mus_IIE,
                                      long_alphas_IIE = alphas_IIE,
                                      long_mus_IDE = mus_IDE,
                                      long_alphas_IDE = alphas_IDE,
                                      # other settings
                                      conf.level,
                                      benchmark_covariates,
                                      control)

    # check the benchmark covariate is omitted in short, but included in the long-short
    # short_XBomitted$med@vars@C
    # short_XBomitted$med_U@vars@C

    # for IIE ------------
    # Gains for alphas
    Den_R2_alpha1 <- mean((alphas_IIE$alpha1 - short_XBomitted$alphas_IIE$alpha1)^2) + short_XBomitted$cal_nu2s_IIE$nu2$nu1
    # alternative
    # Den_R2_alpha1 <- cal_nu2s_IIE$nu2$nu1
    R2_alpha1 <- short_XBomitted$cal_nu2s_IIE$nu2$nu1 / Den_R2_alpha1

    Den_R2_alpha2 <- mean((alphas_IIE$alpha2 - short_XBomitted$alphas_IIE$alpha2)^2) + short_XBomitted$cal_nu2s_IIE$nu2$nu2
    # alternative
    # mean(alphas_IIE$alpha2^2)# mean(short_XBomitted$alphas_IIE$alpha2^2)
    # Den_R2_alpha2 <- cal_nu2s_IIE$nu2$nu2 - 2*mean((alphas_IIE$alpha1 - short_XBomitted$alphas_IIE$alpha1) * short_XBomitted$alphas_IIE$alpha2_target)
    R2_alpha2 <- short_XBomitted$cal_nu2s_IIE$nu2$nu2 / Den_R2_alpha2

    Gain_alpha1 <- (1-R2_alpha1) / R2_alpha1
    Gain_alpha2 <- (1-R2_alpha2) / R2_alpha2

    # Gains for thetas
    # R2_theta1 <- mean((short_XBomitted$long.short_mus_IIE$mu_fit_am_jo - short_XBomitted$mus_IIE$mu_fit_am_jo)^2) / short_XBomitted$cal_sigma2s_IIE$sigma2s$sigma1

    R2_theta1 <- 1 - mean((short_XBomitted$mus_IIE$b_ay_jo - short_XBomitted$long.short_mus_IIE$mu_fit_am_jo)^2) / short_XBomitted$cal_sigma2s_IIE$sigma2s$sigma1

    R2_theta2 <- 1 - cal_sigma2s_IIE$sigma2s$sigma2 / short_XBomitted$cal_sigma2s_IIE$sigma2s$sigma2

    Gain_theta1 <- pmax(R2_theta1/(1-R2_theta1), 0)
    Gain_theta2 <- pmax(R2_theta2/(1-R2_theta2), 0)

    # Correlation
    rho1 <- cor((short_XBomitted$mus_IIE$mu_fit_am_jo - short_XBomitted$long.short_mus_IIE$mu_fit_am_jo), (short_XBomitted$alphas_IIE$alpha1 - alphas_IIE$alpha1))

    rho2 <- cor((short_XBomitted$mus_IIE$mu_fit_ay_jo - mus_IIE$mu_fit_ay_jo), (short_XBomitted$alphas_IIE$alpha2 - alphas_IIE$alpha2))


    benchmark_IIE <- c(Gain_alpha1, Gain_alpha2, Gain_theta1, Gain_theta2, rho1, rho2)
    names(benchmark_IIE) <- c("Gain_alpha1", "Gain_alpha2", "Gain_theta1", "Gain_theta2", "rho1", "rho2")


    # for IDE ------------
    # Gains for alphas
    Den_R2_alpha1 <- mean((alphas_IDE$alpha1 - short_XBomitted$alphas_IDE$alpha1)^2) + short_XBomitted$cal_nu2s_IDE$nu2$nu1
    # alternative
    # Den_R2_alpha1 <- cal_nu2s_IDE$nu2$nu1
    R2_alpha1 <- short_XBomitted$cal_nu2s_IDE$nu2$nu1 / Den_R2_alpha1

    Den_R2_alpha2 <- mean((alphas_IDE$alpha2 - short_XBomitted$alphas_IDE$alpha2)^2) + short_XBomitted$cal_nu2s_IDE$nu2$nu2
    # alternative
    # mean(alphas_IDE$alpha2^2)# mean(short_XBomitted$alphas_IDE$alpha2^2)
    # Den_R2_alpha2 <- cal_nu2s_IDE$nu2$nu2 - 2*mean((alphas_IDE$alpha1 - short_XBomitted$alphas_IDE$alpha1) * short_XBomitted$alphas_IDE$alpha2_target)
    R2_alpha2 <- short_XBomitted$cal_nu2s_IDE$nu2$nu2 / Den_R2_alpha2

    Gain_alpha1 <- (1-R2_alpha1) / R2_alpha1
    Gain_alpha2 <- (1-R2_alpha2) / R2_alpha2

    # Gains for thetas
    # R2_theta1 <- mean((short_XBomitted$long.short_mus_IDE$mu_fit_am_jo - short_XBomitted$mus_IDE$mu_fit_am_jo)^2) / short_XBomitted$cal_sigma2s_IDE$sigma2s$sigma1

    R2_theta1 <- 1 - mean((short_XBomitted$mus_IDE$b_ay_jo - short_XBomitted$long.short_mus_IDE$mu_fit_am_jo)^2) / short_XBomitted$cal_sigma2s_IDE$sigma2s$sigma1

    R2_theta2 <- 1 - cal_sigma2s_IDE$sigma2s$sigma2 / short_XBomitted$cal_sigma2s_IDE$sigma2s$sigma2

    Gain_theta1 <- pmax(R2_theta1/(1-R2_theta1), 0)
    if (mean((short_XBomitted$mus_IDE$b_ay_jo - short_XBomitted$long.short_mus_IDE$mu_fit_am_jo)^2) < 1e-4) {
        Gain_theta1 <- 0
    }
    Gain_theta2 <- pmax(R2_theta2/(1-R2_theta2), 0)

    # Correlation
    rho1 <- cor((short_XBomitted$mus_IDE$mu_fit_am_jo - short_XBomitted$long.short_mus_IDE$mu_fit_am_jo), (short_XBomitted$alphas_IDE$alpha1 - alphas_IDE$alpha1))

    rho2 <- cor((short_XBomitted$mus_IDE$mu_fit_ay_jo - mus_IDE$mu_fit_ay_jo), (short_XBomitted$alphas_IDE$alpha2 - alphas_IDE$alpha2))


    benchmark_IDE <- c(Gain_alpha1, Gain_alpha2, Gain_theta1, Gain_theta2, rho1, rho2)
    names(benchmark_IDE) <- c("Gain_alpha1", "Gain_alpha2", "Gain_theta1", "Gain_theta2", "rho1", "rho2")


    rbind(IIE = benchmark_IIE, IDE = benchmark_IDE)
}
