#' Run sensitivity analysis at specified sensitivity parameter values
#'
#' Compute estimates and confidence intervals for the bounds of the interventional indirect or direct effect at specified sensitivity parameter values.
#'
#' @param effect A \code{string} specifying the name of the effect. Options include  "IIE" (interventional indirect effect) or "IDE" (interventional direct effect)
#' @param cf.y (optional) A \code{numeric} vector of the R-squared sensitivity parameters for the regressions, that is, R^2_\{Y-theta_\{2,s\} ~ theta_2 -theta_\{2,s\} \} and R^2_\{bM-theta_\{1,s\} ~ theta_\{1,u,s\} -theta_\{1,s\} \}.
#' @param cf.m (optional) A \code{numeric} vector of the R-squared sensitivity parameters for the weights, that is, 1-R^2_\{alpha_2 ~ alpha_\{2,s\}\} and 1-R^2_\{alpha_1 ~ alpha_\{1,s\}\}.
#' @param rho (optional) A \code{numeric} vector of the correlation sensitivity parameters, that is, rho_\{1\} and rho_\{2\}.
#' @param short_list An object named "short_list" returned by the function \code{run.sensmed}.
#'
#' @export
#'


med.sens <- function(effect = "IIE", cf.y = c(0.26, 0.26), cf.m = c(0.26, 0.26), rho = c(1, 1), short_list) {

  rho2 <- rho^2
  list2env(short_list, envir = environment())

  ## confidence bounds at given sens params ----
  if (effect == "IIE") {
    sens_IIE <- sens_Mjo(med, eif_IIE, cal_nu2s_IIE, cal_sigma2s_IIE, cf.y, cf.m, rho2, only_YM_confounded = FALSE, estimate = estimates["eif_IIE", ])
    sens_res <- sens_IIE$sens_estimates
  }
  if (effect == "IDE") {
    sens_IDE <- sens_Mjo(med, eif_IDE, cal_nu2s_IDE, cal_sigma2s_IDE, cf.y, cf.m, rho2, only_YM_confounded = FALSE, estimate = estimates["eif_IDE", ])
    sens_res <- sens_IDE$sens_estimates
  }

    rownames(sens_res) <- effect

  sens_res
}
