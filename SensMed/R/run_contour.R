#' Make a contour plot using sensitivity parameter values
#'
#' Compute estimates and confidence intervals for the bounds of the interventional indirect or direct effect at specified sensitivity parameter values.
#'
#' @param effect A \code{string} specifying the name of the effect. Options include  "IIE" (interventional indirect effect) or "IDE" (interventional direct effect).
#' @param target (optional) A \code{string} indicating the values to be plotted, which can be one of the column names of the \code{data.frame} object "sens_res" returned by the function \code{run.sensmed}. Default is \code{target = "theta_m.90%CI_low"}, the lower confidence limit of the lower bound.
#' @param bench.cf.y (optional) A \code{numeric(1)} value to be used for both the R-squared sensitivity parameters for the regressions, to be added as an additional point on the plot.
#' @param bench.cf.m (optional) A \code{numeric(1)} value to be used for both the R-squared sensitivity parameters for the weights, to be added as an additional point on the plot.
#' @param rho (optional) A \code{numeric} vector of the correlation sensitivity parameters, that is, rho_\{1\} and rho_\{2\}.
#' @param short_list An object named "short_list" returned by the function \code{run.sensmed}.
#'
#'
#' @export
#'

run.contour <- function(effect = "IIE", target = "theta_m.90%CI_low", short_list, bench.cf.y = NULL, bench.cf.m = NULL, rho = c(1, 1), cf_grid = NULL, y_axis = NULL, x_axis = NULL, contour_levels = NULL, plot_zero_confounding = FALSE, only_YM_confounded = FALSE) {

    rho2 <- rho^2
    # benchmark -------
    # Benchmark_res <- benchmark.effect(short_list, benchmark_covariates)

    ## Contour plots -----

    list2env(short_list, envir = environment())

    if (effect == "IIE") {
        eif <- eif_IIE; cal_nu2s <- cal_nu2s_IIE; cal_sigma2s <- cal_sigma2s_IIE; estimate <- estimates["eif_IIE", ]
    }
    if (effect == "IDE") {
        eif <- eif_IDE; cal_nu2s <- cal_nu2s_IDE; cal_sigma2s <- cal_sigma2s_IDE; estimate <- estimates["eif_IDE", ]
    }

    sens_grid <- z_grid <- list()
    add_benchmark_max <- add_benchmark_min <- plot_estimate <- list()

    # cf_fix <- c(cf.m2 = 0.2, cf.y2 = 0.1)

    if (is.null(cf_grid)) {
        maxGain <- max(dplyr::select(as.data.frame(Benchmark_res), starts_with("Gain")))
        maxgrid <- min(maxGain + 0.2, 0.99)
        # cf_grid <- expand.grid(cf.y2 = c(0, seq(0.01, maxgrid, length.out = 10)),
        #                        cf.y1 = c(Benchmark_res[effect, "Gain_theta1"]),
        #                        cf.m2 = c(Benchmark_res[effect, "Gain_alpha2"]),
        #                        cf.m1 = c(0, seq(0.01, maxgrid, length.out = 10)))
        cf_grid <- expand.grid(cf.y2 = c(0, seq(0.01, maxgrid, length.out = 10)),
                               cf.m1 = c(0, seq(0.01, maxgrid, length.out = 10)))
        cf_grid$cf.y1 <- cf_grid$cf.y2
        cf_grid$cf.m2 <- cf_grid$cf.m1
        y_axis <- "cf.y2"
        x_axis <- "cf.m1"
    }


    # cf_grid <- expand.grid(cf.y = contour_grid$y,
    #                        cf.m = contour_grid$x)

    if (length(rho2) < 2) {
        rho2 <- c(rho2, rho2)
    }
    ##   contour -----
    sens_grid <- lapply(
        1:nrow(cf_grid),
        \(x=1) sens_Mjo(med, eif, cal_nu2s, cal_sigma2s,
                        cf.y = as.numeric(cf_grid[x, c("cf.y1", "cf.y2")]),
                        cf.m = as.numeric(cf_grid[x, c("cf.m1", "cf.m2")]), rho2, only_YM_confounded, estimate)[["sens_estimates"]]
    ) %>% do.call(rbind, .) %>% as.data.frame()

    # if (is.null(y_axis)) {
    #     y_axis <- names(cf_grid)[apply(cf_grid, 2, \(x) length(unique(x))) > 1][1]
    # }
    # if (is.null(x_axis)) {
    #     x_axis <- names(cf_grid)[apply(cf_grid, 2, \(x) length(unique(x))) > 1][2]
    # }
    sens_grid$cf.y <- sens_grid[[y_axis]]
    sens_grid$cf.m <- sens_grid[[x_axis]]

    z_grid <- sens_grid %>%
        # mutate(cf.y = cf.y1, cf.m = cf.m1) %>%
        dplyr::select(cf.y, cf.m, all_of(target)) %>%
        pivot_wider(names_from = cf.m, values_from = target) %>%
        as.matrix(.)
    z_grid <- z_grid[, -1]

    # bench.cf.y <- Benchmark_res[effect, c("Gain_theta1", "Gain_theta2")]
    # bench.cf.m <- Benchmark_res[effect, c("Gain_alpha1", "Gain_alpha2")]

    add_benchmark_max <- sens_Mjo(med, eif, cal_nu2s, cal_sigma2s, cf.y = as.numeric(bench.cf.y), cf.m = as.numeric(bench.cf.m), rho2, only_YM_confounded, estimate)[["sens_estimates"]]

    # add_benchmark_max$`theta_m.90%CI_low`
    # add_benchmark_min <- sens_Mjo(med, eif, cal_nu2s, cal_sigma2s, cf.y = min(bench.cf.y), cf.m = min(bench.cf.m), rho2, only_YM_confounded, estimate)[["sens_estimates"]]

    plot_estimate <- sens_Mjo(med, eif, cal_nu2s, cal_sigma2s, cf.y = 0, cf.m = 0, rho2, only_YM_confounded, estimate)[["sens_estimates"]]
    plot_estimate <- plot_estimate[[target]]


    if (is.null(contour_levels)) {
        contour_levels <- quantile(z_grid, probs = seq(0.1, 0.9, by=0.05)) %>% round(., 2)
    }
    contour_plot(grid_values.x = unique(sens_grid$cf.m),
                 grid_values.y = unique(sens_grid$cf.y),
                 z_axis = z_grid,
                 levels = contour_levels,
                 labels = NULL,
                 threshold = 0,
                 col.contour = "blue4")

    if (plot_zero_confounding) {
        label.unadjusted <- "No \nConfounding"

        points(0, 0, pch = 17, col = "black", cex = 0.8)
        label.bump.x <- 0.02
        label.bump.y <- 0.03
        text(0 + label.bump.x, 0 + label.bump.y,
             paste0(label.unadjusted,
                    "\n(", round(plot_estimate, 3), ")"),
             cex = 0.6)
    }


    sensemakr::add_bound_to_contour(r2dz.x = add_benchmark_max[[x_axis]],
                                    r2yz.dx = add_benchmark_max[[y_axis]],
                                    bound_value = add_benchmark_max[[target]],
                                    bound_label = "Maximum \nBenchmark",
                                    label.text = TRUE,
                                    label.bump.x = 0.03,
                                    label.bump.y = 0.02,
                                    cex.label.text = 0.7,
                                    round = 2)

    contour_resutls <- list(cf_grid = cf_grid,
                            z_grid = z_grid,
                            add_benchmark_max = add_benchmark_max,
                            add_benchmark_min = add_benchmark_min,
                            plot_estimate = plot_estimate)

    # output ----------------

    contour_resutls
}
