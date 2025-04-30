
shift_data <- function(data, trt, shift = 1, center = NA_character_) {
    names(shift) <- trt
    for (a in trt) {
        data[[a]] <- shift[[a]]
        # data[[glue("{a}0")]] <- 1 - shift[[a]]

        if (!is.na(center)) {
            data[[glue("{a}_cwc")]] <- shift[[a]] - data[[ glue("{a}_clmean") ]]
            # data[[glue("{a}0_cwc")]] <- 1 - shift[[a]] - data[[ glue("{a}_clmean") ]]
        }
    }
    data
}

normalize <- function(x) {
    # if (is_normalized(x)) return(x)
    x / mean(x)
}


estimate_control <- function(crossfit_folds = 4,
                             mlr3superlearner_folds = 10,
                             mlr3superlearner_discrete = TRUE,
                             epochs = 10L,
                             learning_rate = 0.01,
                             batch_size = 64,
                             device = c("cpu"), #, "cuda", "mps"
                             indepM = TRUE,
                             folds_indepM = 4L,
                             riesz_folds = 5L,
                             learners_riesz = c("nn")
                             ) {
    list(
        crossfit_folds = crossfit_folds,
        mlr3superlearner_folds = mlr3superlearner_folds,
        # zprime_folds = zprime_folds,
        mlr3superlearner_discrete = mlr3superlearner_discrete,
        epochs = epochs,
        learning_rate = learning_rate,
        batch_size = as.numeric(batch_size),
        device = torch::torch_device(match.arg(device)),
        indepM = indepM,
        folds_indepM = folds_indepM,
        riesz_folds = riesz_folds,
        learners_riesz = learners_riesz
    )
}

if.binary <- function(x) {
    if ( all(unique(x) %in% c(0,1)) ) {
        TRUE
    } else {
        FALSE
    }
}

bound <- function(vals, tol = 0.01) {
    vals[vals < tol] <- tol
    vals[vals > 1 - tol] <- 1 - tol
    return(vals)
}

bound2 <- function(vals, low = -1, high = 1) {
    vals[vals < low] <- low
    vals[vals > high] <- high
    return(vals)
}

scale_to_unit <- function(vals) {
    vals_scaled <- (vals - min(vals)) / (max(vals) - min(vals))
    return(vals_scaled)
}

bound.ipw <- function(ipw, n) {
    # (Gruber et al., 2022) https://doi.org/10.1093/aje/kwac087
    upper <- sqrt(n)*log(n)/5
    ipw[ipw > upper] <- upper
    ipw
}


# # update interactions
# update.interactions <- function(data, tt, R, ttR) {
#
#     data[, ttR] <- data[[tt]] * data[[R]]
#
#     data
# }



one_hot_encode <- function(data_tmp) {
    # tmp <- data[, vars, drop = FALSE]
    # as.data.frame(model.matrix(~ ., data = tmp))[, -1, drop = FALSE]
    as.data.frame(model.matrix(~ ., data = data_tmp))[, -1, drop = FALSE]
}
as_torch <- function(data, device) {
    torch::torch_tensor(as.matrix(data), dtype = torch::torch_float(), device = device)
}
