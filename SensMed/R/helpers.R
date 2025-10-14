

#'
#' @export
#'

estimate_control <- function(crossfit_folds = 5,
                             mlr3superlearner_folds = 5
                             ) {
    list(
        crossfit_folds = crossfit_folds,
        mlr3superlearner_folds = mlr3superlearner_folds
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


normalize <- function(x) {
    # if (is_normalized(x)) return(x)
    x / mean(x)
}


