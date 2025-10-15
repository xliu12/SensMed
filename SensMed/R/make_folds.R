


make.fold_K <- function(data_in, Snames=NULL, cv_folds=4) {

    if (cv_folds <=1) {
        folds1 <- origami::make_folds(data_in,
                                      # fold_fun = origami::folds_vfold,
                                      V = 1)
        folds1[[1]]$training_set <- folds1[[1]]$validation_set
        folds <- folds1


    }

    if (cv_folds > 1) {
        if (is.null(Snames)) {
            folds <- origami::make_folds(data_in,
                                         # fold_fun = origami::folds_vfold,
                                         V = cv_folds)

        }

        if (!is.null(Snames)) {

            K <- match(data_in$S, unique(data_in$S))
            # folds <- origami::make_folds(data_in, V = cv_folds, strata_ids = K)

            data_in$id <- 1:nrow(data_in)
            fold_K <- lapply(unique(K), FUN = function(k=1) {

                if (nrow(data_in[K==k, ]) >= 1) {
                    fk <- origami::make_folds(data_in[K==k, ],
                                              #fold_fun = origami::folds_vfold,
                                              V = cv_folds)
                    fold_k <- fk
                    v <- 1
                    for(v in 1:length(fk)) {
                        fold_k[[v]]$validation_set <- data_in$id[K==k][fk[[v]]$validation_set]
                        fold_k[[v]]$training_set <- data_in$id[K==k][fk[[v]]$training_set]
                    }
                }

                return(fold_k)
            })

            folds <- origami::make_folds(data_in,
                                         fold_fun = origami::folds_vfold,
                                         V = cv_folds)

            for(v in 1:cv_folds) {
                folds[[v]]$validation_set <- unlist(lapply(1:length(fold_K), FUN = function(k=1) {
                    fold_K[[k]][[v]]$validation_set
                }))
                folds[[v]]$training_set <- unlist(lapply(1:length(fold_K), FUN = function(k=1) {
                    fold_K[[k]][[v]]$training_set
                }))
            }


        }

    }
  folds
}

# modified based on R package crumble https://github.com/nt-williams/crumble
combine_folds <- function(pred_folds, folds) {
        ind <- Reduce(c, lapply(folds, function(x) x[["validation_set"]]))
        name1 <- names(pred_folds[[1]])[1]
        # ind[order(ind)]
        out <- list(valid_row = ind[order(ind)])

        for (name1 in names(pred_folds[[1]]) ) {
            vals <- Reduce(rbind, lapply(1:length(pred_folds), \(x) do.call(cbind, pred_folds[[2]][[name1]]) ))
            out[[name1]] <- vals[order(ind)]
        }

        out
    }

combine_folds_theta <- function(vals_folds, folds) {
  ind <- Reduce(c, lapply(folds, function(x) x[["validation_set"]]))

  if(is.data.frame(vals_folds[[1]])) {
      vals <- do.call(rbind, vals_folds)
      if (is.null(dim(vals))) {
          all_folds <- vals[order(ind)]
      } else {
          all_folds <- vals[order(ind), ]
      }
      names(all_folds) <- names(vals_folds[[1]])
  } else {
      all_folds <- vector("list", length = length( vals_folds[[1]] ))
      i <- 1
      for (i in 1:length( vals_folds[[1]] ) ) {
          vals <- do.call(rbind, lapply(vals_folds, \(x) x[[ i ]]))
          vals <- drop(vals)
          if (is.null(dim(vals))) {
              all_folds[[ i ]] <- vals[order(ind)]
          } else {
              all_folds[[ i ]] <- vals[order(ind), ]
          }
          names(all_folds)[i] <- names(vals_folds[[1]])[i]
      }
  }

  all_folds
}

combine_folds_alpha <- function(vals_folds, folds) {
  ind <- Reduce(c, lapply(folds, function(x) x[["validation_set"]]))

  all_folds <- vector("list", length = length( vals_folds[[1]] ))

  i <- 1
  for (i in 1:length( vals_folds[[1]] ) ) {
    vals <- do.call(rbind, lapply(vals_folds, \(x) do.call(cbind, x[[ i ]]) ))
    all_folds[[ i ]] <- data.frame(vals[order(ind), ])
    names(all_folds)[i] <- names(vals_folds[[1]])[i]
  }


  all_folds
}
