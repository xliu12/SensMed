# from R package crumble 
draw_indepM <- function(med, data, which.M, folds_indepM = 5, Sname_dummies = NULL) {
    AC <- one_hot_encode(data[, c(med@vars@A, med@vars@C)])
    if (!is.null(Sname_dummies)) {
        AC <- one_hot_encode(data[, c(med@vars@A, med@vars@C, Sname_dummies)])
    }
    L <- length(med@vars@M)
    indepM_dat <- data.frame(matrix(NA, nrow = nrow(AC), ncol = L))
    names(indepM_dat) <- med@vars@M
    # if (nrow(AC) > 500) {
    #     folds_indepM <- 4 # round(nrow(AC) / 500)
    # }
    drawM_folds <- origami::make_folds(data, V = folds_indepM)
    if (folds_indepM == 1) {
        drawM_folds[[1]]$training_set <- drawM_folds[[1]]$validation_set
    }
    
    # l <- 1
    for (l in which.M) {
        Ml <- data[, med@vars@M[l], drop = FALSE]
        AC_l <- AC
        # AC_l <- AC[sample.int(nrow(AC), replace = FALSE), ]
        
        permute <- function(i) {
            zp <- data.frame(matrix(NA, nrow = nrow(AC_l), ncol = ncol(Ml)))
            names(zp) <- names(Ml)
            P <- linear_permutation(AC_l[i, ])
            for (z in names(zp)) {
                zp[i, z] <- as.vector(P %*% as.matrix(Ml[i, z]))
            }
            zp
        }
        
        permuted <- vector("list", folds_indepM)
        i <- 1
        cli::cli_progress_step("Permuting Z-prime variables... {i}/{folds_indepM} tasks")
        for (i in seq_along(drawM_folds)) {
            permuted[[i]] <- permute(drawM_folds[[i]]$validation_set) |>
                as.list()
            cli::cli_progress_update()
        }
        
        cli::cli_progress_done()
        
        indepMl <- revert_list(permuted) |>
            lapply(\(x) Reduce(data.table::fcoalesce, x)) |>
            data.frame()
        
        indepM_dat[, l] <- indepMl
    }
    
    indepM_dat
}

# Not used modified based on R package crumble 
draw_indepM.within <- function(med, which.M, folds_indepM = 1) {
    L <- length(which.M)
    indepM_dat <- data.frame(matrix(NA, nrow = nrow(med@data), ncol = L))
    names(indepM_dat) <- med@vars@M[which.M]
    
    j <- 1
    cli::cli_progress_step("Permuting to create independent mediator variables... {j}/{length(unique(med@data$S))} clusters")
    
    folds_indepM <- 1
    for (j in unique(med@data$S)) {
        AC <- one_hot_encode(med@data[med@data$S==j, c(med@vars@A, med@vars@C)])
        # if (nrow(AC) > 500) {
        #     folds_indepM <- 4
        # }
        drawM_folds <- origami::make_folds(AC, V = folds_indepM)
        if (folds_indepM == 1) {
            drawM_folds[[1]]$training_set <- drawM_folds[[1]]$validation_set
        }
        
        l <- 2
        for (l in which.M) {
            Ml <- med@data[med@data$S==j, med@vars@M[l], drop = FALSE]
            # AC_l <- AC[sample.int(nrow(AC), replace = FALSE), ]
            
            permute <- function(i) {
                zp <- data.frame(matrix(NA, nrow = nrow(AC), ncol = ncol(Ml)))
                names(zp) <- names(Ml)
                P <- linear_permutation(AC[i, ])
                for (z in names(zp)) {
                    zp[i, z] <- as.vector(P %*% as.matrix(Ml[i, z]))
                }
                zp
            }
            
            permuted <- vector("list", length(drawM_folds))
            i <- 1
            
            for (i in seq_along(drawM_folds)) {
                permuted[[i]] <- permute(drawM_folds[[i]]$validation_set) |>
                    as.list()
                cli::cli_progress_update()
            }
            
            
            indepMl <- revert_list(permuted) |>
                lapply(\(x) Reduce(data.table::fcoalesce, x)) |>
                data.frame()
            
            indepM_dat[med@data$S==j, ] <- indepMl
        }
        cli::cli_progress_update()
    }
    cli::cli_progress_done()
    
    if (length(grep("_cwc$", med@vars@M)) > 0) {
        indepM_cwc <- indepM_dat
    } else {
        indepM_cwc <- indepM_dat - med@data[, glue("{med@vars@M[which.M]}_clmean")]
    }
    
    names(indepM_cwc) <- glue("{med@vars@M[which.M]}_cwc")
    data.frame(indepM_dat, indepM_cwc)
    # indepM_dat
}




linear_permutation <- function(data) {
    D <- dist(data)
    D <- D / max(D)
    d  <- as.vector(t(as.matrix(D)))
    n <- nrow(data)
    rows <- c(as.numeric(gl(n, n, n^2)), as.numeric(gl(n, n, n^2)) + n)
    cols <- c(1:(n^2), unlist(lapply(0:(n - 2), function(j) j + seq(1, n^2, n))), (0:(n - 1)*(n + 1) + 1))
    A <- Matrix::sparseMatrix(i = rows, j = cols, x = 1)
    b <- Matrix::Matrix(Matrix::sparseVector(i = 1:(2*n - 1), x = 1, length = 2*n), ncol = 1)
    matrix(Rsymphony::Rsymphony_solve_LP(d, A, dir = rep("==", nrow(b)), rhs = b)$solution, n, n)
}

sample_M <- function(M) { # univariate
    set.seed(12345)
    vals_M <- unique(M)
    vals_M[sample.int(length(vals_M), size = length(M), replace = TRUE)]
}
revert_list <- function(ls) { # @Josh O'Brien
    # get sub-elements in same order
    x <- lapply(ls, `[`, names(ls[[1]]))
    # stack and reslice
    apply(do.call(rbind, x), 2, as.list)
}


add_indepM <- function(med, control, which.M, cwc = FALSE, Sname_dummies) {
    indepM_dat <- draw_indepM(med, med@data, which.M, control$folds_indepM, Sname_dummies)
    # med@data_M1_indepM2 <- med@data
    # med@data_M2_indepM1 <- med@data
    # med@data_0_M1_indepM2 <- med@data_0
    # med@data_0_M2_indepM1 <- med@data_0
    # med@data_1_M1_indepM2 <- med@data_1
    # med@data_1_M2_indepM1 <- med@data_1
    med@data_indepM <- med@data
    med@data_0_indepM <- med@data_0
    med@data_1_indepM <- med@data_1
    med@data_indepM[, med@vars@M[which.M]] <- indepM_dat[, which.M]
    med@data_0_indepM[, med@vars@M[which.M]] <- indepM_dat[, which.M]
    med@data_1_indepM[, med@vars@M[which.M]] <- indepM_dat[, which.M]
    
    if (cwc) {
        med@data_indepM[, glue("{med@vars@M[which.M]}_cwc")] <- indepM_dat[, which.M] - med@data_indepM[, glue("{med@vars@M[which.M]}_clmean")]
        med@data_0_indepM[, glue("{med@vars@M[which.M]}_cwc")] <- indepM_dat[, which.M] - med@data_0_indepM[, glue("{med@vars@M[which.M]}_clmean")]
        med@data_1_indepM[, glue("{med@vars@M[which.M]}_cwc")] <- indepM_dat[, which.M] - med@data_1_indepM[, glue("{med@vars@M[which.M]}_clmean")]
    }
    med
}


add_indepM.within <- function(med, control, which.M, if_indepM = FALSE) {
    if (if_indepM == FALSE) {
        newM_dat <- draw_indepM.within(med, which.M, control$folds_indepM)
        med@data_indepM <- med@data
        med@data_0_indepM <- med@data_0
        med@data_1_indepM <- med@data_1
        # names(newM_dat)
        Mother <- c(med@vars@M[which.M], glue("{med@vars@M[which.M]}_cwc"))
        # Mother <- med@vars@M[which.M] 
        med@data_indepM[, Mother] <- newM_dat[, Mother]
        med@data_0_indepM[, Mother] <- newM_dat[, Mother]
        med@data_1_indepM[, Mother] <- newM_dat[, Mother]  
    }
    if (if_indepM == TRUE) {
        med@data_indepM <- med@data
        med@data_0_indepM <- med@data_0
        med@data_1_indepM <- med@data_1
    }
    med
}
