
# from crumble https://github.com/nt-williams/crumble

estimate_alpha_each <- function(med, m_each, folds, nn_module, control, which.M, benchmark_covar = NULL) {
    alpha_each_list <- vector("list", control$crossfit_folds)
    i <- 1
    cli::cli_progress_step("Computing alpha n density ratios... {i}/{control$crossfit_folds} folds")
    for (i in seq_along(folds)) {
        train <- training(med, folds, i)
        valid <- validation(med, folds, i)
        
        m_each1 <- m_each[[1]]
        alpha_each_list[[i]] <- lapply(
            m_each,
            \(m_each1) alpha_each(train, valid, med@vars, nn_module, m_each1, control, which.M, benchmark_covar)
        )
        
        names(alpha_each_list[[i]]) <- map_chr(m_each, ~paste(glue( '{names(.)}{gsub("data_", "", .)}' ), collapse = "_") )
        #unlist(lapply(m_each, \(x) paste0(gsub("data_", "", x), collapse = "")))
        cli::cli_progress_update()
    }
    
    cli::cli_progress_done()
    
    alpha_allfolds <- combine_folds_alpha(alpha_each_list, folds)
    alpha_allfolds
}


alpha_each <- function(train, valid, vars, architecture, m_each1, control, which.M, benchmark_covar = NULL) {
    
    .f1 <- \(alpha, dl) alpha(dl[[ m_each1["am"] ]])
    .f2 <- \(alpha, dl) alpha(dl[[ m_each1["az"] ]])
    .f3 <- \(alpha, dl) alpha(dl[[ m_each1["ay"] ]])
    alpha1 <- hat.alpha(
        train = train,
        valid = valid,
        vars = c(vars@A, vars@C) ,
        architecture = architecture,
        .f = .f1,
        target = m_each1["am"],
        weights = NULL,
        control = control
    )
    if (m_each1["az"] == "") {
        alpha2 <- alpha1
    }
    if (m_each1["az"] != "") {
        alpha2 <- hat.alpha(
            train = train,
            valid = valid,
            vars = c(vars@A, vars@C, vars@M[-which.M]) ,
            architecture = architecture,
            .f = .f2,
            target = m_each1["az"],
            weights = alpha1$train,
            control = control
        )
    }
    
    if (length(benchmark_covar) > 0) {
        shortC <- setdiff(vars@C, benchmark_covar)
    } else {
        shortC <- vars@C
    }
    
    alpha3 <- hat.alpha(
        train = train,
        valid = valid,
        vars = c(vars@A, shortC, vars@M) ,
        architecture = architecture,
        .f = .f3,
        target = m_each1["ay"],
        weights = alpha2$train,
        control = control
    )
    # paste(substr(m_each1, 6, 6), collapse = "")
    alpha_list_valid <- list(#zmy = paste(substr(m_each1, 6, 6), collapse = ""),
         alpha1 = alpha1$valid,
         alpha2 = alpha2$valid,
         alpha3 = alpha3$valid,
         alpha1_target = alpha1$valid_target,
         alpha2_target = alpha2$valid_target,
         alpha3_target = alpha3$valid_target)
    alpha_list_valid
}


