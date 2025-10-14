eif_each <- function(med, mus1, alpha_each1, m_each1) {

    Y <- med@data[[med@vars@Y]]

    eif_each <- alpha_each1$alpha2 * (Y - mus1$mu_fit_ay_jo) +
        alpha_each1$alpha1 * (mus1$b_ay_jo - mus1$mu_fit_am_jo) +
        mus1$b_am_jo

    eif_each1 <- eif_each
    if (!is.null(m_each1)) {
        eif_each1 <- data.frame(eif_each)
        names(eif_each1) <- paste0("y", paste(substr(m_each1, 6, 6), collapse = ""))
    }
    eif_each1
}


eif_each_cal <- function(med, mus, alphas_each, m_each) {
    if (!is.null(m_each)) {
        eifs <- lapply(1:length(m_each), \(x = 1) eif_each(med, mus[[x]], alphas_each[[x]], m_each[[x]]))
        eif_df <- data.frame(do.call(cbind, eifs))
    } else {
        eif_df <- eif_each(med, mus, alphas_each, m_each)
    }

    eif_df
}

