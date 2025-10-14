
inference.iid <- function(eifs, conf.level = 0.95) {
  n <- nrow(eifs)

  estSE <- data.frame(
    estimate = colMeans(eifs),
    std_error = map_dbl(1:ncol(eifs), ~sqrt( var(eifs[[.]])/n ))
  ) %>%
    dplyr::mutate(
      p_value = pnorm(abs(estimate) / std_error, lower.tail = FALSE) * 2,
      ci_low = estimate - qnorm(1-(1-conf.level)/2) * std_error,
           ci_high = estimate + qnorm(1-(1-conf.level)/2)*std_error )

  names(estSE)[names(estSE)%in% c("ci_low", "ci_high")] <- c(glue("{scales::percent(conf.level)} CI_low"),glue("{scales::percent(conf.level)} CI_high"))

  estSE
}

inference.cl <- function(eifs, S = NA_character_, average = "individual", conf.level = 0.95, nboot=2) {

  if (!is.na(S[1])) {
    eifs_df <- eifs %>%
      mutate(S = S)

    eif_clmean <- eifs_df %>% group_by(S) %>% summarise( across(everything(), mean) )
    sizes <- eifs_df %>% group_by(S) %>% summarise( nj = n() ) %>% pull(nj)
    K <- length(unique(S))
    mean_size <- mean(sizes)

    if (average=="individual") {
        estimate <- colMeans(eifs)
        std_error <- apply(eif_clmean[, -1], 2, \(x) sqrt(var(x*sizes/mean_size)/(K-2)) )
        # variance <- var(eif_clmean[1:K, estimand] * sizes/mean_size)/K
    }

    if (average=="cluster") {
        # estimate <- sum(eif_clmean[1:K, estimand])/K # mean(eif_clmean[, estimand])
        estimate <- colMeans(eif_clmean[, -1])
        std_error <- apply(eif_clmean[, -1], 2, \(x) sqrt(var(x)/(K-2)) )
        # variance <- var(eif_clmean[1:K, estimand])/K
    }

    each.boot <- function(b=1) {
        ind <- sample(S, size = K, replace = TRUE) # K: number of clusters
        eifs_b <- lapply(ind, \(j) eifs_df[eifs_df$S == j, ]) %>% do.call(rbind, .)
        if (average=="individual") {
            estb <- colMeans(select(eifs_b, !S))
        }
        if (average=="cluster") {
            eif_b_clmean <- eifs_b %>% group_by(S) %>% summarise( across(everything(), mean) )
            estb <- colMeans(eif_b_clmean[, -1, drop=FALSE])
        }
        estb
    }

    boots_r <- replicate(nboot, each.boot())
    if(is.matrix(boots_r)) {
        boots <- matrix(boots_r, nrow=1)
    } else {
        boots <- boots_r
    }

    # std_error <- apply(boots, 1, sd, na.rm=TRUE)
    sd.nonextreme <- function(v) {
        # v1 <- v[v > quantile(v, 0.001) & v < quantile(v, 0.999)]
        v1 <- v
        sd(v1, na.rm = TRUE)
    }

    estSE <- data.frame(
      estimate = estimate,
      std_error = std_error,
      CLBstd_error = apply(boots, 1, sd.nonextreme)
    ) %>%
      mutate(
          CLBci_low = apply(boots, 1, quantile, (1-conf.level)/2, na.rm=TRUE),
          CLBci_high = apply(boots, 1, quantile, 1-(1-conf.level)/2, na.rm=TRUE),
          CLBNci_low = estimate - qt(1-(1-conf.level)/2, df = K-2)*CLBstd_error,
          CLBNci_high = estimate + qt(1-(1-conf.level)/2, df = K-2)*CLBstd_error,
        # p_value = pnorm(abs(estimate) / std_error, lower.tail = FALSE) * 2,
          p_value = pt(abs(estimate) / std_error, df = K-2, lower.tail = FALSE) * 2,
        ci_low = estimate - qt(1-(1-conf.level)/2, df = K-2)*std_error,
        ci_high = estimate + qt(1-(1-conf.level)/2, df = K-2)*std_error )
  }

  if (is.na(S[1])) {
    n <- nrow(eifs)

    estSE <- data.frame(
      estimate = colMeans(eifs),
      std_error = map_dbl(1:ncol(eifs), ~sqrt( var(eifs[[.]])/n ))
    ) %>%
      mutate(
        p_value = pnorm(abs(estimate) / std_error, lower.tail = FALSE) * 2,
        ci_low = estimate - qnorm(1-(1-conf.level)/2) * std_error,
        ci_high = estimate + qnorm(1-(1-conf.level)/2)*std_error )
  }
  names(estSE)[names(estSE)%in% c("ci_low", "ci_high")] <- c(glue("{scales::percent(conf.level)} CI_low"),glue("{scales::percent(conf.level)} CI_high"))

  estSE

}
