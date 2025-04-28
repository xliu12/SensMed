
inference.iid <- function(eifs, conf.level = 0.95) {
  n <- nrow(eifs)
  
  estSE <- data.frame(
    estimate = colMeans(eifs),
    std_error = map_dbl(1:ncol(eifs), ~sqrt( var(eifs[[.]])/n ))
  ) %>% 
    mutate(
      p_value = pnorm(abs(estimate) / std_error, lower.tail = FALSE) * 2,
      ci_low = estimate - qnorm(1-(1-conf.level)/2) * std_error, 
           ci_high = estimate + qnorm(1-(1-conf.level)/2)*std_error )
  
  names(estSE)[names(estSE)%in% c("ci_low", "ci_high")] <- c(glue("{scales::percent(conf.level)} CI_low"),glue("{scales::percent(conf.level)} CI_high"))
  
  estSE
}

inference.cl <- function(eifs, S = NA_character_, average = "individual", conf.level = 0.95) {

  if (!is.na(S[1])) {
    eifs_df <- eifs %>% 
      mutate(S = S) 
    eif_clmean <- eifs_df %>% group_by(S) %>% summarise( across(everything(), mean) )
    sizes <- eifs_df %>% group_by(S) %>% summarise( nj = n() ) %>% pull(nj)
    K <- length(unique(S))
    mean_size <- mean(sizes)
    
    if (average=="individual") {
      estimate <- colMeans(eifs)
      std_error <- apply(eif_clmean[, -1], 2, \(x) sqrt(var(x*sizes/mean_size)/K) )  
      # variance <- var(eif_clmean[1:K, estimand] * sizes/mean_size)/K
    }
    
    if (average=="cluster") {
      # estimate <- sum(eif_clmean[1:K, estimand])/K # mean(eif_clmean[, estimand])
      estimate <- colMeans(eif_clmean[, -1])
      std_error <- apply(eif_clmean[, -1], 2, \(x) sqrt(var(x)/K) )  
      # variance <- var(eif_clmean[1:K, estimand])/K
    }
    
    estSE <- data.frame(
      estimate = estimate,
      std_error = std_error
    ) %>% 
      mutate(
        p_value = pnorm(abs(estimate) / std_error, lower.tail = FALSE) * 2,
        ci_low = estimate - qnorm(1-(1-conf.level)/2)*std_error, 
        ci_high = estimate + qnorm(1-(1-conf.level)/2)*std_error )
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
