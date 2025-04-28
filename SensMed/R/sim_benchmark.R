


gendata <- function(n = 1e4) {
    
    set.seed(123)
    
    # long ---------------------
    p.A <- function(X) {
        plogis(1 * X[, 1]) 
    }
    p.M1 <- function(A, X) {
        plogis(1 * X[, 1] + 0.5*(A - 0.5))
    }
    p.M2 <- function(A, X) {
        plogis(1 * X[, 1] + 0.5*(A - 0.5))
    }
    theta.3 <- function(A, M2, M1, X) {
        (1*X[, 1] + 0.5*(A-0.5) + 0.5*(M1-0.5) + 0.5*(M2-0.5))
    }
    
    X <- replicate(3, rbinom(n, 1, 0.5) - 0.5 )
    A <- rbinom(n, 1, prob = p.A(X))
    M1 <- rbinom(n, 1, prob = p.M1(A, X))
    M2 <- rbinom(n, 1, prob = p.M2(A, X))
    Y <- rnorm(n, mean = theta.3(A, M2, M1, X), 
               sd = 1)
    dat <- data.frame(id = 1:n,
        X = X,
        A = A,
        M1 = M1,
        M2 = M2,
        Y = Y
    )
    
    m_each1 <- list(
        IIE_M1_m = c(az = "data_0_indepM", am = "data_0_indepM", ay = "data_1_indepM"),
        IIE_M1_p = c(az = "data_0_indepM", am = "data_1_indepM", ay = "data_1_indepM")
    )
    
    ay <- 1
    az <- 0
    
    Ypseudo_3 <- theta.3(A = rep(ay, n), M2, M1, X) 
    
    theta.2 <- function(A, M1, X) {
        theta.3(ay, M2=p.M2(A, X), M1, X)
    }
    
    # az <- 0
    Ypseudo_2 <- theta.2(A = rep(az, n), M1, X)
    
    theta.1 <- function(A, X) {
        theta.2(A, M1=p.M1(A, X), X)
        # theta.3(ay, M2=p.M2(az, X), M1=p.M1(A, X), X)
    }
    
    
    alpha.1 <- function(A, X, am1) {
        p_am1 <- p.A(X)*(am1 == 1) + (1-p.A(X))*(am1 == 0)
        1*(A==am1) / p_am1
    }
    # mean((1*(A==1) / p.A(X)) * theta.1(A, X, ay, az)) - mean(theta.1(A = rep(1, n), X, ay, az))
    # mean(alpha.1(A, X, 1) * theta.1(A, X, ay, az)) -  
    #     (mean(theta.1(A = rep(1, n), X, ay, az)) )
    
    alpha.2 <- function(A, M1, X, az = 0, am1) {
        p_az <- p.A(X)*(az == 1) + (1-p.A(X))*(az == 0)
        ipwA <- 1*(A==az) / p_az
        
        # am <- c(p=1, m=0)
        pM1_a <- p.M1(A=rep(am1,n), X) * (M1 == 1) + (1 - p.M1(A=rep(am1,n), X)) * (M1 == 0)
        pM1_A <- p.M1(A, X) * (M1 == 1) + (1 - p.M1(A, X)) * (M1 == 0)
        
        rmpw <- pM1_a / pM1_A
        
        ipwA * rmpw
    }
    
    # mean(alpha.2(A, M1, X, az, am1=1) * theta.2(A, M1, X, ay))
    alpha.3 <- function(A, M2, M1, X, ay = 1, az = 0, am1 = 1) {
        p_ay <- p.A(X)*(ay == 1) + (1-p.A(X))*(ay == 0)
        ipwA <- 1*(A==ay) / p_ay
        
        pM1_a <- p.M1(A=rep(am1,n), X) * (M1 == 1) + (1 - p.M1(A=rep(am1,n), X)) * (M1 == 0)
        pM1_A <- p.M1(A, X) * (M1 == 1) + (1 - p.M1(A, X)) * (M1 == 0)
        
        pM2_a <- p.M2(A=rep(az,n), X) * (M2 == 1) + (1 - p.M2(A=rep(az,n), X)) * (M2 == 0)
        pM2_A <- p.M2(A, X) * (M2 == 1) + (1 - p.M2(A, X)) * (M2 == 0)
        
        ipwA * (pM1_a/pM1_A) * (pM2_a/pM2_A)
    }
    
    # mean(alpha.3(A, M2, M1, X, ay, az, am1 = 1) * theta.3(A, M2, M1, X))
    # mean(theta.1(1, X, ay, az))
    
    # true estimate ---------------
    mu_list <- map(c(0, 1), \(x) data.frame(mu_fit3 = theta.3(A, M2, M1, X),
                                            mu_fit2 = theta.2(A, M1, X),
                                            mu_fit1 = theta.1(A, X),
                                            b3 = Ypseudo_3,
                                            b2 = Ypseudo_2,
                                            b1 = theta.1(A = x, X) ) ) 
    
    # alpha_list <- data.frame(alpha3 = alpha.3(A, M2, M1, X, ay, az, am1),
    #                    alpha2 = alpha.2(A, M1, X, az, am1),
    #                    alpha1 = alpha.1(A, X, am1) )
    
    alpha_list <- map(c(0, 1), \(x) data.frame(
        alpha3 = alpha.3(A, M2, M1, X, ay, az, am1=x),
        alpha2 = alpha.2(A, M1, X, az, am1=x),
        alpha1 = alpha.1(A, X, am1=x) )) # %>% reduce(`-`)
    
    
    med <- med_data(
        data = dat,
        vars = med_vars(
            A = "A",
            Y = "Y",
            M = c("M1", "M2"),
            C = names(dat)[grep("^X", names(dat))],
            cluster_id = NA_character_,
            Wj = NA_character_,
            S = NA_character_,
            id = "id"
        ),
        weights = rep(1, nrow(dat)),
        a0 = 0,
        a1 = 0
    )
    
    eif <- eif_each_cal(med, mu_list, alpha_list, m_each1)
    eif[["IIE_M1"]] <- eif[[2]] - eif[[1]]
    
    estimate <- inference.cl(eif, S = med@vars@cluster_id, conf.level = 0.95)
    
    out <- mget(ls(envir = environment()))
    
    return(out)
}


par.fit <- function(med, benchmark = "X.1") {
    X <- med@data[, med@vars@C]
    if (!is.na(benchmark)) {
        Xs <- med@data[, setdiff(med@vars@C, benchmark)]
    } else{
        Xs <- med@data[, med@vars@C]
    }
    
    A <- med@data[, med@vars@A]
    M1 <- med@data[, med@vars@M[1]]
    M2 <- med@data[, med@vars@M[2]]
    Y <- med@data[, med@vars@Y]
    
    fit_A <- glm(A ~ ., data = data.frame(A, Xs), family = binomial())
    fit_M1 <- glm(M1 ~ ., data = data.frame(M1, A, Xs), family = binomial())
    fit_M2 <- glm(M2 ~ ., data = data.frame(M2, A, Xs), family = binomial())
    
    fit_A_l <- glm(A ~ ., data = data.frame(A, X), family = binomial())
    fit_M1_l <- glm(M1 ~ ., data = data.frame(M1, A, X), family = binomial())
    fit_M2_l <- glm(M2 ~ ., data = data.frame(M2, A, X), family = binomial())
    
    Yfamily <- ifelse(length(unique(Y)) > 2, "gaussian", "binomial")
    
    theta.3l <- function(A, M2, M1, X) {
        fit_Yl <- glm(Y ~ ., data = data.frame(Y, M2, M1, A, X), family = Yfamily)
        
        predict(fit_Yl, type = "response")
    }
    
    theta.3s <- function(A, M2, M1, Xs) {
        fit_Ys <- glm(Y ~ ., data = data.frame(Y, M2, M1, A, Xs), family = Yfamily)
        
        predict(fit_Ys, type = "response")
    }
    
    Ypseudo3 <- theta.3s(A = rep(ay, n), M2, M1, Xs) 
    
    
    theta.2s <- function(A, M1, Xs) {
        fit_Ypseudo3_s <- glm(Ypseudo3 ~ ., 
                              data = data.frame(Ypseudo3, M1, A, Xs), family = "gaussian")
        predict(fit_Ypseudo3_s, type = "response")
    }
    theta.2l.3s <- function(A, M1, X) {
        fit_Ypseudo3_l <- glm(Ypseudo3 ~ ., 
                              data = data.frame(Ypseudo3, M1, A, X), family = "gaussian")
        predict(fit_Ypseudo3_l, type = "response")
    }
    # az <- 0
    Ypseudo2 <- theta.2s(A = rep(az, n), M1, Xs)
    
    
    theta.1s <- function(A, Xs) {
        fit_Ypseudo2_s <- glm(Ypseudo2 ~ ., 
                            data = data.frame(Ypseudo2, A, Xs), family = "gaussian")
        
        predict(fit_Ypseudo2_s, type = "response")
    }
    theta.1l.2s <- function(A, X) {
        fit_Ypseudo2_l <- glm(Ypseudo2 ~ ., 
                              data = data.frame(Ypseudo2, A, X), family = "gaussian")
        
        predict(fit_Ypseudo2_l, type = "response")
    }
    
    # am <- c(p = 1, m = 0)
    
    pAs <- predict(fit_A, type = "response")
    p.M1s <- function(A, Xs) {
        predict(fit_M1, newdata = data.frame(A, Xs), type = "response")
    }
    p.M2s <- function(A, Xs) {
        predict(fit_M2, newdata = data.frame(A, Xs), type = "response")
    }
    
    pA_l <- predict(fit_A_l, type = "response")
    p.M1_l <- function(A, X) {
        predict(fit_M1_l, newdata = data.frame(A, X), type = "response")
    }
    p.M2_l <- function(A, X) {
        predict(fit_M2_l, newdata = data.frame(A, X), type = "response")
    }
    
    alpha.1s <- function(A, Xs, am1) {
        p_am1 <- pAs*(am1 == 1) + (1-pAs)*(am1 == 0)
        1*(A==am1) / p_am1
    }
    
    alpha.1l <- function(A, X, am1=1) {
        p_am1 <- pA_l*(am1 == 1) + (1-pA_l)*(am1 == 0)
        1*(A==am1) / p_am1
    }
    
    # mean(alpha.1l(A, X, am1=1)^2)
    # mean((1*(A==1) / p.A(X)) * theta.1(A, X, ay, az)) - mean(theta.1(A = rep(1, n), X, ay, az))
    # mean(alpha.1(A, X, 1) * theta.1(A, X, ay, az)) -  
    #     (mean(theta.1(A = rep(1, n), X, ay, az)) )
    
    alpha.2s <- function(A, M1, Xs, az = 0, am1) {
        p_az <- pAs*(az == 1) + (1-pAs)*(az == 0)
        ipwA <- 1*(A==az) / p_az
        
        # am <- c(p=1, m=0)
        pM1_a <- p.M1s(A=rep(am1,n), Xs) * (M1 == 1) + (1 - p.M1s(A=rep(am1,n), Xs)) * (M1 == 0)
        pM1_A <- p.M1s(A, Xs) * (M1 == 1) + (1 - p.M1s(A, Xs)) * (M1 == 0)
        
        rmpw <- pM1_a / pM1_A
        
        ipwA * rmpw
    }
    
    alpha.2l <- function(A, M1, X, az = 0, am1=1) {
        p_azl <- pA_l*(az == 1) + (1-pA_l)*(az == 0)
        ipwAl <- 1*(A==az) / p_azl
        
        # am <- c(p=1, m=0)
        pM1_al <- p.M1_l(A=rep(am1,n), X) * (M1 == 1) + (1 - p.M1_l(A=rep(am1,n), X)) * (M1 == 0)
        pM1_Al <- p.M1_l(A, X) * (M1 == 1) + (1 - p.M1_l(A, X)) * (M1 == 0)
        
        rmpwl <- pM1_al / pM1_Al
        
        ipwAl * rmpwl
    }
    
   mean(theta.2s(A, M1, Xs) * (alpha.2l(A, M1, X, az = 0, am1=1) - alpha.2s(A, M1, Xs, az = 0, am1=1))) 
   
   mean((alpha.2l(A, M1, X, az = 0, am1=1) - alpha.2s(A, M1, Xs, az = 0, am1=1))^2) 
   
   mean((alpha.2s(A, M1, Xs, az = 0, am1=0))^2) 
   
   mean((alpha.2l(A, M1, X, az = 0, am1=0))^2) - mean((alpha.2s(A, M1, Xs, az = 0, am1=0))^2) 
    
   mean((alpha.1l(A, X, am1=1) - alpha.1s(A, Xs, am1=1))^2) 
   mean((alpha.1l(A, X, am1=1))^2) - mean((alpha.1s(A, Xs, am1=1))^2) 
   
   
    # mean(alpha.2(A, M1, X, az, am1=1) * theta.2(A, M1, X, ay))
    alpha.3s <- function(A, M2, M1, Xs, ay = 1, az = 0, am1 = 1) {
        p_ay <- pAs*(ay == 1) + (1-pAs)*(ay == 0)
        ipwA <- 1*(A==ay) / p_ay
        
        pM1_a <- p.M1s(A=rep(am1,n), Xs) * (M1 == 1) + (1 - p.M1s(A=rep(am1,n), Xs)) * (M1 == 0)
        pM1_A <- p.M1s(A, Xs) * (M1 == 1) + (1 - p.M1s(A, Xs)) * (M1 == 0)
        
        pM2_a <- p.M2s(A=rep(az,n), Xs) * (M2 == 1) + (1 - p.M2s(A=rep(az,n), Xs)) * (M2 == 0)
        pM2_A <- p.M2s(A, Xs) * (M2 == 1) + (1 - p.M2s(A, Xs)) * (M2 == 0)
        
        ipwA * (pM1_a/pM1_A) * (pM2_a/pM2_A)
    }
    
    alpha.3l <- function(A, M2, M1, X, ay = 1, az = 0, am1 = 1) {
        p_ay <- pA_l*(ay == 1) + (1-pA_l)*(ay == 0)
        ipwA <- 1*(A==ay) / p_ay
        
        pM1_a <- p.M1_l(A=rep(am1,n), X) * (M1 == 1) + (1 - p.M1_l(A=rep(am1,n), X)) * (M1 == 0)
        pM1_A <- p.M1_l(A, X) * (M1 == 1) + (1 - p.M1_l(A, X)) * (M1 == 0)
        
        pM2_a <- p.M2_l(A=rep(az,n), X) * (M2 == 1) + (1 - p.M2_l(A=rep(az,n), X)) * (M2 == 0)
        pM2_A <- p.M2_l(A, X) * (M2 == 1) + (1 - p.M2_l(A, X)) * (M2 == 0)
        
        ipwA * (pM1_a/pM1_A) * (pM2_a/pM2_A)
    }
    
    # short estimate -------------
    m_each1
    
    mu_s_list <- map(c(0, 1), \(x) data.frame(mu_fit3 = theta.3s(A, M2, M1, Xs),
                                              mu_fit2 = theta.2s(A, M1, Xs),
                                              mu_fit1 = theta.1s(A, Xs),
                                              b3 = Ypseudo_3,
                                              b2 = Ypseudo_2,
                                              b1 = theta.1s(A = x, Xs) )) 
    
    mu_long_short <- data.frame(mu_fit3 = theta.3l(A, M2, M1, X),
                               mu_fit2 = theta.2l.3s(A, M1, X),
                               mu_fit1 = theta.1l.2s(A, X) )
    
    # alpha_s_list <- list(alpha3 = alpha.3s(A, M2, M1, Xs, ay, az, am1),
    #                    alpha2 = alpha.2s(A, M1, Xs, az, am1),
    #                    alpha1 = alpha.1s(A, Xs, am1) )
    
    alpha_s_list <- map(c(0, 1), \(x) data.frame(
        alpha3 = alpha.3s(A, M2, M1, Xs, ay, az, am1=x),
        alpha2 = alpha.2s(A, M1, Xs, az, am1=x),
        alpha1 = alpha.1s(A, Xs, am1=x) ) ) #%>%  reduce(`-`)
    
    alpha_l_list <- map(c(0, 1), \(x) data.frame(
        alpha3 = alpha.3l(A, M2, M1, X, ay, az, am1=x),
        alpha2 = alpha.2l(A, M1, X, az, am1=x),
        alpha1 = alpha.1l(A, X, am1=x) ) ) #%>%  reduce(`-`)
    
    eif_l <- eif_each_cal(med, mu_list, alpha_l_list, m_each1)
    eif_l[["IIE_M1"]] <- eif_l[[2]] - eif_l[[1]]
    estimate_l <- inference.cl(eif_l, S = med@vars@cluster_id, conf.level = 0.95)
    
    eif_s <- eif_each_cal(med, mu_s_list, alpha_s_list, m_each1)
    eif_s[["IIE_M1"]] <- eif_s[[2]] - eif_s[[1]]
    estimate_s <- inference.cl(eif_s, S = med@vars@cluster_id, conf.level = 0.95)
    
    - estimate_s$estimate + estimate_l$estimate
    # sensitivity ------------------
    map(c(1,2), \(x=1) {
        bias3 <- mean((mu_s_list[[x]]$mu_fit3 - mu_long_short$mu_fit3) * (alpha_s_list[[x]]$alpha3 - alpha_l_list[[x]]$alpha3))
        bias2 <- mean((mu_s_list[[x]]$mu_fit2 - mu_long_short$mu_fit2) * (alpha_s_list[[x]]$alpha2 - alpha_l_list[[x]]$alpha2))
        bias1 <- mean((mu_s_list[[x]]$mu_fit1 - mu_long_short$mu_fit1) * (alpha_s_list[[x]]$alpha1 - alpha_l_list[[x]]$alpha1))
        bias <- bias3 + bias2 + bias1
    })
    
    
    
    map(c(1,2), \(x) {
        bias3 <- mean((mu_s_list[[x]]$mu_fit3 - mu_list[[x]]$mu_fit3) * (alpha_s_list[[x]]$alpha3 - alpha_list[[x]]$alpha3))
        bias2 <- mean((mu_s_list[[x]]$mu_fit2 - mu_list[[x]]$mu_fit2) * (alpha_s_list[[x]]$alpha2 - alpha_list[[x]]$alpha2))
        bias1 <- mean((mu_s_list[[x]]$mu_fit1 - mu_list[[x]]$mu_fit1) * (alpha_s_list[[x]]$alpha1 - alpha_list[[x]]$alpha1))
        bias <- bias3 + bias2 + bias1
    })
    
    nu2_s <- map(c(1:2), \(i) map(alpha_s_list[[i]], \(x) mean(x^2)))
    # nu2 <- map(c(1:2), \(i) map(alpha_list[[i]], \(x) mean(x^2)))
    
    one_minus_R2_alpha <- map(1:2, \(i) map(1:3, \(x) mean((alpha_l_list[[i]][[x]] - alpha_s_list[[i]][[x]])^2) / mean(alpha_s_list[[i]][[x]]^2) ))
    
    alpha_l_list_diff <- reduce(alpha_l_list, `-`)
    alpha_s_list_diff <- reduce(alpha_s_list, `-`)
    map(1:3, \(x) mean((alpha_l_list_diff[[x]] - alpha_s_list_diff[[x]])^2) / mean(alpha_s_list_diff[[x]]^2) )
    
    
    resid_3 <- theta.3l(A, M2, M1, X) - theta.3s(A, M2, M1, Xs)
    resid_2 <- theta.2l.3s(A, M1, X) - theta.2s(A, M1, Xs)
    resid_1 <- theta.1l.2s(A, X) - theta.1s(A, Xs)
    
    Ypseudo_list <- data.frame(Ypseudo4 = Y, Ypseudo3 = Ypseudo3, Ypseudo2 = Ypseudo2)
    map(1:3, \(x=1) {
        resid_x_num <- mu_long_short[[glue("mu_fit{x}")]] - mu_s_list[[1]][[glue("mu_fit{x}")]]
        resid_x_den <- Ypseudo_list[[glue("Ypseudo{x+1}")]] - mu_s_list[[1]][[glue("mu_fit{x}")]]
        mean(resid_x_num^2) / mean(resid_x_den^2)
    })
    
    
}
