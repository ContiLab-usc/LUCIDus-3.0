

Mstep_GtoX <- function(G, r, selectG, penalty, K, N) {
  nOmics <- length(K)
  dimG <- dim(G)
  # store multinomial logistic regression model with corresponding coefficients
  fit <- vector(mode = "list", length = nOmics)
  Beta <- vector(mode = "list", length = nOmics)

  # if 2 omics layers
  # how to do it, not selected at all for G to be excluded?? !!!!!
  if(nOmics == 2) {
    for(i in 1:nOmics) {
      r_margin <- t(sapply(1:N, function(j) {
        marginSums(lastInd(r,j), margin = i)
      }))
      if(selectG){
        if(dimG < 2) {
          stop("At least 2 exposure variables are needed to conduct variable selection")
        }
        tryLasso <- try(glmnet(as.matrix(G),
                               as.matrix(r_margin),
                               family = "multinomial",
                               lambda = penalty))
        if("try-error" %in% class(tryLasso)){
          breakdown <- TRUE
          print(paste("lasso failed"))
        }
      }else{
        invisible(capture.output(temp_fit <- nnet::multinom(r_margin ~ G)))
        fit[[i]] <- temp_fit
        Beta[[i]] <- coef(temp_fit)
      }
    }
  }

  # if 3 omics layers
  if(nOmics == 3 | nOmics == 4 |nOmics == 5) {
    for(i in 1:nOmics) {
      r_margin <- t(sapply(1:N, function(j) {
        marginSums(lastInd(r,j), margin = i)
      }))
      invisible(capture.output(temp_fit <- nnet::multinom(r_margin ~ G)))
      fit[[i]] <- temp_fit
      Beta[[i]] <- coef(temp_fit)
    }
  }


  for (i in 1:nOmics){
    colnames(Beta[[i]])[2:length(colnames(Beta[[i]]))] = colnames(G)
  }

  return(list(fit = fit,
              Beta = Beta))
}


Mstep_XtoZ <- function(Z, r, selectZ, penalty.mu, penalty.cov, K, modelNames, N, na_pattern, mu){
  nOmics <- length(K)
  # store GMM model with corresponding model
  Mu <- vector(mode = "list", length = nOmics)
  Sigma <- vector(mode = "list", length = nOmics)
  
  for(i in 1:nOmics) {
    r_margin <- t(sapply(1:N, function(j) {
      marginSums(lastInd(r,j), margin = i)
    }))
    r_margin <- round(r_margin, digits = 8)
    
    if (selectZ){
      temp_result = Mstep_Z_parallel(Z = Z[[i]], 
                                     r = r_margin, 
                                     selectZ = selectZ, penalty.mu = penalty.mu, penalty.cov = penalty.cov,
                                     model.name = modelNames[i], K = K[i], ind.na = na_pattern[[i]]$indicator_na, mu = mu[[i]])
      
      
      Mu[[i]] <- temp_result$mu
      Sigma[[i]] <- temp_result$sigma
    }else{
      temp_fit <- mstep(data = Z[[i]][na_pattern[[i]]$indicator_na != 3, ],
                        G = K[i],
                        z = r_margin[na_pattern[[i]]$indicator_na != 3, ],
                        modelName = modelNames[i])
      Mu[[i]] <- temp_fit$parameters$mean
      Sigma[[i]] <- temp_fit$parameters$variance$sigma
    }
  }
  
  
  
  
  
  
  return(list(Mu = Mu,
              Sigma = Sigma))
}






Mstep_XtoY <- function(Y, r, K, N, family, CoY) {


  # if 2 omics layers
  if(length(K) == 2) {
    r_matrix <- t(sapply(1:N, function(i) {
      c(rowSums(lastInd(r,i)), colSums(lastInd(r,i)))
    }))
    r_fit <- r_matrix[, -c(1, K[1] + 1)]

    if(family == "gaussian") {
      if(is.null(CoY)) {
        fit <- lm(Y ~ r_fit)
        mu <- as.numeric(coef(fit))
        sd <- sd(resid(fit))
      }else{
        Set0 <- as.data.frame(cbind(Y, r_fit))
        Set0 <- cbind(Set0, CoY)
        colnames(Set0) <- c("Y", paste0("LC", 1:ncol(r_fit)), colnames(CoY))
        fit <- glm(as.formula(paste("Y~", paste(colnames(Set0)[-1], collapse = "+"))), data = Set0, family = gaussian)
        beta_f <- coef(fit)
        mu <- as.numeric(beta_f)
        sd <- sd(resid(fit))
      }}

    if(family == "binomial") {
      if(is.null(CoY)) {
        fit <- glm(Y ~ r_fit, family = "binomial")
        mu <- as.numeric(coef(fit))
      }else{
        Set0 <- as.data.frame(cbind(Y, r_fit, CoY))
        colnames(Set0) <- c("Y", paste0("LC", 1:ncol(r_fit)), colnames(CoY))
        fit <- glm(as.formula(paste("Y~", paste(colnames(Set0)[-1], collapse = "+"))), data = Set0, family ="binomial")
        mu <- as.numeric(coef(fit))
      }

      mu_array <- vec_to_array(K = K, mu = mu)
      p <- exp(mu_array) / sum(exp(mu_array))
      sd <- NULL
      fit <- fit
      mu <- p
      sd <- NULL
    }

    if(any(is.na(mu))) {
      na_index <- which(is.na(mu))
      if(na_index <= K[1]) {
        stop("no cluster strucutre is defined for Z1")
      } else{
        stop("no cluster structure is defined for Z2")
      }
    }
  }


  # if 3 omics layers
  if(length(K) == 3) {
    r_matrix <- t(sapply(1:N, function(i) {
      c(marginSums(lastInd(r,i), margin = 1),
        marginSums(lastInd(r,i), margin = 2),
        marginSums(lastInd(r,i), margin = 3))
    }))
    r_fit <- r_matrix[, -c(1, K[1] + 1, K[1] + K[2] + 1)]
    #We actually use ip as the input for the fit of gamma
    #but we would remove reference ip columnds
    #so the coefficient we get for gamma are as if we fit categorical predictors
    #but indeed we use ips for more information
    if(family == "gaussian") {
      if(is.null(CoY)) {
        fit <- lm(Y ~ r_fit)
        mu <- as.numeric(coef(fit))
        sd <- sd(resid(fit))
      }else{
        Set0 <- as.data.frame(cbind(Y, r_fit))
        Set0 <- cbind(Set0, CoY)
        colnames(Set0) <- c("Y", paste0("LC", 1:ncol(r_fit)), colnames(CoY))
        fit <- glm(as.formula(paste("Y~", paste(colnames(Set0)[-1], collapse = "+"))), data = Set0, family = gaussian)

        beta_f <- coef(fit)
        mu <- as.numeric(beta_f)
        sd <- sd(resid(fit))
      }}

    if(family == "binomial") {
      if(is.null(CoY)) {
        fit <- glm(Y ~ r_fit, family = "binomial")
        mu <- as.numeric(coef(fit))
      }else{
        Set0 <- as.data.frame(cbind(Y, r_fit, CoY))
        colnames(Set0) <- c("Y", paste0("LC", 1:ncol(r_fit)), colnames(CoY))
        fit <- glm(as.formula(paste("Y~", paste(colnames(Set0)[-1], collapse = "+"))), data = Set0, family ="binomial")

        mu <- as.numeric(coef(fit))

      }

      mu_array <- vec_to_array(K = K, mu = mu)
      p <- exp(mu_array) / sum(exp(mu_array))
      sd <- NULL
      fit <- fit
      mu <- p
      sd <- NULL
    }

    if(any(is.na(mu))) {
      na_index <- which(is.na(mu))
      if(na_index <= K[1]) {
        stop("no cluster strucutre is defined for Z1")
      } else{
        stop("no cluster structure is defined for Z2")
      }
    }
  }


  # if 4 omics layers
  if(length(K) == 4) {
    r_matrix <- t(sapply(1:N, function(i) {
      c(marginSums(lastInd(r,i), margin = 1),
        marginSums(lastInd(r,i), margin = 2),
        marginSums(lastInd(r,i), margin = 3),
        marginSums(lastInd(r,i), margin = 4))
    }))
    r_fit <- r_matrix[, -c(1, K[1] + 1, K[1] + K[2] + 1, K[1] + K[2] + K[3] + 1)]

    if(family == "gaussian") {
      if(is.null(CoY)) {
        fit <- lm(Y ~ r_fit)
        mu <- as.numeric(coef(fit))
        sd <- sd(resid(fit))
      }else{
        fit <- lm(Y ~ r_fit + CoY)
        beta_f <- coef(fit)
        mu <- as.numeric(beta_f)
        sd <- sd(resid(fit))
      }}

    if(family == "binomial") {
      if(is.null(CoY)) {
        fit <- glm(Y ~ r_fit, family = "binomial")
        mu <- as.numeric(coef(fit))
      }else{
        fit <- glm(Y ~ r_fit + CoY, family = "binomial")
        mu <- as.numeric(coef(fit))
      }

      mu_array <- vec_to_array(K = K, mu = mu)
      p <- exp(mu_array) / sum(exp(mu_array))
      sd <- NULL
      fit <- fit
      mu <- p
      sd <- NULL
    }

    if(any(is.na(mu))) {
      na_index <- which(is.na(mu))
      if(na_index <= K[1]) {
        stop("no cluster strucutre is defined for Z1")
      } else{
        stop("no cluster structure is defined for Z2")
      }
    }
  }


  # if 5 omics layers
  if(length(K) == 5) {
    r_matrix <- t(sapply(1:N, function(i) {
      c(marginSums(lastInd(r,i), margin = 1),
        marginSums(lastInd(r,i), margin = 2),
        marginSums(lastInd(r,i), margin = 3),
        marginSums(lastInd(r,i), margin = 4),
        marginSums(lastInd(r,i), margin = 5))
    }))
    r_fit <- r_matrix[, -c(1,
                           K[1] + 1,
                           K[1] + K[2] + 1,
                           K[1] + K[2] + K[3] + 1,
                           K[1] + K[2] + K[3] + K[4] + 1)]

    if(family == "gaussian") {
      if(is.null(CoY)) {
        fit <- lm(Y ~ r_fit)
        mu <- as.numeric(coef(fit))
        sd <- sd(resid(fit))
      }else{
        fit <- lm(Y ~ r_fit + CoY)
        beta_f <- summary(fit)$coefficients[, 1]
        mu <- as.numeric(beta_f)
        sd <- sd(resid(fit))
      }}

    if(family == "binomial") {
      if(is.null(CoY)) {
        fit <- glm(Y ~ r_fit, family = "binomial")
        mu <- as.numeric(coef(fit))
      }else{
        fit <- glm(Y ~ r_fit + CoY, family = "binomial")
        mu <- as.numeric(coef(fit))
      }

      mu_array <- vec_to_array(K = K, mu = mu)
      p <- exp(mu_array) / sum(exp(mu_array))
      sd <- NULL
      fit <- fit
      mu <- p
      sd <- NULL
    }

    if(any(is.na(mu))) {
      na_index <- which(is.na(mu))
      if(na_index <= K[1]) {
        stop("no cluster strucutre is defined for Z1")
      } else{
        stop("no cluster structure is defined for Z2")
      }
    }
  }



  Delta <- list(mu = mu,
                sd = sd,
                K = K)
  return(list(fit = fit,
              Gamma = Delta))
}


#modified Mstep_Z_parallel with LASSO penalty
Mstep_Z_parallel <- function(Z, r, selectZ, penalty.mu, penalty.cov,
                             model.name, K, ind.na, mu) {
  mu = t(mu)
  dz <- Z[ind.na != 3, ]
  #cat("dz",dz)
  dr <- r[ind.na != 3, ]
  #cat("dr",dr)
  Q <- ncol(Z)
  #cat("Q",Q)
  new_sigma <- array(rep(0, Q^2 * K), dim = c(Q, Q, K))
  #cat("new_sigma",new_sigma)
  new_mu <- matrix(rep(0, Q * K), nrow = K)
  #cat("new_mu",new_mu)
  
  k <- 1
  while(k <= K){
    #estimate E(S_k) to be used by glasso
    Z_mu <- t(t(dz) - mu[k, ])
    
    E_S <- (matrix(colSums(dr[, k] * t(apply(Z_mu, 1, function(x) return(x %*% t(x))))), Q, Q)) / sum(dr[, k])
    
    #use glasso and E(S_k) to estimate new_sigma and new_sigma_inv
    l_cov <- try(glasso(E_S, penalty.cov))
    if("try-error" %in% class(l_cov)){
      print(paste("glasso failed, restart lucid \n"))
      break
    }
    else{
      new_sigma[, , k] <- l_cov$w
      # function to calculate mean
      new_mu[k, ] <- est_mu(j = k,
                            rho = penalty.mu,
                            z = dz,
                            r = dr,
                            mu = mu[k, ],
                            wi = l_cov$wi)
    }
    k <- k + 1
  }
  if("try-error" %in% class(l_cov)){
    return(structure(list(mu = NULL,
                          sigma = NULL)))
  }else{
    return(structure(list(mu = t(new_mu),
                          sigma = new_sigma)))
  }
}


