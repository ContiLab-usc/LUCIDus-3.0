
####################utility functions for the EM algorithm of LUCID in parallel####################
#' Log-Sum-Exp Trick
#'
#' Computes the log-sum-exp trick for a vector.
#'
#' @param vec A numeric vector.
#' @return The result of the log-sum-exp trick.
#'
#'
LogSumExp <- function(vec) {
  max_vec <- max(vec)
  trick <- max_vec + log(sum(exp(vec - max_vec)))
  return(trick)
}

#' Get the Last 'n' Elements of an Array
#'
#' Retrieves the last 'n' elements from a vector or array.
#'
#' @param x A vector or array.
#' @param n The number of elements to retrieve.
#' @return A vector or array containing the last 'n' elements.
#'

lastInd <- function(x, n){
  d <- dim(x)
  d.new <- d[-length(d)]
  block.size <- prod(d.new)
  res <- x[(block.size * (n - 1) + 1):(block.size * n)]
  array(res, dim = d.new)
}

#' Check Validity of 'K'
#'
#' Checks if the elements in 'K' are integers and greater than or equal to 2.
#'
#' @param K A vector of integers.
#'

check_K <- function(K) {
  for(x in K) {
    if(x != as.integer(x)) {
      stop("K should be a vector of integer")
    }
  }
  if(min(K) < 2) {
    stop("each element in K should be greater or equal than 2")
  }
}


#' Initialize Beta
#'
#' Initializes the Beta parameter for the EM algorithm.
#'
#' @param K A vector specifying the number of clusters for each omics layer.
#' @param nG The number of exposures.
#'

initialize_Beta <- function(K, nG) {
  nOmics <- length(K)
  Beta <- vector(mode = "list", length = nOmics)
  for(i in 1:nOmics) {
    # row represent cluster, column represents variable
    Beta[[i]] <- matrix(runif((nG + 1) * (K[i] - 1), min = -1, max = 1),
                        nrow = K[i] - 1)
  }
  return(Beta)
}

#' Initialize Mu
#'
#' Initializes the Mu parameter for the EM algorithm.
#'
#' @param K A vector specifying the number of clusters for each omics layer.
#' @param nZ A vector specifying the number of variables for each omics layer.
#'

initialize_Mu <- function(K, nZ) {
  nOmics <- length(K)
  Mu <- vector(mode = "list", length = nOmics)
  for(i in 1:nOmics) {
    # column represents cluster, row represents variable
    Mu[[i]] <- matrix(runif(K[i] * nZ[i], min = -1, max = 1),
                      nrow = nZ[i])
  }
  return(Mu)
}


#' Initialize Sigma
#'
#' Initializes the Sigma parameter for the EM algorithm.
#'
#' @param K A vector specifying the number of clusters for each omics layer.
#' @param nZ A vector specifying the number of variables for each omics layer.
#'

initialize_Sigma <- function(K, nZ) {
  nOmics <- length(K)
  Sigma <- vector(mode = "list", length = nOmics)
  for(i in 1:nOmics) {
    # each element is an array of nZ[i] x nZ[i] x K[i]
    Sigma[[i]] <- array(0, dim = c(nZ[i], nZ[i], K[i]))
    for(j in 1:K[i]) {
      Sigma[[i]][, , j] <- diag(nZ[i])
    }
  }
  return(Sigma)
}



#' Initialize Mu and Sigma
#'
#' Initializes the Mu and Sigma parameters for the EM algorithm.
#'
#' @param K A vector specifying the number of clusters for each omics layer.
#' @param Z A list of matrices representing the data for each omics layer.
#' @param modelNames A list of model names for each omics layer.
#' @param na_pattern A list of NA patterns for each omics layer.
#'

initialize_Mu_Sigma <- function(K, Z, modelNames, na_pattern) {
  nOmics <- length(K)
  Mu <- vector(mode = "list", length = nOmics)
  Sigma <- vector(mode = "list", length = nOmics)
  z <- vector(mode = "list", length = nOmics)
  for(i in 1:nOmics) {
    # don't initial beta for na_pattern[[i]]$indicator_na == 3
    temp_fit <- Mclust(data = Z[[i]][na_pattern[[i]]$indicator_na != 3, ],
                       G = K[i],
                       modelNames = modelNames[i])
    Mu[[i]] <- temp_fit$parameters$mean
    Sigma[[i]] <- temp_fit$parameters$variance$sigma
    z[[i]] <- temp_fit$z
  }
  return(list(Mu = Mu,
              Sigma = Sigma,
              z = z))
}

#' Initialize Delta
#'
#' Initializes the Delta parameter for the EM algorithm.
#'
#' @param K A vector specifying the number of clusters for each omics layer.
#' @param CoY A data frame for covariates (optional).
#' @param family The distribution family ("gaussian" or "binomial").
#' @param z A list of matrices representing the clustering assignments for each omics layer.
#' @param Y The outcome variable.
#'

initialize_Delta <- function(K, CoY, family = c("gaussian", "binomial"),
                             z, Y) {
  family <- match.arg(family)
  if(family == "gaussian") {

    # if 2 omics layers
    if(length(K) == 2) {
      r_matrix <- cbind(z[[1]], z[[2]])
      r_fit <- r_matrix[, -c(1, K[1] + 1)]

      if(is.null(CoY)) {
        fit <- lm(Y ~ r_fit)
        mu <- as.numeric(coef(fit))
        sd <- sd(resid(fit))
        x <- list(mu = mu,
                  sd = sd,
                  K = K)
      }else{
        Set0 <- as.data.frame(cbind(Y, r_fit))
        Set0 <- cbind(Set0, CoY)
        colnames(Set0) <- c("Y", paste0("LC", 1:ncol(r_fit)), colnames(CoY))
        fit <- glm(as.formula(paste("Y~", paste(colnames(Set0)[-1], collapse = "+"))), data = Set0, family = gaussian)

        beta_f <- coef(fit)
        mu <- as.numeric(beta_f)
        sd <- sd(resid(fit))
        x <- list(mu = mu,
                  sd = sd,
                  K = K)
      }}

    # if 3 omics layers
    if(length(K) == 3) {
      r_matrix <- cbind(z[[1]], z[[2]], z[[3]])
      r_fit <- r_matrix[, -c(1, K[1] + 1, K[1] + K[2] + 1)]

      if(is.null(CoY)) {
        fit <- lm(Y ~ r_fit)
        mu <- as.numeric(coef(fit))
        sd <- sd(resid(fit))
        x <- list(mu = mu,
                  sd = sd,
                  K = K)
      }else{
        Set0 <- as.data.frame(cbind(Y, r_fit))
        Set0 <- cbind(Set0, CoY)
        colnames(Set0) <- c("Y", paste0("LC", 1:ncol(r_fit)), colnames(CoY))
        fit <- glm(as.formula(paste("Y~", paste(colnames(Set0)[-1], collapse = "+"))), data = Set0, family = gaussian)

        beta_f <- coef(fit)
        mu <- as.numeric(beta_f)
        sd <- sd(resid(fit))
        x <- list(mu = mu,
                  sd = sd,
                  K = K)
      }}

    # if 4 omics layers
    if(length(K) == 4) {
      r_matrix <- cbind(z[[1]], z[[2]], z[[3]], z[[4]])
      r_fit <- r_matrix[, -c(1, K[1] + 1, K[1] + K[2] + 1, K[1] + K[2] + K[3] + 1)]

      if(is.null(CoY)) {
        fit <- lm(Y ~ r_fit)
        mu <- as.numeric(coef(fit))
        sd <- sd(resid(fit))
        x <- list(mu = mu,
                  sd = sd,
                  K = K)
      }else{
        fit <- lm(Y ~ r_fit + CoY)
        beta_f <- summary(fit)$coefficients[, 1]
        mu <- as.numeric(beta_f)
        sd <- sd(resid(fit))
        x <- list(mu = mu,
                  sd = sd,
                  K = K)
      }}

    # if 5 omics layers
    if(length(K) == 5) {
      r_matrix <- cbind(z[[1]], z[[2]], z[[3]], z[[4]], z[[5]])
      r_fit <- r_matrix[, -c(1, K[1] + 1, K[1] + K[2] + 1, K[1] + K[2] + K[3] + 1, K[1] + K[2] + K[3] + K[4] + 1)]

      if(is.null(CoY)) {
        fit <- lm(Y ~ r_fit)
        mu <- as.numeric(coef(fit))
        sd <- sd(resid(fit))
        x <- list(mu = mu,
                  sd = sd,
                  K = K)
      }else{
        fit <- lm(Y ~ r_fit + CoY)
        beta_f <- summary(fit)$coefficients[, 1]
        mu <- as.numeric(beta_f)
        sd <- sd(resid(fit))
        x <- list(mu = mu,
                  sd = sd,
                  K = K)
      }}

  }


  if(family == "binomial") {


    # if 2 omics layers
    if(length(K) == 2) {
      r_matrix <- cbind(z[[1]], z[[2]])
      r_fit <- r_matrix[, -c(1, K[1] + 1)]

      if(is.null(CoY)) {
        fit <- glm(Y ~ r_fit, family = "binomial")
        b <- as.numeric(coef(fit))
      }else{

        Set0 <- as.data.frame(cbind(Y, r_fit, CoY))
        colnames(Set0) <- c("Y", paste0("LC", 1:ncol(r_fit)), colnames(CoY))
        fit <- glm(as.formula(paste("Y~", paste(colnames(Set0)[-1], collapse = "+"))), data = Set0, family ="binomial")
        b <- as.numeric(coef(fit))

      }

      b_array <- vec_to_array(K = K, mu = b)
      p <- 1 / (1 + exp(-b_array))
      x <- list(mu = p,
                K = K)
    }


    # if 3 omics layers
    if(length(K) == 3) {
      r_matrix <- cbind(z[[1]], z[[2]], z[[3]])
      r_fit <- r_matrix[, -c(1, K[1] + 1, K[1] + K[2] + 1)]

      if(is.null(CoY)) {
        fit <- glm(Y ~ r_fit, family = "binomial")
        b <- as.numeric(coef(fit))
      }else{
        Set0 <- as.data.frame(cbind(Y, r_fit, CoY))
        colnames(Set0) <- c("Y", paste0("LC", 1:ncol(r_fit)), colnames(CoY))
        fit <- glm(as.formula(paste("Y~", paste(colnames(Set0)[-1], collapse = "+"))), data = Set0, family ="binomial")
        b <- as.numeric(coef(fit))
      }

      b_array <- vec_to_array(K = K, mu = b)
      p <- 1 / (1 + exp(-b_array))
      x <- list(mu = p,
                K = K)
    }

    # if 4 omics layers
    if(length(K) == 4) {
      r_matrix <- cbind(z[[1]], z[[2]], z[[3]], z[[4]])
      r_fit <- r_matrix[, -c(1, K[1] + 1, K[1] + K[2] + 1, K[1] + K[2] + K[3] + 1)]

      if(is.null(CoY)) {
        fit <- glm(Y ~ r_fit, family = "binomial")
        b <- as.numeric(coef(fit))
      }else{
        fit <- glm(Y ~ r_fit + CoY, family = "binomial")
        b <- as.numeric(coef(fit))
      }

      b_array <- vec_to_array(K = K, mu = b)
      p <- 1 / (1 + exp(-b_array))
      x <- list(mu = p,
                K = K)
    }

    # if 5 omics layers
    if(length(K) == 5) {
      r_matrix <- cbind(z[[1]], z[[2]], z[[3]], z[[4]], z[[4]])
      r_fit <- r_matrix[, -c(1, K[1] + 1, K[1] + K[2] + 1, K[1] + K[2] + K[3] + 1, K[1] + K[2] + K[3] + K[4] + 1)]

      if(is.null(CoY)) {
        fit <- glm(Y ~ r_fit, family = "binomial")
        b <- as.numeric(coef(fit))
      }else{
        fit <- glm(Y ~ r_fit + CoY, family = "binomial")
        b <- as.numeric(coef(fit))
      }

      b_array <- vec_to_array(K = K, mu = b)
      p <- 1 / (1 + exp(-b_array))
      x <- list(mu = p,
                K = K)
    }


  }
  return(x)
}


#' Indicator Function
#'
#' Computes an indicator function.
#'
#' @param x A numeric value.
#' @return 1 if x > 1, otherwise 0.
#'

indicator <- function(x) {
  m <- 0
  if(x > 1) {
    m <- 1
  }
  return(m)
}

#' Transform Mu to an Array
#'
#' Transforms a vector of Mu values to an array based on the number of clusters.
#'
#' @param K A vector specifying the number of clusters for each omics layer.
#' @param mu A vector of Mu values.
#' @return A multi-dimensional array of Mu values.
#'

vec_to_array <- function(K, mu) {
  res <- array(data = rep(0, prod(K)),
               dim = K)

  # if nK = 2, transform mu to a matrix
  if(length(K) == 2) {
    for(i in 1:K[1]) {
      for(j in 1:K[2]) {
        res[i, j] <- mu[1] + indicator(i) * mu[i] + indicator(j) * mu[K[1] + j - 1]
      }
    }
  }

  # if nK = 3, transform mu to a 3d array
  if(length(K) == 3) {
    for(i in 1:K[1]) {
      for(j in 1:K[2]) {
        for (k in 1:K[3]) {
          res[i, j, k] <- mu[1] + indicator(i) * mu[i] + indicator(j) * mu[K[1] + j - 1] + indicator(k) * mu[K[1] + K[2] + k - 2]
        }
      }
    }
  }

  # if nK = 4, transform mu to a 4d array
  if(length(K) == 4) {
    for(i in 1:K[1]) {
      for(j in 1:K[2]) {
        for (k in 1:K[3]) {
          for(l in 1:K[4]) {
            res[i, j, k, l] <- mu[1] + indicator(i) * mu[i] + indicator(j) * mu[K[1] + j - 1] + indicator(k) * mu[K[1] + K[2] + k - 2] + indicator(l) * mu[K[1] + K[2] + K[3] + l - 3]
          }
        }
      }
    }
  }

  # if nK = 5, transform mu to a 5d array
  if(length(K) == 5) {
    for(i in 1:K[1]) {
      for(j in 1:K[2]) {
        for (k in 1:K[3]) {
          for(l in 1:K[4]) {
            for(m in 1:K[5]) {
              res[i, j, k, l, m] <- mu[1] + indicator(i) * mu[i]
              + indicator(j) * mu[K[1] + j - 1]
              + indicator(k) * mu[K[1] + K[2] + k - 2] +
                indicator(l) * mu[K[1] + K[2] + K[3] + l - 3] +
                indicator(m) * mu[K[1] + K[2] + K[3] + K[4] + m - 4]
            }
          }
        }
      }
    }
  }


  return(res)
}
