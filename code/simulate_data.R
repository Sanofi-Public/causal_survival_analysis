
############ RCT 
# RCT1:  Random treatment assignment + independent censoring
# RCT2:  Random treatment assignment + dependent censoring (conditional on X 
# and A)
simulate_data_RCT <- function(n, mu = c(1, 1, -1, 1), 
                              sigma = diag(4), 
                              colnames_cov = c("X1", "X2", "X3", "X4"),
                              tau, 
                              coefT0 = 0.01,
                              parsS = c(0.5, 0.5, -0.5, 0.5), 
                              coefC = 0.03,
                              parsC = c(0.7, 0.3, -0.25, -0.1), 
                              parsC_A = c(-0.2), 
                              scenario = "RCT2",
                              mis_specification="none") {
  
  # Generate X from a multivariate normal distribution
  X <- MASS::mvrnorm(n, mu, sigma)
  X <- as.data.frame(X)
  colnames(X) <- colnames_cov
  
  # Treatment variable selection: all X
  X_treatment <- as.matrix(X)
  
  # Propensity score: constant for random assignment
  e <- rep(0.5, n)
  
  # Random treatment assignment
  A <- sapply(e, FUN = function(p) rbinom(1, 1, p))
  
  # Outcome variable selection: all X
  X_outcome <- as.matrix(X)
  
  # Simulate the outcome using the cumulative hazard inversion method
  epsilon <- runif(n, min = 1e-8, max = 1)
  T0 <- -log(epsilon) / (coefT0 * exp(X_outcome %*% parsS))
  
  if (scenario == "RCT1") {
    # Simulate independent censoring time
    epsilon <- runif(n, min = 1e-8, max = 1)
    C <- -log(epsilon) / coefC
  }
  else if (scenario == "RCT2") {
    # Simulate dependent censoring time
    X_censoring <- as.matrix(cbind(X,A))
    parsC <- c(parsC,parsC_A)
    
    epsilon <- runif(n, min = 1e-8, max = 1)
    C <- -log(epsilon) / (coefC * exp(rowSums(X_censoring %*% diag(parsC))))
  }
  # T(1) = T(0) + 10
  T1 <- T0 + 10
  
  # True survival time
  T_true <- A * T1 + (1 - A) * T0
  
  # Observed time
  T_obs <- pmin(T_true, C)
  
  # Status indicator
  status <- as.numeric(T_true <= C)
  censor.status <- as.numeric(T_true > C)
  
  # Restricted survival time
  T_obs_tau <- pmin(T_obs, tau)
  status_tau <- as.numeric((T_obs > tau) | (T_obs <= tau & status == 1))
  
  # Combine all data into a single data frame
  data_target_population <- data.frame(X, tau, A, T0, T1, C, T_obs, T_obs_tau, 
                                       status, censor.status, status_tau, e)
  
  return(data_target_population)
}


# Obs1:  Treatment assignment dependent on X + independent censoring
# Obs2:  Treatment assignment dependent on X + dependent censoring (conditional 
# on X)

# Function to simulate observational data for two scenarios: Obs1 and Obs2
simulate_data_obs <- function(n, 
                              mu = c(1, 1, -1, 1), 
                              sigma = diag(4), 
                              colnames_cov = c("X1", "X2", "X3", "X4"),
                              tau,
                              coefT0 = 0.01, 
                              parsS = c(0.5, 0.5, -0.5, 0.5),
                              parsA = c(-1, -1, -2.5, -1), 
                              parsC_A = c(0), 
                              coefC = 0.03,
                              parsC = c(0.7, 0.3, -0.25, -0.1), 
                              scenario = "Obs2") {
  
  # Generate covariates X from a multivariate normal distribution
  X <- mvrnorm(n, mu, sigma)
  X <- as.data.frame(X)
  colnames(X) <- colnames_cov
  
  # Propensity score model based on X
  e <- rowSums(as.matrix(X) %*% diag(parsA))
  e <- plogis(e)  # Transform to probability scale
  
  # Treatment assignment based on the propensity score
  A <- sapply(e, FUN = function(p) rbinom(n = 1, size = 1, prob = p))
  
  # Outcome model based on X
  X_outcome <- as.matrix(X)
  epsilon <- runif(n, min = 0.00000001, max = 1)
  T0 <- -log(epsilon) / (coefT0 * exp(X_outcome %*% parsS))
  
  # Define treatment effect (shift in survival time due to treatment)
  T1 <- T0 + 10
  
  if (scenario == "Obs1") {
    # Scenario 1: Independent censoring
    C <- -log(runif(n, min = 0.00000001, max = 1)) / coefC
    
  } else if (scenario == "Obs2") {
    # Scenario 2: Dependent censoring based on X
    X_censoring <- as.matrix(cbind(X,A))
    parsC <- c(parsC,parsC_A)
    
    C <- -log(runif(n, min = 0.00000001, max = 1)) / 
      (coefC * exp(rowSums(X_censoring %*% diag(parsC))))
    
  } else {
    stop("Invalid scenario. Choose 'Obs1' or 'Obs2'.")
  }
  
  # Determine the true survival time based on treatment
  T_true <- A * T1 + (1 - A) * T0
  
  # Observed time is the minimum of the true survival time and censoring time
  T_obs <- pmin(T_true, C)
  
  # Status indicator: 1 if the event (death) occurred, 0 if censored
  status <- as.numeric(T_true <= C)
  
  # Restricted survival time (censored at tau)
  T_obs_tau <- pmin(T_obs, tau)
  status_tau <- as.numeric((T_obs > tau) | (T_obs <= tau & status == 1))
  
  # Compile the simulated data into a data frame
  DATA_target_population <- data.frame(X, tau, A, T0, T1, C, T_obs, T_obs_tau, 
                                       status, status_tau, e)
  
  return(DATA_target_population)
}


# DGP for misspecification 
simulate_data_mis <- function(n, 
                              mu = c(0.5, 0.5, 0.7, 0.5),
                              sigma =  matrix(c(1, 0, 0, 0, 
                                                0, 1, 0, 0, 
                                                0, 0, 1, 0,
                                                0, 0, 0, 1), 
                                              nrow = 4, byrow = TRUE),
                              colnames_cov = c("X1", "X2", "X3", "X4"),
                              parsA =  c(0.05, -0.1, 0.5, -0.1),
                              tau){
  
  # Generate X from a multivariate normal distribution
  X <- MASS::mvrnorm(n, mu, sigma)
  X <- as.data.frame(X)
  colnames(X) <- colnames_cov
  
  # Treatment variable selection: all X
  X_treatment <- as.matrix(X)
  
  # Propensity score model based on X
  e <- parsA[1]*X_treatment[, "X1"]^2 + parsA[2]*X_treatment[, "X2"]^2 + 
    parsA[3]*X_treatment[, "X3"]^2 + parsA[4]*X_treatment[, "X4"]^2-
    X_treatment[, "X1"]*X_treatment[, "X2"] +
    X_treatment[, "X1"]*X_treatment[, "X4"]
  
  # Logistic regression
  e <- plogis(e)
  
  # Treatment assignment based on the propensity score
  A <- sapply(e, FUN = function(p) rbinom(n = 1, size = 1, prob = p))
  
  # Outcome variable selection: all X
  X_outcome <- as.matrix(X)
  
  lambda <- exp(0.2*X[,1]^2 + 0.3*X[,2]^2 + 0.1*X[,3]^2 + 0.1*X[,4]^2 + 
                  X[,1] * X[,2] + X[,3] * X[,4])
  # Simulate the outcome using the cumulative hazard inversion method
  epsilon <- runif(n, min = 1e-8, max = 1)
  T0 <- -log(epsilon) / lambda
  
  # Simulate independent censoring time
  censoring_lambda <- exp(0.05*X[,1]^2 + 0.05*X[,2]^2-0.1*X[,3]^2 + 0.1*X[,4]^2 + 
                            X[,3] * X[,1] - X[,2]*X[,4])
  epsilon <- runif(n, min = 1e-8, max = 1)
  C <- -log(epsilon) / censoring_lambda
  
  
  # T(1) = T(0) + 1
  T1 <- T0 + 1
  
  # True survival time
  T_true <- A * T1 + (1 - A) * T0
  
  # Observed time
  T_obs <- pmin(T_true, C)
  
  # Status indicator
  status <- as.numeric(T_true <= C)
  censor.status <- as.numeric(T_true > C)
  
  # Restricted survival time
  T_obs_tau <- pmin(T_obs, tau)
  status_tau <- as.numeric((T_obs > tau) | (T_obs <= tau & status == 1))
  # Compile the simulated data into a data frame
  DATA_target_population <- data.frame(X, tau, A, T0, T1, C, T_obs, T_obs_tau, 
                                       status, status_tau, censor.status, e)
  
  return(DATA_target_population)
}
