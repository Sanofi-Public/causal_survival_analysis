# Naive estimator
Naive <- function(data, tau) {
  # Remove censored observations
  data<- data[data$status == 1, ]
  # Compute the restricted survival time
  data$T_obs_tau <- ifelse(data$T_obs >= tau, tau, data$T_obs)
  
  # Compute the difference of the restricted mean survival time of treated 
  # and control
  mean_naive <- mean(data$T_obs_tau[data$A == 1]) - 
    mean(data$T_obs_tau[data$A == 0])
  
  return(mean_naive)
}


# Kaplan-Meier estimator handmade implementation
# The database 'data' must be in the same form as that shown in 
# notation (Table 1) and with the same variable names (status, T_obs) 
Kaplan_meier_handmade <- function(data, 
                                  status = data$status, 
                                  T_obs = data$T_obs) {
  # Sort unique observed times
  Y.grid <- sort(unique(T_obs))
  
  # Initialize vectors for number of events, number at risk, and survival 
  # probability
  d <- rep(NA, length(Y.grid))  # Number of events at time Y.grid[i]
  n <- rep(NA, length(Y.grid))  # Number at risk just before time Y.grid[i]
  S <- rep(NA, length(Y.grid))  # Survival probability at time Y.grid[i]
  
  # Loop over each unique observed time
  for (i in 1:length(Y.grid)) {
    d[i] <- sum(T_obs == Y.grid[i] & status == 1, na.rm = TRUE)  # Count events
    n[i] <- sum(T_obs >= Y.grid[i])  # Count at risk
    
    # Calculate survival probability
    S[i] <- cumprod(1 - d / n)[i]
  }
  
  # Create a dataframe with the results
  df <- data.frame(d = d, n = n, S = S, T = Y.grid)
  
  return(df)
}


# Function to calculate RMST (Restricted Mean Survival Time):
# Two possibilities for computing RMST: 
# - in using directly S_A1 and S_A0 (survival function of treated and control)
# - in using the dataframe and the function computes the survival functions
RMST_1 <- function(data = NULL, A1 = 1, A0 = 0, tau, S_A1 = NULL, S_A0 = NULL) {
  if (is.null(S_A1) & is.null(S_A0)) {
    # Subset data for treatment groups
    data1 <- data[data$A == A1,]
    data0 <- data[data$A == A0,]
    
    # Calculate Kaplan-Meier survival estimates
    S_A1 <- Kaplan_meier_handmade(data1, status = data1$status, 
                                  T_obs = data1$T_obs)
    S_A0 <- Kaplan_meier_handmade(data0, status = data0$status, 
                                  T_obs = data0$T_obs)
    
    # Restrict observations to those less than or equal to tau
    Y.grid1 <- data1$T_obs[data1$T_obs <= tau]
    Y.grid0 <- data0$T_obs[data0$T_obs <= tau]
  } else {
    # Restrict observations to those less than or equal to tau
    Y.grid1 <- S_A1$T[S_A1$T <= tau]
    Y.grid0 <- S_A0$T[S_A0$T <= tau]
  }
  
  # Filter survival estimates to restricted observations
  S_A1 <- S_A1 %>%
    dplyr::filter(T %in% Y.grid1)
  S_A0 <- S_A0 %>%
    dplyr::filter(T %in% Y.grid0)
  
  # Check if there is any event at tau for S_A1
  if (!any(S_A1$T == tau)) {
    new_row <- tibble(T = tau, S = S_A1$S[nrow(S_A1)])
    S_A1 <- dplyr::bind_rows(S_A1, new_row)
  }
  
  # Check if there is any event at tau for S_A0
  if (!any(S_A0$T == tau)) {
    new_row <- tibble(T = tau, S = S_A0$S[nrow(S_A0)])
    S_A0 <- dplyr::bind_rows(S_A0, new_row)
  }
  
  # Calculate integrals from 0 to tau of survival probabilities
  intA1 <- integral_rectangles(S_A1$T, S_A1$S)
  intA0 <- integral_rectangles(S_A0$T, S_A0$S)
  RMST1 <- intA1 - intA0
  
  return(list(RMST=RMST1, intA1=intA1,intA0=intA0))
}


# Alternative code to estimate Kaplan-Meier estimator with survival package
# instead of handmade KM
RMST_alternative <- function(data, A1 = 1, A0 = 0, tau){
  # Estimate Kaplan-Meier estimator with survfit function on data subset
  # Groupe A = 0
  fit0 <- survfit(Surv(T_obs, status) ~ 1, data = data[data$A == A0,]) 
  # Groupe A = 1
  fit1 <- survfit(Surv(T_obs, status) ~ 1, data = data[data$A == A1,])  
  
  # Estimate the RMST with rmean
  summary_fit0 <- summary(fit0, rmean = tau)  # RMST for A = 0
  summary_fit1 <- summary(fit1, rmean = tau)  # RMST for A = 1
  
  # Extract the RMST from the summary objects
  rmst0 <- summary_fit0$table["rmean"][[1]]
  rmst1 <- summary_fit1$table["rmean"][[1]]
  
  # Compute the difference of RMST between the two groups
  difference_rmst <- rmst1 - rmst0
  return(difference_rmst)
}

# Kaplan-Meier adjusted
# Times of event 
# Failures:  1 if event, 0 if censored
# Variable:  1 if treated, 0 if control
# Weights:  Weight of the individual
adjusted.KM <- function(times, failures, variable, weights = NULL) {
  # Sanity checks
  if (sum(times < 0) > 0) {
    stop("Error: times must be positive")
  }
  if (!is.null(weights) && sum(weights < 0, na.rm = TRUE) > 0) {
    stop("Error: weights must be superior to 0")
  }
  if (sum(failures != 0 & failures != 1) > 0) {
    stop("Error: failures must be a vector of 0 or 1")
  }
  # If 'weights' is NULL, initialize 'w' with ones of the same length as 'times', 
  # otherwise use 'weights'
  w <- if (is.null(weights)) rep(1, length(times)) else weights
  
  # Create a DataFrame 'data' with columns t (times), f (failures), 
  # v (stratification variable: often treatment variable), and w (weights)
  data <- data.frame(t = times, f = failures, v = variable, w = w)
  
  # Remove rows from the DataFrame where the stratification variable is NA
  data <- data[!is.na(data$v),]
  
  # Initialize an empty DataFrame to store the Kaplan-Meier results
  table_KM <- data.frame(times = NULL, n.risk = NULL, n.event = NULL, 
                         survival = NULL, variable = NULL)
  
  # Loop over each unique value of the stratification variable
  for (i in unique(variable)) {
    # Subset the data for the current stratification variable value
    d <- data[data$v == i,]
    
    # Create a sorted vector of unique event times, including time 0 and the 
    # maximum time
    tj <- c(0, sort(unique(d$t[d$f == 1])), max(d$t))
    
    # Calculate the number of events at each time point
    dj <- sapply(tj, function(x) {
      sum(d$w[d$t == x & d$f == 1])
    })
    
    # Calculate the number of individuals at risk at each time point
    nj <- sapply(tj, function(x) {
      sum(d$w[d$t >= x])
    })
    
    # Compute the cumulative product for the survival probabilities
    st <- cumprod((nj - dj) / nj)
    
    # Append the results to the Kaplan-Meier table
    table_KM <- rbind(table_KM, data.frame(T = tj, n = nj, d = dj, 
                                           S = st, variable = i))
  }
  return(table_KM)
}


# IPCW Kaplan-Meier estimator with restricted tau
IPCW_Kaplan_meier <- function(data, tau, 
                              X.names.censoring, 
                              nuisance_censoring = "cox", 
                              n.folds = NULL) {
  
  # Compute of truncated T_obs, status and censored status
  data$T_obs_tau <- ifelse(data$T_obs >= tau, tau, data$T_obs)
  data$censor.status_tau <- 1 - as.numeric((data$T_obs >= tau) | 
                                             (data$T_obs < tau & data$status == 1))
  data$status_tau <- as.numeric((data$T_obs >= tau) | 
                                  (data$T_obs < tau & data$status == 1))
  Y.grid <- sort(unique(data$T_obs_tau))
  
  # Estimate probability of remaining uncensored based on nuisance model 
  S_C_hat <- estimate_survival_function(data = data, X.names = X.names.censoring,
                                        Y.grid = Y.grid, T_obs = "T_obs_tau",
                                        status = "censor.status_tau",
                                        type_of_model = nuisance_censoring,
                                        n.folds = n.folds)
  
  # Select the probability of censoring for each observed T_obs_tau from the 
  # curve
  data$S_C <- S_C_hat$S_hat[cbind(1:nrow(data), match(data$T_obs_tau, Y.grid))]
  
  # Compute IPC weights
  data$weights <- data$status_tau / data$S_C
  
  # Compute the adjusted IPCW Kaplan-Meier
  S <- adjusted.KM(times = data$T_obs, failures = data$status, 
                   variable = data$A, weights = data$weights)
  
  # Compute differenceof RMST between the two groups
  RMST <- RMST_1(S_A1 = S[S$variable == 1,], S_A0 = S[S$variable == 0,], tau = tau)
  
  return(list(RMST = RMST$RMST,
              intA1 = RMST$intA1,
              intA0 = RMST$intA0,
              weights = data$weights))
}

# Alternative code to estimate IPCW Kaplan-Meier, IPTW Kaplan-Meier or 
# IPTW-IPCW Kaplan-Meier estimator with survival package instead of using 
# handmade adjusted.KM function (the weights need to be calculated before).

# Weights0 corresponds to weights of the control and weights1 of treated
Adjusted_Kaplan_meier_alternative <- function(data, A1 = 1, A0 = 0, tau, 
                                              weights0, weights1){
  # Estimate Kaplan-Meier estimator with survfit function on data subset 
  # Groupe A = 0
  fit0 <- survfit(Surv(T_obs, status) ~ 1, data = data[data$A == A0,], weights = weights0)  
  # Groupe A = 1
  fit1 <- survfit(Surv(T_obs, status) ~ 1, data = data[data$A == A1,], weights = weights1)  
  
  # Estimate the RMST with rmean
  summary_fit0 <- summary(fit0, rmean = tau)  # RMST for A = 0
  summary_fit1 <- summary(fit1, rmean = tau)  # RMST for A = 1
  
  # Extract the RMST from the summary objects
  rmst0 <- summary_fit0$table["rmean"][[1]]
  rmst1 <- summary_fit1$table["rmean"][[1]]
  
  # Compute the difference in RMST between the two groups
  difference_rmst <- rmst1 - rmst0
  return(difference_rmst)
}

# Compute the Restricted Mean Survival Time (RMST) difference
BJ <- function(data, tau, X.names.outcome = c("X1", "X2", "X3", "X4"),
               nuisance = "cox", n.folds = NULL) {
  # Truncate observed times at tau
  data$T_obs_tau <- ifelse(data$T_obs >= tau, tau, data$T_obs)
  Y.grid <- sort(unique(data$T_obs_tau))
  
  # Censoring status at tau
  data$status_tau <- as.numeric((data$T_obs >= tau) | 
                                  (data$T_obs < tau & data$status == 1))
  
  # Compute Q_t for all time points
  Q_t <- Q_t_hat(data, tau, X.names.outcome, nuisance, n.folds)
  data$Q_y <- Q_Y(data, tau, Q_t)
  
  # Split data by treatment group
  data_treated <- data %>% dplyr::filter(A == 1)
  data_not_treated <- data %>% dplyr::filter(A == 0)
  
  # Calculate Restricted Survival Time (RST) for each group
  data_treated$RST <- data_treated$status_tau * data_treated$T_obs_tau + 
    (1 - data_treated$status_tau) * data_treated$Q_y
  
  data_not_treated$RST <- data_not_treated$status_tau * data_not_treated$T_obs_tau + 
    (1 - data_not_treated$status_tau) * data_not_treated$Q_y
  
  # Calculate RMST difference between treated and not treated
  RMST <- mean(data_treated$RST) - mean(data_not_treated$RST)
  
  # Return RMST and other relevant metrics
  return(list(
    RMST = RMST, 
    ATE_treated = mean(data_treated$RST), 
    ATE_not_treated = mean(data_not_treated$RST)
  ))
}


# Function to calculate IPTW Kaplan-Meier
IPTW_Kaplan_meier <- function(data, tau, X.names.propensity, 
                              nuisance_propensity = "glm", n.folds = NULL) {
  # Estimate propensity scores
  data$e_hat <- estimate_propensity_score(
    data,
    treatment_covariates = X.names.propensity,
    type_of_model = nuisance_propensity,
    n.folds = n.folds)
  
  # Truncate observed times at tau
  data$T_obs_tau <- pmin(data$T_obs, tau)
  
  # Define censoring status at tau
  data$status_tau <- as.numeric((data$T_obs >= tau) | 
                                  (data$T_obs < tau & data$status == 1))
  
  # Calculate weights
  data$weights <- (data$A) * (1 / data$e_hat) + (1 - data$A) / (1 - data$e_hat)
  
  # Adjusted Kaplan-Meier estimator
  S <- adjusted.KM(
    times = data$T_obs, 
    failures = data$status,
    variable = data$A, 
    weights = data$weights)
  
  # Calculate RMST from the adjusted survival curves
  RMST <- RMST_1(S_A1 = S[S$variable == 1,], 
                 S_A0 = S[S$variable == 0,], 
                 tau = tau)
  
  return(list("intA0" = RMST$intA0, "intA1" = RMST$intA1, "RMST" = RMST$RMST))
}




# Function to estimate the g-formula Two-learner.
g_formula_T_learner <- function(data, 
                                X.names.outcome, 
                                tau, 
                                nuisance_survival = "cox", 
                                n.folds = NULL) {
  # Compute min(T_obs,tau)
  data$T_obs_tau <- ifelse(data$T_obs >= tau, tau, data$T_obs)
  
  # Y.grid is the grid of time points where we want to estimate the 
  # survival function.
  Y.grid <- sort(unique(data$T_obs_tau))
  
  S_hat <- estimate_survival_function(data, X.names.outcome, 
                                      Y.grid, 
                                      type_of_model = nuisance_survival,
                                      T_obs = "T_obs", 
                                      status = "status", 
                                      n.folds = n.folds)
  
  # Compute the area under each survival curve up to max(Y.grid) = tau.
  E_hat1 <- expected_survival(S_hat$S_hat1, Y.grid)
  E_hat0 <- expected_survival(S_hat$S_hat0, Y.grid)
  
  # Calculate the mean difference.
  theta_g_formula <- mean(E_hat1 - E_hat0)
  
  return(theta_g_formula)
}

# Function to estimate the g-formula Single-learner.
g_formula_S_learner <- function(data, 
                                X.names.outcome, 
                                tau, 
                                nuisance_survival = "cox", 
                                n.folds = NULL) {
  # Compute min(T_obs,tau)
  data$T_obs_tau <- ifelse(data$T_obs >= tau, tau, data$T_obs)
  
  # Y.grid is the grid of time points where we want to estimate the 
  # survival function.
  Y.grid <- sort(unique(data$T_obs_tau))
  
  S_hat <- estimate_survival_function(data, X.names.outcome, 
                                      Y.grid, 
                                      type_of_model = nuisance_survival,
                                      learner = "S-learner",
                                      T_obs = "T_obs", 
                                      status = "status", 
                                      n.folds = n.folds)
  
  # Compute the area under each survival curve until max(Y.grid) = tau.
  E_hat1 <- expected_survival(S_hat$S_hat1, Y.grid)
  E_hat0 <- expected_survival(S_hat$S_hat0, Y.grid)
  
  # Calculate the mean difference.
  theta_g_formula <- mean(E_hat1 - E_hat0)
  
  return(theta_g_formula)
}


IPTW_IPCW_Kaplan_meier <- function(data, 
                                   X.names.propensity, 
                                   X.names.censoring, 
                                   tau,
                                   nuisance_propensity = "glm",
                                   nuisance_censoring = "cox",
                                   n.folds = NULL) {
  # Censoring time to tau if observed time exceeds tau
  data$T_obs_tau <- ifelse(data$T_obs >= tau, tau, data$T_obs)
  
  # Create censoring status for tau
  data$censor.status_tau <- 1 - as.numeric((data$T_obs >= tau) | 
                                             (data$T_obs < tau & data$status == 1))
  
  # Create status at tau
  data$status_tau <- as.numeric((data$T_obs >= tau) | 
                                  (data$T_obs < tau & data$status == 1))
  
  # Grid of unique observed times truncated at tau
  Y.grid <- sort(unique(data$T_obs_tau))
  
  # Estimate propensity scores
  data$e_hat <- estimate_propensity_score(data,
                                          treatment_covariates = X.names.propensity,
                                          type_of_model = nuisance_propensity,
                                          n.folds = n.folds)
  
  # Estimate survival function for censoring
  S_C_hat <- estimate_survival_function(data, X.names = X.names.censoring,
                                        Y.grid = Y.grid, T_obs = "T_obs_tau",
                                        status = "censor.status_tau",
                                        type_of_model = nuisance_censoring,
                                        n.folds = n.folds)
  
  # Get estimated survival probabilities for censoring
  data$S_C <- S_C_hat$S_hat[cbind(1:nrow(data), match(data$T_obs_tau, Y.grid))]
  
  # Calculate weights
  data$weights <- data$status_tau / data$S_C * 
    (data$A * (1 / data$e_hat) + 
       (1 - data$A) * (1 / (1 - data$e_hat)))
  
  # Compute adjusted Kaplan-Meier estimator
  S <- adjusted.KM(times = data$T_obs, 
                   failures = data$status, 
                   variable = data$A, 
                   weights = data$weights)
  
  # Compute Restricted Mean Survival Time (RMST)
  RMST <- RMST_1(S_A1 = S[S$variable == 1, ], 
                 S_A0 = S[S$variable == 0, ],
                 tau = tau)
  
  # Return RMST and ATE for treated and not treated groups
  return(list(RMST = RMST$RMST, ATE_treated = RMST$intA1, 
              ATE_not_treated = RMST$intA0))
}



IPTW_BJ <- function(data, 
                    X.names.propensity,
                    X.names.outcome, 
                    tau,
                    nuisance_propensity = "glm",
                    nuisance = "cox",
                    n.folds = NULL) {
  # Minimum of T_obs and tau
  data$T_obs_tau <- ifelse(data$T_obs >= tau, tau, data$T_obs)
  
  # Grid of unique observed times truncated at tau
  Y.grid <- sort(unique(data$T_obs_tau))
  
  # Indicator for min(T, tau) < C
  data$status_tau <- as.numeric((data$T_obs >= tau) | 
                                  (data$T_obs < tau & data$status == 1))
  
  
  # Estimate propensity scores
  data$e_hat <- estimate_propensity_score(data,
                                          treatment_covariates = X.names.propensity,
                                          type_of_model = nuisance_propensity,
                                          n.folds = n.folds)
  
  
  # Estimation of Q_s
  Q_t <- Q_t_hat(data, tau, X.names.outcome, nuisance, n.folds)
  data$Q_y <-  Q_Y(data,tau,Q_t)
  
  # BJ transformation
  data$Y <-  data$status_tau * data$T_obs_tau + 
    (1 - data$status_tau) * data$Q_y
  
  # IPTW on BJ transformation 
  data$RST <- data$Y * (data$A/data$e_hat-(1-data$A)/(1-data$e_hat))
  
  RMST <- mean(data$RST)
  
  # Return RMST
  return(RMST)
}


# DR censoring transformation
AIPCW <-function(data,
                 tau,
                 X.names.censoring = c("X1","X2","X3","X4"),
                 X.names.outcome = c("X1","X2","X3","X4"),
                 nuisance_Qt = "cox",
                 nuisance_censoring = "cox", 
                 n.folds = NULL, 
                 h_C_hat = NULL,
                 method_aipw = 1) {
  
  # Truncate observed times at tau
  data$T_obs_tau <- pmin(data$T_obs, tau)
  
  # Define status at tau
  data$status_tau <-  as.numeric((data$T_obs > tau) | 
                                   (data$T_obs <= tau &  data$status == 1 ))  
  
  data$censor.status_tau <- 1- as.numeric(
    (data$T_obs > tau) | (data$T_obs <= tau &  data$status == 1 ))
  
  Y.grid <- sort(unique(data$T_obs_tau))
  
  # Estimate survival function for censoring
  S_C_hat <- estimate_survival_function(data = data,X.names.censoring,
                                        type_of_model = nuisance_censoring,
                                        n.folds = n.folds,
                                        Y.grid = Y.grid,
                                        T_obs = "T_obs_tau",
                                        status = "censor.status_tau")
  
  Y.index <- findInterval(data$T_obs_tau, Y.grid)
  
  data$S_C_hat_T_obs_tau <- S_C_hat$S_hat[cbind(seq_along(Y.index), Y.index)]
  
  
  if (is.null(h_C_hat)) {
    h_C_hat <- estimate_hazard_function(S_C_hat$S_hat,Y.grid)
  } 
  
  # Compute Q.t.hat
  Q.t.hat <- Q_t_hat(data = data,
                     X.names = X.names.outcome,
                     tau = tau,
                     nuisance = nuisance_Qt,
                     n.folds = n.folds)
  
  # Compute Q.Y.hat
  data$Q.Y.hat <- Q_Y(data = data, tau, Q.t.hat)
  
  # Compute first term
  data$first_term <- (data$T_obs_tau * data$status_tau) / 
    data$S_C_hat_T_obs_tau
  
  # Compute second term
  data$second_term <- (data$Q.Y.hat * (1 - data$status_tau)) / 
    data$S_C_hat_T_obs_tau
  
  Y.diff <- diff(c(0, Y.grid))
  
  # Compute integrand for the third term
  integrand <- sweep( ( (h_C_hat) / S_C_hat$S_hat )* (Q.t.hat), 2, Y.diff, "*")
  
  # Compute third term
  data$third_term <- integrate(integrand, Y.grid, data$T_obs_tau)
  
  # Compute pseudo outcome
  pseudo_outcome <- data$first_term + data$second_term - data$third_term
  
  return(pseudo_outcome) 
}


AIPTW_AIPCW <- function(data, 
                        tau, 
                        X.names.propensity = c("X1", "X2", "X3", "X4"),
                        X.names.censoring = c("X1", "X2", "X3", "X4"),
                        X.names.outcome = c("X1", "X2", "X3", "X4"),
                        nuisance_propensity = "glm",
                        nuisance_regression = "cox",
                        nuisance_censoring = "cox",
                        nuisance_Qt = "cox",
                        n.folds = NULL) {
  
  # Estimate propensity scores
  data$e_hat <- estimate_propensity_score(
    data = data, 
    treatment_covariates = X.names.propensity, 
    type_of_model = nuisance_propensity, 
    n.folds = n.folds
  )
  
  # Prepare data for censoring model
  data$T_obs_tau <- ifelse(data$T_obs >= tau, tau, data$T_obs)
  
  data$censor.status_tau <- 1 - as.numeric((data$T_obs >= tau) | 
                                             (data$T_obs < tau & data$status == 1))
  
  data$status_tau <- as.numeric((data$T_obs >= tau) | 
                                  (data$T_obs < tau & data$status == 1))
  
  # Create unique time grid
  Y.grid <- sort(unique(data$T_obs_tau))
  
  S_hat <- estimate_survival_function(data, X.names.outcome, 
                                      type_of_model = nuisance_regression, 
                                      Y.grid = Y.grid,
                                      T_obs= "T_obs", 
                                      status = "status", 
                                      n.folds = n.folds)
  
  # Compute area under the survival curve up to tau
  data$E_hat1 <- expected_survival(S_hat$S_hat1, Y.grid)
  data$E_hat0 <- expected_survival(S_hat$S_hat0, Y.grid)
  
  # Compute IPW-weighted residuals
  data$IPW_res <- data$E_hat1 * (1 - data$A / data$e_hat) - 
    data$E_hat0 * (1 - (1 - data$A) / (1 - data$e_hat))
  
  # Compute AIPCW weights
  TDR <- AIPCW(
    data = data, 
    tau = tau,
    X.names.censoring = X.names.censoring,
    X.names.outcome = X.names.outcome,
    nuisance_Qt = nuisance_Qt, 
    nuisance_censoring = nuisance_censoring, 
    n.folds = n.folds
  )
  
  data$TDR <- TDR
  
  # Compute AIPCW-weighted residuals
  data$AIPCW_w <- data$TDR * (data$A / data$e_hat - 
                                (1 - data$A) / (1 - data$e_hat))
  
  # Compute regression residuals
  data$reg <- data$E_hat1 - data$E_hat0
  data$reg_res <- data$A / data$e_hat * (data$TDR - data$E_hat1) - 
    (1 - data$A) / (1 - data$e_hat) * (data$TDR - data$E_hat0)
  
  # Compute estimators
  # na.rm = TRUE to remove NA for the mean calculation
  AIPTW_AIPCW_IPW_res <- mean(data$AIPCW_w + data$IPW_res, na.rm = TRUE)
  AIPTW_AIPCW_reg_res <- mean(data$reg + data$reg_res, na.rm = TRUE)
  
  return(list(AIPTW_AIPCW_reg_res = AIPTW_AIPCW_reg_res, 
              AIPTW_AIPCW_IPW_res = AIPTW_AIPCW_IPW_res))
}


