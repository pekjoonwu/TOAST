## Three Outcome Approach for deciSion-making in clinical Trials (TOAST)
## This R code includes all the basic estimation codes for TOAST framework
## Package: TOAST
## Version: 1.0.0
## Authors: Peijun Wu, Sheng Qiu, Wen-Chi Wu, Qing Li


#' TOAST for continuous type of endpoints
#' @param n1_prev Numeric. Sample size from control group in the previous studies
#' @param n2_prev Numeric. Sample size from treatment group in the previous studies
#' @param mean_1 Numeric. Sample mean from control group in the previous studies
#' @param mean_2 Numeric. Sample mean from treatment group in the previous studies
#' @param std_1 Numeric. Std for the mean parameter from control group in the previous studies
#' @param std_2 Numeric. Std for the mean parameter from treatment group in the previous studies
#' @param Method Character. Method chosen for confidence interval calculation - empirical or bootstrap
#' @param nrep Numeric. Replicates to evaluate the decision once
#' @param NREP Numeric. Replicates to calculate the decision proportions
#' @param sample_size Numeric. Target Sample size in the new studies
#' @param LRV Numeric. Lower reference value (treatment difference) defined prior to the studies
#' @param TV Numeric. Target value (treatment difference) defined prior to the studies
#' @param FGR Numeric. False go risk level
#' @param FSR Numeric. False stop risk level
#' @param seed Numeric. Random seed
#' @return A list containing the population parameters/ TPP values/ number of simulations/ decision results
#' @examples
#' library(TOAST)
#' ## Example code for one continuous primary endpoint
#' cont_examples <- three_outcomes_cont(n1_prev = 50, n2_prev = 50,
#' mean_1 = 1.2, mean_2 = 1.45, std_1 = 0.3, std_2 = 0.3,
#' Method = "empirical", nrep = 500, NREP = 1000, sample_size = 100,
#' FGR = 0.2, FSR = 0.1,
#' LRV = 0.2, TV = 0.25, seed = 2024)
#'
#' @export
three_outcomes_cont <- function(n1_prev, n2_prev, mean_1, mean_2, std_1, std_2,
                                Method = NULL, nrep, NREP, sample_size, LRV, TV,
                                FGR = 0.2, FSR = 0.1, seed) {
  ## Init
  decision_vecs <- c()
  set.seed(seed)
  if (is.null(Method)) {
    Method <- "empirical"
  }

  ## Initializes the progress bar to show the progess of simulations
  pb <- txtProgressBar(min = 0,      ## Minimum value of the progress bar
                       max = NREP, ## Maximum value of the progress bar
                       style = 3,    ## Progress bar style (also available style = 1 and style = 2)
                       width = 50,   ## Progress bar width. Defaults to getOption("width")
                       char = "=")   ## Character used to create the bar

  ## Main loop for each decision
  for (idx in 1:NREP) {
    ## 1. Generate samples for the population parameters
    pop_samples_1 <- rnorm(nrep, mean = mean_1, sd = std_1)
    pop_samples_2 <- rnorm(nrep, mean = mean_2, sd = std_2)

    ## 2. Generate data using the population samples
    ind_samples_1 <- sapply(pop_samples_1, function(x) {
      return(rnorm(sample_size, mean = x, sd = std_1 * sqrt(n1_prev)))
    })
    ind_samples_2 <- sapply(pop_samples_2, function(x) {
      return(rnorm(sample_size, mean = x, sd = std_2 * sqrt(n2_prev)))
    })

    ## 3. Collect the sample mean estimate for each replicates
    mean_est_1 <- colMeans(ind_samples_1)
    mean_est_2 <- colMeans(ind_samples_2)

    ## 3. Calculate the confidence intervals for the treatment difference
    confidence_interval <- function(ind_samples_1, ind_samples_2, interval,
                                    method = c("empirical", "bootstrap")) {
      if (method == "empirical") {
        mean_est_1 <- colMeans(ind_samples_1)
        mean_est_2 <- colMeans(ind_samples_2)
        ## CLT to get the confidence intervals
        vec_sd <- sd(mean_est_2 - mean_est_1)
        n <- length(mean_est_2)
        vec_mean <- mean(mean_est_2 - mean_est_1)
        # Error according to normal distribution - should use one sided critical value
        error <- qt(interval, df = n - 1) * vec_sd / sqrt(n)
        # Confidence interval as a vector
        result <- c("lower" = vec_mean - error, "upper" = vec_mean + error)
        return(result)
      } else if (method == "bootstrap") {
        # Step 2: Define the Statistic Function
        diff_in_mean <- function(data, indices) {
          # Extract control and treatment data
          control_data <- data[indices, 1:nrow(ind_samples_1)]
          treatment_data <- data[indices, (nrow(ind_samples_1) + 1):ncol(data)]

          # Calculate rates for control and treatment groups
          mean_control <- mean(rowMeans(control_data))
          mean_treatment <- mean(rowMeans(treatment_data))

          # Calculate the difference in rates
          return(mean_treatment - mean_control)
        }

        data <- cbind(t(ind_samples_1), t(ind_samples_2))

        # The indices are for selecting which simulations to include in the resampling
        ## Boot can only accept matrix type of input, so need to combine control and treatment into one big matrix
        bootstrap_results <- boot::boot(data = data, statistic = diff_in_mean, R = 1000)

        # Step 4: Construct the Confidence Interval
        boot_ci <- boot::boot.ci(bootstrap_results, type = "perc", conf = interval)

        # Output the bootstrap confidence interval
        return(boot_ci$percent[4:5])
      }
    }

    ## We use different quantiles to calculat the cutoff based on the risk levels
    twenty_percentile <- confidence_interval(ind_samples_1, ind_samples_2,
                                             interval = 1 - FGR, method = Method)
    nighty_percentile <- confidence_interval(ind_samples_1, ind_samples_2,
                                             interval = 1 - FSR, method = Method)

    ## 4. Make a decision based on the decision criteria
    if (nighty_percentile[2] <= TV & twenty_percentile[1] <= LRV) {
    # if (nighty_percentile[2] <= TV) {
      decision_vecs <- c(decision_vecs, "stop")
    } else if (nighty_percentile[2] > TV & twenty_percentile[1] > LRV) {
      decision_vecs <- c(decision_vecs, "go")
    } else {
      decision_vecs <- c(decision_vecs, "consider")
    }

    setTxtProgressBar(pb, idx)
    rm(ind_samples_1)
    rm(ind_samples_2)
  }
  ## Return the results with decision proportions
  res <- list(c(mean_1, std_1, mean_2, std_2), c(LRV, TV), sample_size, decision_vecs)
  names(res) <- c("Population_params", "TPP_values", "sample_size", "decision")
  return(res)
}

#' TOAST for binary type of endpoints
#' @param n1_prev Numeric. Sample size from control group in the previous studies
#' @param n2_prev Numeric. Sample size from treatment group in the previous studies
#' @param p_1 Numeric. Response rate from control group in the previous studies
#' @param p_2 Numeric. Response rate from treatment group in the previous studies
#' @param Method Character. Method chosen for confidence interval calculation - empirical or bootstrap
#' @param nrep Numeric. Replicates to evaluate the decision once
#' @param NREP Numeric. Replicates to calculate the decision proportions
#' @param sample_size Numeric. Target Sample size in the new studies
#' @param LRV Numeric. Lower reference value (treatment difference) defined prior to the studies
#' @param TV Numeric. Target value (treatment difference) defined prior to the studies
#' @param FGR Numeric. False go risk level
#' @param FSR Numeric. False stop risk level
#' @param seed Numeric. Simulation seed for replication purpose
#' @return A list containing the population parameters/ TPP values/ number of simulations/ decision results
#' @examples
#' library(TOAST)
#' ## Example code for one binary primary endpoint
#' binary_examples <- three_outcomes_binary(n1_prev = 60, n2_prev = 60, p_1 = 0.1, p_2 = 0.2,
#'    nrep = 500, NREP = 1000, sample_size = 75,
#'    LRV = 0.093, TV = 0.1, Method = "empirical",
#'    FGR = 0.2, FSR = 0.1, seed = 2024)
#' @export
three_outcomes_binary <- function(n1_prev, n2_prev, p_1, p_2, Method = NULL,
                                  nrep, NREP, sample_size, LRV, TV,
                                  FGR = 0.2, FSR = 0.1, seed) {
  ## Init
  decision_vecs <- c()
  set.seed(seed)
  if (is.null(Method)) {
    Method <- "empirical"
  }

  ## Initializes the progress bar to show the progess of simulations
  pb <- txtProgressBar(min = 0,      ## Minimum value of the progress bar
                       max = NREP, ## Maximum value of the progress bar
                       style = 3,    ## Progress bar style (also available style = 1 and style = 2)
                       width = 50,   ## Progress bar width. Defaults to getOption("width")
                       char = "=")   ## Character used to create the bar

  ## Main loop for each decision
  for (idx in 1:NREP) {
    ## 1. Generate samples for the population parameters
    response_1 <- round(n1_prev * p_1)
    response_2 <- round(n2_prev * p_2)
    pop_sample_1 <- rbeta(nrep, response_1, n1_prev - response_1)
    pop_sample_2 <- rbeta(nrep, response_2, n2_prev - response_2)

    ## 2. Generate data using the population samples
    ind_samples_1 <- sapply(pop_sample_1, function(x) {
      return(rbinom(1, sample_size, prob = x))
    })
    ind_samples_2 <- sapply(pop_sample_2, function(x) {
      return(rbinom(1, sample_size, prob = x))
    })

    ## 3. Calculate the confidence intervals for the treatment difference
    confidence_interval <- function(sample_dat_1, sample_dat_2, interval,
                                    method = c("empirical", "bootstrap")) {
      if (method == "empirical") {
        mean_est_1 <- sample_dat_1/sample_size
        mean_est_2 <- sample_dat_2/sample_size
        ## CLT to get the confidence intervals
        vec_sd <- sd(mean_est_2 - mean_est_1)
        n <- length(mean_est_2)
        vec_mean <- mean(mean_est_2 - mean_est_1)
        # Error according to normal distribution
        error <- qnorm(interval) * vec_sd / sqrt(n)
        # Confidence interval as a vector
        result <- c("lower" = vec_mean - error, "upper" = vec_mean + error)
        return(result)
      } else if (method == "bootstrap") {
        diff_in_proportions <- function(data, indices) {
          control_sample <- data[indices, 1]
          treatment_sample <- data[indices, 2]
          p0_hat <- mean(control_sample) / sample_size
          p1_hat <- mean(treatment_sample) / sample_size
          return(p1_hat - p0_hat)
        }

        # Combine data into a matrix for bootstrapping
        data <- cbind(as.numeric(sample_dat_1), as.numeric(sample_dat_2))

        # Perform bootstrapping
        boot_results <- boot::boot(data = data, statistic = diff_in_proportions, R = 1000)

        # Calculate the bootstrap confidence interval using the percentile method
        boot_ci <- boot::boot.ci(boot_results, type = "perc", conf = interval)

        # Output the bootstrap confidence interval
        return(boot_ci$percent[4:5])
      }
    }
    twenty_percentile <- confidence_interval(ind_samples_1, ind_samples_2,
                                             interval = 1 - FGR, method = Method)
    nighty_percentile <- confidence_interval(ind_samples_1, ind_samples_2,
                                             interval = 1 - FSR, method = Method)

    ## 4. Make a decision based on the decision criteria
    if (nighty_percentile[2] <= TV & twenty_percentile[1] <= LRV) {
    # if (nighty_percentile[2] <= TV) {
      decision_vecs <- c(decision_vecs, "stop")
    } else if (nighty_percentile[2] > TV & twenty_percentile[1] > LRV) {
      decision_vecs <- c(decision_vecs, "go")
    } else {
      decision_vecs <- c(decision_vecs, "consider")
    }

    setTxtProgressBar(pb, idx)
    rm(ind_samples_1)
    rm(ind_samples_2)
  }
  ## Return the results with decision proportions
  res <- list(c(p_1, p_2), c(LRV, TV), sample_size, decision_vecs)
  names(res) <- c("Population_params", "TPP_values", "sample_size", "decision")
  return(res)
}


#' TOAST for count type of endpoints
#' @param shape_1 Numeric. Shape parameter to generate Poisson rate from control group in the previous studies
#' @param shape_2 Numeric. Shape parameter to generate Poisson rate from treatment group in the previous studies
#' @param rate_1 Numeric. Rate parameter to generate Poisson rate from control group in the previous studies
#' @param rate_2 Numeric. Rate parameter to generate Poisson rate from treatment group in the previous studies
#' @param Method Character. Method chosen for confidence interval calculation - empirical or bootstrap
#' @param nrep Numeric. Replicates to evaluate the decision once
#' @param NREP Numeric. Replicates to calculate the decision proportions
#' @param sample_size Numeric. Target Sample size in the new studies
#' @param LRV Numeric. Lower reference value (treatment difference) defined prior to the studies
#' @param TV Numeric. Target value (treatment difference) defined prior to the studies
#' @param FGR Numeric. False go risk level
#' @param FSR Numeric. False stop risk level
#' @param seed Numeric. Simulation seed for replication purpose
#' @return A list containing the population parameters/ TPP values/ number of simulations/ decision results
#' @examples
#' library(TOAST)
#' ## Example code for one count primary endpoint
#' count_examples <- three_outcomes_count(shape_1 = 7, shape_2 = 4, rate_1 = 8.8, rate_2 = 2,
#'    nrep = 500, NREP = 1000, sample_size = 100,
#'    LRV = 1.0, Method = "empirical",
#'    FGR = 0.2, FSR = 0.1, TV = 1.2, seed = 2024)
#' @export
three_outcomes_count <- function(shape_1, shape_2, rate_1, rate_2, Method = NULL,
                                 nrep, NREP, sample_size, LRV, TV,
                                 FGR = 0.2, FSR = 0.1, seed) {
  ## Init
  decision_vecs <- c()
  set.seed(seed)
  if (is.null(Method)) {
    Method <- "empirical"
  }

  ## Initializes the progress bar to show the progess of simulations
  pb <- txtProgressBar(min = 0,      ## Minimum value of the progress bar
                       max = NREP, ## Maximum value of the progress bar
                       style = 3,    ## Progress bar style (also available style = 1 and style = 2)
                       width = 50,   ## Progress bar width. Defaults to getOption("width")
                       char = "=")   ## Character used to create the bar

  ## Main loop for each decision
  for (idx in 1:NREP) {
    ## 1. Generate samples for the population parameters
    pop_sample_1 <- rgamma(n = nrep, shape = shape_1, rate = rate_1)
    pop_sample_2 <- rgamma(n = nrep, shape = shape_2, rate = rate_2)

    ## 2. Generate data using the population samples
    ind_samples_1 <- sapply(pop_sample_1, function(x) {
      return(rpois(n = sample_size, lambda = x))
    })
    ind_samples_2 <- sapply(pop_sample_2, function(x) {
      return(rpois(n = sample_size, lambda = x))
    })

    ## 3. Collect the sample mean estimate for each replicates
    mean_est_1 <- colMeans(ind_samples_1)
    mean_est_2 <- colMeans(ind_samples_2)

    ## 3. Calculate the confidence intervals for the treatment difference
    confidence_interval <- function(ind_samples_1, ind_samples_2, interval,
                                    method = c("empirical", "bootstrap")) {
      if (method == "empirical") {
        mean_est_1 <- colMeans(ind_samples_1)
        mean_est_2 <- colMeans(ind_samples_2)
        ## CLT to get the confidence intervals
        vec_sd <- sd(mean_est_2 - mean_est_1)
        n <- length(mean_est_2)
        vec_mean <- mean(mean_est_2 - mean_est_1)
        # Error according to normal distribution - should use one sided critical value
        error <- qnorm(interval) * vec_sd / sqrt(n)
        # Confidence interval as a vector
        result <- c("lower" = vec_mean - error, "upper" = vec_mean + error)
        return(result)
      } else if (method == "bootstrap") {
        # Step 2: Define the Statistic Function
        diff_in_rates <- function(data, indices) {
          # Extract control and treatment data
          control_data <- data[indices, 1:nrow(ind_samples_1)]
          treatment_data <- data[indices, (nrow(ind_samples_1) + 1):ncol(data)]

          # Calculate rates for control and treatment groups
          rate_control <- mean(rowMeans(control_data))
          rate_treatment <- mean(rowMeans(treatment_data))

          # Calculate the difference in rates
          return(rate_treatment - rate_control)
        }

        data <- cbind(t(ind_samples_1), t(ind_samples_2))

        # The indices are for selecting which simulations to include in the resampling
        ## Boot can only accept matrix type of input, so need to combine control and treatment into one big matrix
        bootstrap_results <- boot::boot(data = data, statistic = diff_in_rates, R = 1000)

        # Step 4: Construct the Confidence Interval
        boot_ci <- boot::boot.ci(bootstrap_results, type = "perc", conf = interval)

        # Output the bootstrap confidence interval
        return(boot_ci$percent[4:5])
      }
    }
    twenty_percentile <- confidence_interval(ind_samples_1, ind_samples_2,
                                             interval = 1 - FGR, method = Method)
    nighty_percentile <- confidence_interval(ind_samples_1, ind_samples_2,
                                             interval = 1 - FSR, method = Method)

    ## 4. Make a decision based on the decision criteria
    if (nighty_percentile[2] <= TV & twenty_percentile[1] <= LRV) {
      # if (nighty_percentile[2] <= TV) {
      decision_vecs <- c(decision_vecs, "stop")
    } else if (nighty_percentile[2] > TV & twenty_percentile[1] > LRV) {
      decision_vecs <- c(decision_vecs, "go")
    } else {
      decision_vecs <- c(decision_vecs, "consider")
    }

    setTxtProgressBar(pb, idx)
    rm(ind_samples_1)
    rm(ind_samples_2)
  }

  ## Return the results with decision proportions
  res <- list(c(shape_1, shape_2, rate_1, rate_2), c(LRV, TV), sample_size, decision_vecs)
  names(res) <- c("Population_params", "TPP_values", "sample_size", "decision")
  return(res)
}

#' TOAST for count type of endpoints using negative binomial distribution
#' @param shape_1 Numeric. Shape parameter to generate NB mean from control group in the previous studies
#' @param shape_2 Numeric. Shape parameter to generate NB mean from treatment group in the previous studies
#' @param rate_1 Numeric. Rate parameter to generate NB mean from control group in the previous studies
#' @param rate_2 Numeric. Rate parameter to generate NB mean from treatment group in the previous studies
#' @param theta_1 Numeric. Variance parameter to generate NB observations from control group in the previous studies
#' @param theta_2 Numeric. Variance parameter to generate NB observations from treatment group in the previous studies
#' @param nrep Numeric. Replicates to evaluate the decision once
#' @param NREP Numeric. Replicates to calculate the decision proportions
#' @param sample_size Numeric. Target Sample size in the new studies
#' @param LRV Numeric. Lower reference value (treatment difference) defined prior to the studies
#' @param TV Numeric. Target value (treatment difference) defined prior to the studies
#' @param FGR Numeric. False go risk level
#' @param FSR Numeric. False stop risk level
#' @param seed Numeric. Simulation seed for replication purpose
#' @return A list containing the population parameters/ TPP values/ number of simulations/ decision results
#' @examples
#' library(TOAST)
#' ## Example code for one negative binomial primary endpoint
#' count_examples <- three_outcomes_count_negbin(shape_1 = 7, shape_2 = 5, rate_1 = 3.5, rate_2 = 2,
#'    theta_1  = 1, theta_2 = 1, Method = "empirical",
#'    nrep = 500, NREP = 1000, sample_size = 100, LRV = 0.35,
#'    FGR = 0.2, FSR = 0.1, TV = 0.5, seed = 2024)
#' @export
three_outcomes_count_negbin <- function(shape_1, shape_2, rate_1, rate_2,
                                        theta_1, theta_2, Method = NULL,
                                        nrep, NREP, sample_size, LRV, TV,
                                        FGR = 0.2, FSR = 0.1, seed) {
  ## Init
  decision_vecs <- c()
  set.seed(seed)
  if (is.null(Method)) {
    Method <- "empirical"
  }

  ## Initializes the progress bar to show the progess of simulations
  pb <- txtProgressBar(min = 0,      ## Minimum value of the progress bar
                       max = NREP, ## Maximum value of the progress bar
                       style = 3,    ## Progress bar style (also available style = 1 and style = 2)
                       width = 50,   ## Progress bar width. Defaults to getOption("width")
                       char = "=")   ## Character used to create the bar

  ## Main loop for each decision
  for (idx in 1:NREP) {
    ## 1. Generate samples for the population parameters
    pop_sample_1 <- rgamma(n = nrep, shape = shape_1, rate = rate_1)
    pop_sample_2 <- rgamma(n = nrep, shape = shape_2, rate = rate_2)

    ## 2. Generate data using the population samples
    ind_samples_1 <- sapply(pop_sample_1, function(x) {
      return(MASS::rnegbin(sample_size, mu = x, theta = theta_1))
    })
    ind_samples_2 <- sapply(pop_sample_2, function(x) {
      return(MASS::rnegbin(sample_size, mu = x, theta = theta_1))
    })

    ## 3. Collect the sample mean estimate for each replicates
    mean_est_1 <- colMeans(ind_samples_1)
    mean_est_2 <- colMeans(ind_samples_2)

    ## 3. Calculate the confidence intervals for the treatment difference
    confidence_interval <- function(ind_samples_1, ind_samples_2, interval,
                                    method = c("empirical", "bootstrap")) {
      if (method == "empirical") {
        mean_est_1 <- colMeans(ind_samples_1)
        mean_est_2 <- colMeans(ind_samples_2)
        ## CLT to get the confidence intervals
        vec_sd <- sd(mean_est_2 - mean_est_1)
        n <- length(mean_est_2)
        vec_mean <- mean(mean_est_2 - mean_est_1)
        # Error according to normal distribution - should use one sided critical value
        error <- qnorm(interval) * vec_sd / sqrt(n)
        # Confidence interval as a vector
        result <- c("lower" = vec_mean - error, "upper" = vec_mean + error)
        return(result)
      } else if (method == "bootstrap") {
        # Step 2: Define the Statistic Function
        diff_in_rates <- function(data, indices) {
          # Extract control and treatment data
          control_data <- data[indices, 1:nrow(ind_samples_1)]
          treatment_data <- data[indices, (nrow(ind_samples_1) + 1):ncol(data)]

          # Calculate rates for control and treatment groups
          rate_control <- mean(rowMeans(control_data))
          rate_treatment <- mean(rowMeans(treatment_data))

          # Calculate the difference in rates
          return(rate_treatment - rate_control)
        }

        data <- cbind(t(ind_samples_1), t(ind_samples_2))

        # The indices are for selecting which simulations to include in the resampling
        ## Boot can only accept matrix type of input, so need to combine control and treatment into one big matrix
        bootstrap_results <- boot::boot(data = data, statistic = diff_in_rates, R = 1000)

        # Step 4: Construct the Confidence Interval
        boot_ci <- boot::boot.ci(bootstrap_results, type = "perc", conf = interval)

        # Output the bootstrap confidence interval
        return(boot_ci$percent[4:5])
      }
    }
    twenty_percentile <- confidence_interval(ind_samples_1, ind_samples_2,
                                             interval = 1 - FGR, method = Method)
    nighty_percentile <- confidence_interval(ind_samples_1, ind_samples_2,
                                             interval = 1 - FSR, method = Method)

    ## 4. Make a decision based on the decision criteria
    if (nighty_percentile[2] <= TV & twenty_percentile[1] <= LRV) {
      # if (nighty_percentile[2] <= TV) {
      decision_vecs <- c(decision_vecs, "stop")
    } else if (nighty_percentile[2] > TV & twenty_percentile[1] > LRV) {
      decision_vecs <- c(decision_vecs, "go")
    } else {
      decision_vecs <- c(decision_vecs, "consider")
    }

    setTxtProgressBar(pb, idx)
    rm(ind_samples_1)
    rm(ind_samples_2)
  }
  ## Return the results with decision proportions
  res <- list(c(shape_1, shape_2, rate_1, rate_2, theta_1, theta_2), c(LRV, TV), sample_size, decision_vecs)
  names(res) <- c("Population_params", "TPP_values", "sample_size", "decision")
  return(res)
}




#' TOAST for time to event type of endpoints (survival data) using Exponential distribution
#' @param hazard_rate_1 Numeric. 1 over Mean estimate for time to event from control group in the previous studies
#' @param HR Numeric. HR comparing the treatment group to the control group
#' @param rate_1 Numeric. Rate parameter to control the failure rate (gamma distribution) from control group in the previous studies
#' @param rate_2 Numeric. Shape parameter to control the failure rate (gamma distribution) from treatment group in the previous studies
#' @param censoring_factor Numeric. Censoring factor - the mean of censoring time would be the original mean times this factor
#' @param nrep Numeric. Replicates to evaluate the decision once
#' @param NREP Numeric. Replicates to calculate the decision proportions
#' @param sample_size Numeric. Target Sample size in the new studies
#' @param LRV Numeric. Lower reference value (treatment difference) defined prior to the studies
#' @param TV Numeric. Target value (treatment difference) defined prior to the studies
#' @param FGR Numeric. False go risk level
#' @param FSR Numeric. False stop risk level
#' @param seed Numeric. Simulation seed for replication purpose
#' @return A list containing the population parameters/ TPP values/ number of simulations/ decision results
#' @examples
#' library(TOAST)
#' ## Example code for one time to event primary endpoint
#' three_outcomes_survival(hazard_rate_1 = 1/60, HR = 0.6,
#'    rate_1 = 1, rate_2 = 1, censoring_factor = 5,
#'    nrep = 500, NREP = 1000, sample_size = 100,
#'    FGR = 0.2, FSR = 0.1, LRV = 0.8, TV = 0.6, seed = 2024)
#' @export
three_outcomes_survival <- function(hazard_rate_1, HR,
                                    rate_1, rate_2, censoring_factor = 2,
                                    nrep, NREP, sample_size, LRV, TV,
                                    FGR = 0.2, FSR = 0.1, seed) {
  ## Init
  decision_vecs <- c()
  set.seed(seed)

  ## Initializes the progress bar to show the progess of simulations
  pb <- txtProgressBar(min = 0,      ## Minimum value of the progress bar
                       max = NREP, ## Maximum value of the progress bar
                       style = 3,    ## Progress bar style (also available style = 1 and style = 2)
                       width = 50,   ## Progress bar width. Defaults to getOption("width")
                       char = "=")   ## Character used to create the bar

  ## Main loop for each decision
  for (idx in 1:NREP) {
    ## 1. Calculate the rate parameters based on the input rate parameters from the control group
    # scale_1 <- mean_TTE_1 / gamma(1 + 1 / shape_1)
    # scale_2 <- scale_1 * (1 / HR)^(1 / shape_1)
    hazard_rate_2 <- hazard_rate_1 * HR

    ## 2. Calculate the gamma shape parameters based on the Weibull scale parameters and input rate parameters
    gamma_shape_1 <- 1/hazard_rate_1 * rate_1
    gamma_shape_2 <- 1/hazard_rate_2 * rate_2

    ## 3. Generate samples for the population parameters
    pop_sample_1 <- rgamma(n = nrep, shape = gamma_shape_1, rate = rate_1)
    pop_sample_2 <- rgamma(n = nrep, shape = gamma_shape_2, rate = rate_2)

    ## 4. Generate data using the population samples (Weibull distribution)
    ind_samples_1 <- sapply(pop_sample_1, function(x) {
      return(rexp(sample_size, rate = 1/x))
    })
    ind_samples_2 <- sapply(pop_sample_2, function(x) {
      return(rexp(sample_size, rate = 1/x))
    })

    ## 5. Simulate the censoring time based on the simulated parameters
    ## and determine the observed times
    censoring_time_1 <- sapply(pop_sample_1, function(x) {
      return(rexp(sample_size, rate = 1/censoring_factor * hazard_rate_1))
    })
    censoring_time_2 <- sapply(pop_sample_2, function(x) {
      return(rexp(sample_size, rate = 1/censoring_factor * hazard_rate_2))
    })

    observed_time_1 <- sapply(1:ncol(ind_samples_1), function(x) {
      return(pmin(ind_samples_1[, x], censoring_time_1[, x]))
    })
    observed_time_2 <- sapply(1:ncol(ind_samples_2), function(x) {
      return(pmin(ind_samples_2[, x], censoring_time_2[, x]))
    })
    event_indicator_1 <- sapply(1:ncol(ind_samples_1), function(x) {
      return(ind_samples_1[, x] <= censoring_time_1[, x])
    })
    event_indicator_2 <- sapply(1:ncol(ind_samples_2), function(x) {
      return(ind_samples_2[, x] <= censoring_time_2[, x])
    })


    ## 6. Fit Cox model to estimate the treatment effect
    cox_results <- sapply(1:ncol(ind_samples_1), function(x) {
      time_to_event_data_tmp <- data.frame(
        observed_time = c(observed_time_1[, x], observed_time_2[, x]),
        event = c(event_indicator_1[, x], event_indicator_2[, x]),
        group = factor(rep(c("A", "B"), times = c(sample_size, sample_size)))
      )
      ## Fit Cox model
      cox_fit <- survival::coxph(Surv(observed_time, event) ~ group, data = time_to_event_data_tmp)
      return(exp(coef(cox_fit)["groupB"]))
    })

    ## 7. Calculate the one sided confidence interval for the harzard ratio using the estimated treatment effect
    confidence_interval <- function(hazard_ratio_vecs, interval) {
      ## CLT to get the confidence intervals
      log_hazard_ratios <- log(hazard_ratio_vecs)
      vec_sd <- sd(log_hazard_ratios)
      n <- length(log_hazard_ratios)
      vec_mean <- mean(log_hazard_ratios)
      # Error according to normal distribution - should use one sided critical value
      error <- qnorm(interval) * vec_sd / sqrt(n)
      # Confidence interval as a vector
      result <- c("lower" = exp(vec_mean - error), "upper" = exp(vec_mean + error))
      return(result)
    }
    twenty_percentile <- confidence_interval(cox_results, interval = 1 - FGR)
    nighty_percentile <- confidence_interval(cox_results, interval = 1 - FSR)

    ## 8. Make a decision based on the decision criteria
    if (nighty_percentile[1] >= TV & twenty_percentile[2] >= LRV) {
      decision_vecs <- c(decision_vecs, "stop")
    } else if (nighty_percentile[1] < TV & twenty_percentile[2] < LRV) {
      decision_vecs <- c(decision_vecs, "go")
    } else {
      decision_vecs <- c(decision_vecs, "consider")
    }

    setTxtProgressBar(pb, idx)
    rm(ind_samples_1)
    rm(ind_samples_2)
  }
  ## Return the results with decision proportions
  res <- list(c(hazard_rate_1, HR, rate_1, rate_2), c(LRV, TV), sample_size, decision_vecs)
  names(res) <- c("Population_params", "TPP_values", "sample_size", "decision")
  return(res)
}
