## Three Outcome Approach for deciSion-making in clinical Trials (TOAST)
## This R code includes functions and simulation codes to perform three outcome approach prediction
## for co-primary endpoints
## Package: TOAST
## Version: 1.0.0
## Authors: Peijun Wu, Sheng Qiu, Wen-Chi Wu, Qing Li

#' Three outcome approach decision making for continuous + continuous type of endpoints
#' @param n1_prev Numeric. Sample size from control group in the previous studies
#' @param n2_prev Numeric. Sample size from treatment group in the previous studies
#' @param end_1_m1 Numeric. Sample mean for endpoints one from control group in the previous studies
#' @param end_1_std1 Numeric. Std for endpoints one from control group in the previous studies
#' @param end_2_m1 Numeric. Sample mean for endpoints two from control group in the previous studies
#' @param end_2_std1 Numeric. Std for endpoints two from control group in the previous studies
#' @param end_1_m2 Numeric. Sample mean for endpoints one from treatment group in the previous studies
#' @param end_1_std2 Numeric. Std for endpoints one from treatment group in the previous studies
#' @param end_2_m2 Numeric. Sample mean for endpoints two from treatment group in the previous studies
#' @param end_2_std2 Numeric. Std for endpoints two from treatment group in the previous studies
#' @param nrep Numeric. Replicates to evaluate the decision once
#' @param NREP Numeric. Replicates to calculate the decision proportions
#' @param sample_size Numeric. Target Sample size in the new studies
#' @param independent Logical. TRUE if two endpoints are independent or FALSE if they are correlated
#' @param correlation Numeric. correlation strength between two endpoints if they are correlated
#' @param strategy Character. Two values: conservative or aggressive for decision criteria
#' @param LRV_1 Numeric. Lower reference value (treatment difference) for endpoint one
#' @param TV_1 Numeric. Target value (treatment difference) for endpoint one
#' @param LRV_2 Numeric. Lower reference value (treatment difference) for endpoint two
#' @param TV_2 Numeric. Target value (treatment difference) for endpoint two
#' @param FGR Numeric. False go risk level
#' @param FSR Numeric. False stop risk level
#' @param seed Numeric. Random seed
#' @examples
#' library(TOAST)
#' ## Example code for two continuous primary endpoint
#' cont_examples_ind <- three_outcomes_cont_coprim(n1_prev = 50, n2_prev = 50,
#'                         end_1_m1 = 1.2, end_1_std1 = 0.3,
#'                         end_2_m1 = 0.6, end_2_std1 = 0.3,
#'                         end_1_m2 = 1.45, end_1_std2 = 0.3,
#'                         end_2_m2 = 0.8, end_2_std2 = 0.3,
#'                         nrep = 300, NREP = 200, sample_size = 100, strategy = "conservative",
#'                         LRV_1 = 0.2, LRV_2 = 0.15, independent = T,
#'                         TV_1 = 0.25, TV_2 = 0.2, FGR = 0.2, FSR = 0.1, seed = 2024)
#' cont_examples_cor <- three_outcomes_cont_coprim(n1_prev = 50, n2_prev = 50,
#'                         end_1_m1 = 1.2, end_1_std1 = 0.3,
#'                         end_2_m1 = 0.6, end_2_std1 = 0.3,
#'                         end_1_m2 = 1.45, end_1_std2 = 0.3,
#'                         end_2_m2 = 0.8, end_2_std2 = 0.3,
#'                         nrep = 300, NREP = 200, sample_size = 100, strategy = "conservative",
#'                         LRV_1 = 0.2, LRV_2 = 0.15, independent = F, correlation = 0.3,
#'                         TV_1 = 0.25, TV_2 = 0.2, FGR = 0.2, FSR = 0.1, seed = 2024)
#'
#' @export
three_outcomes_cont_coprim <- function(n1_prev, n2_prev,
                                           end_1_m1, end_1_std1, end_2_m1, end_2_std1,
                                            end_1_m2, end_1_std2, end_2_m2, end_2_std2,
                                            nrep, NREP, sample_size,
                                           independent = T, correlation = NULL,
                                           strategy = "conservative",
                                           LRV_1, LRV_2,
                                           TV_1, TV_2,
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
    ## Different simulation schemes for correlated or independent endpoints pair
    if (independent == T) {
      ## 1. Generate samples for the population parameters
      control_group_pop_samples <- list(end_points_1 = rnorm(nrep, mean = end_1_m1, sd = end_1_std1),
                                        end_points_2 = rnorm(nrep, mean = end_2_m1, sd = end_2_std1))
      trt_group_pop_samples <- list(end_points_1 = rnorm(nrep, mean = end_1_m2, sd = end_1_std2),
                                    end_points_2 = rnorm(nrep, mean = end_2_m2, sd = end_2_std2))

      ## 2. Generate data using the population samples
      control_group_e1_samples <- sapply(control_group_pop_samples[['end_points_1']], function(x) {
        return(rnorm(sample_size, mean = x, sd = end_1_std1 * sqrt(n1_prev)))
      })
      control_group_e2_samples <- sapply(control_group_pop_samples[['end_points_2']], function(x) {
        return(rnorm(sample_size, mean = x, sd = end_2_std1 * sqrt(n1_prev)))
      })
      trt_group_e1_samples <- sapply(trt_group_pop_samples[['end_points_1']], function(x) {
        return(rnorm(sample_size, mean = x, sd = end_1_std2 * sqrt(n2_prev)))
      })
      trt_group_e2_samples <- sapply(trt_group_pop_samples[['end_points_2']], function(x) {
        return(rnorm(sample_size, mean = x, sd = end_2_std2 * sqrt(n2_prev)))
      })

      ## 3. Calculate the confidence intervals for the treatment difference
      confidence_interval <- function(ind_samples_1, ind_samples_2, interval) {
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
      }
      # For endpoints 1
      twenty_percentile_e1 <- confidence_interval(control_group_e1_samples,
                                                  trt_group_e1_samples, interval = 1 - FGR)
      nighty_percentile_e1 <- confidence_interval(control_group_e1_samples,
                                                  trt_group_e1_samples, interval = 1 - FSR)

      # For endpoints 2
      twenty_percentile_e2 <- confidence_interval(control_group_e2_samples,
                                                  trt_group_e2_samples, interval = 1 - FGR)
      nighty_percentile_e2 <- confidence_interval(control_group_e2_samples,
                                                  trt_group_e2_samples, interval = 1 - FSR)

      ## 4. Make a decision based on the decision criteria
      # For endpoints 1
      if (nighty_percentile_e1[2] <= TV_1 & twenty_percentile_e1[1] <= LRV_1) {
        e1_decision <- "stop"
      } else if (nighty_percentile_e1[2] > TV_1 & twenty_percentile_e1[1] > LRV_1) {
        e1_decision <- "go"
      } else {
        e1_decision <- "consider"
      }
      # For endpoints 2
      if (nighty_percentile_e2[2] <= TV_2 & twenty_percentile_e2[1] <= LRV_2) {
        e2_decision <- "stop"
      } else if (nighty_percentile_e2[2] > TV_2 & twenty_percentile_e2[1] > LRV_2) {
        e2_decision <- "go"
      } else {
        e2_decision <- "consider"
      }
      ## Different decision criteria
      if (strategy == "conservative") {
        decision_vecs <- c(decision_vecs,
                           ifelse(e1_decision == "go" & e2_decision == "go", "go",
                                  ifelse(e1_decision == "stop" | e2_decision == "stop",
                                         "stop", "consider")))
      } else {
        decision_vecs <- c(decision_vecs,
                           ifelse(e1_decision == "stop" | e2_decision == "stop", "stop",
                                  ifelse(e1_decision == "go" | e2_decision == "go", "go",
                                         "consider")))
      }

      setTxtProgressBar(pb, idx)
      rm(control_group_e1_samples, control_group_e2_samples,
         trt_group_e1_samples, trt_group_e2_samples)
    } else {
      ## 1. Generate samples for the population parameters
      control_group_pop_samples <- list(end_points_1 = rnorm(nrep, mean = end_1_m1, sd = end_1_std1),
                                        end_points_2 = rnorm(nrep, mean = end_2_m1, sd = end_2_std1))
      trt_group_pop_samples <- list(end_points_1 = rnorm(nrep, mean = end_1_m2, sd = end_1_std2),
                                    end_points_2 = rnorm(nrep, mean = end_2_m2, sd = end_2_std2))

      ## 2. Generate data using the population samples - correlated
      # Build covariance matrix
      cov_matrix_control <- matrix(c(end_1_std1 ^ 2 * n1_prev,
                                correlation * end_1_std1 * sqrt(n1_prev) * end_2_std1 * sqrt(n1_prev),
                                correlation * end_1_std1 * sqrt(n1_prev) * end_2_std1 * sqrt(n1_prev),
                                end_2_std1 ^ 2 * n1_prev), nrow = 2)
      cov_matrix_trt <- matrix(c(end_1_std2 ^ 2 * n1_prev,
                                correlation * end_1_std2 * sqrt(n1_prev) * end_2_std2 * sqrt(n1_prev),
                                correlation * end_1_std2 * sqrt(n1_prev) * end_2_std2 * sqrt(n1_prev),
                                end_2_std2 ^ 2 * n1_prev), nrow = 2)
      # Generate samples for each group separately
      control_group_samples <- lapply(1:nrep, function(x) {
                MASS::mvrnorm(n = sample_size, mu = c(control_group_pop_samples[['end_points_1']][x],
                                        control_group_pop_samples[['end_points_2']][x]),
                                Sigma = cov_matrix_control)
      })

      trt_group_samples <- lapply(1:nrep, function(x) {
                MASS::mvrnorm(n = sample_size, mu = c(trt_group_pop_samples[['end_points_1']][x],
                                              trt_group_pop_samples[['end_points_2']][x]),
                      Sigma = cov_matrix_trt)
      })

      # Reframe the data structure
      control_group_e1_samples <- do.call(cbind, lapply(control_group_samples, function(x) {
        return(x[, 1, drop = F])
      }))
      control_group_e2_samples <- do.call(cbind, lapply(control_group_samples, function(x) {
        return(x[, 2, drop = F])
      }))
      trt_group_e1_samples <- do.call(cbind, lapply(trt_group_samples, function(x) {
        return(x[, 1, drop = F])
      }))
      trt_group_e2_samples <- do.call(cbind, lapply(trt_group_samples, function(x) {
        return(x[, 2, drop = F])
      }))

      ## 3. Calculate the confidence intervals for the treatment difference
      confidence_interval <- function(ind_samples_1, ind_samples_2, interval) {
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
      }
      # For endpoints 1
      twenty_percentile_e1 <- confidence_interval(control_group_e1_samples,
                                                  trt_group_e1_samples, interval = 1 - FGR)
      nighty_percentile_e1 <- confidence_interval(control_group_e1_samples,
                                                  trt_group_e1_samples, interval = 1 - FSR)

      # For endpoints 2
      twenty_percentile_e2 <- confidence_interval(control_group_e2_samples,
                                                  trt_group_e2_samples, interval = 1 - FGR)
      nighty_percentile_e2 <- confidence_interval(control_group_e2_samples,
                                                  trt_group_e2_samples, interval = 1 - FSR)

      ## 4. Make a decision based on the decision criteria
      # For endpoints 1
      if (nighty_percentile_e1[2] <= TV_1 & twenty_percentile_e1[1] <= LRV_1) {
        e1_decision <- "stop"
      } else if (nighty_percentile_e1[2] > TV_1 & twenty_percentile_e1[1] > LRV_1) {
        e1_decision <- "go"
      } else {
        e1_decision <- "consider"
      }
      # For endpoints 2
      if (nighty_percentile_e2[2] <= TV_2 & twenty_percentile_e2[1] <= LRV_2) {
        e2_decision <- "stop"
      } else if (nighty_percentile_e2[2] > TV_2 & twenty_percentile_e2[1] > LRV_2) {
        e2_decision <- "go"
      } else {
        e2_decision <- "consider"
      }
      ## Different decision criteria
      if (strategy == "conservative") {
        decision_vecs <- c(decision_vecs,
                           ifelse(e1_decision == "go" & e2_decision == "go", "go",
                                  ifelse(e1_decision == "stop" | e2_decision == "stop",
                                         "stop", "consider")))
      } else {
        decision_vecs <- c(decision_vecs,
                           ifelse(e1_decision == "stop" | e2_decision == "stop", "stop",
                                  ifelse(e1_decision == "go" | e2_decision == "go", "go",
                                         "consider")))
      }

      setTxtProgressBar(pb, idx)
      rm(control_group_samples, trt_group_samples, control_group_e1_samples, control_group_e2_samples,
         trt_group_e1_samples, trt_group_e2_samples)
    }

  }
  ## Return the results with decision proportions
  res <- list(sample_size, decision_vecs)
  names(res) <- c("sample_size", "decision")
  return(res)
}

#' Three outcome approach decision making for binary type of endpoints
#' The limits are determined by the CLT only
#' @param n1_prev Numeric. Sample size from control group in the previous studies
#' @param n2_prev Numeric. Sample size from treatment group in the previous studies
#' @param end_1_p1 Numeric. Response rate from control group in the previous studies for endpoints one
#' @param end_2_p1 Numeric. Response rate from control group in the previous studies for endpoints two
#' @param end_1_p2 Numeric. Response rate from treatment group in the previous studies for endpoints one
#' @param end_2_p2 Numeric. Response rate from treatment group in the previous studies for endpoints two
#' @param nrep Numeric. Replicates to evaluate the decision once
#' @param NREP Numeric. Replicates to calculate the decision proportions
#' @param sample_size Numeric. Target Sample size in the new studies
#' @param independent Logical. TRUE if two endpoints are independent or FALSE if they are correlated
#' @param correlation Numeric. correlation strength between two endpoints if they are correlated
#' @param strategy Character. Two values: conservative or aggressive for decision criteria
#' @param LRV_1 Numeric. Lower reference value (treatment difference) for endpoint one
#' @param TV_1 Numeric. Target value (treatment difference) for endpoint one
#' @param LRV_2 Numeric. Lower reference value (treatment difference) for endpoint two
#' @param TV_2 Numeric. Target value (treatment difference) for endpoint two
#' @param FGR Numeric. False go risk level
#' @param FSR Numeric. False stop risk level
#' @param seed Numeric. Random seed
#' @examples
#' library(TOAST)
#' ## Example code for two binary primary endpoint
#' binary_examples_ind <- three_outcomes_binary_coprim(n1_prev = 50, n2_prev = 50,
#'                         end_1_p1 = 0.1, end_2_p1 = 0.23,
#'                         end_1_p2 = 0.2, end_2_p2 = 0.35,
#'                         nrep = 300, NREP = 200, sample_size = 55, strategy = "conservative",
#'                         LRV_1 = 0.05, LRV_2 = 0.08, independent = T,
#'                         TV_1 = 0.1, TV_2 = 0.12, FGR = 0.2, FSR = 0.1, seed = 2024)
#' binary_examples_cor <- three_outcomes_binary_coprim(n1_prev = 50, n2_prev = 50,
#'                         end_1_p1 = 0.1, end_2_p1 = 0.23,
#'                         end_1_p2 = 0.2, end_2_p2 = 0.35,
#'                         nrep = 300, NREP = 200, sample_size = 55, strategy = "conservative",
#'                         LRV_1 = 0.05, LRV_2 = 0.08, independent = F, correlation = 0.3,
#'                         TV_1 = 0.1, TV_2 = 0.12, FGR = 0.2, FSR = 0.1, seed = 2024)
#'
#' @export
three_outcomes_binary_coprim <- function(n1_prev, n2_prev,
                                         end_1_p1, end_2_p1,
                                         end_1_p2, end_2_p2,
                                         nrep, NREP, sample_size,
                                         independent = T, correlation = NULL,
                                         strategy = "conservative",
                                         LRV_1, LRV_2,
                                         TV_1, TV_2,
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
    ## Different simulation schemes for correlated or independent endpoints pair
    if (independent == T) {
      ## 1. Generate samples for the population parameters
      control_response_e1 <- round(n1_prev * end_1_p1)
      control_response_e2 <- round(n1_prev * end_2_p1)
      trt_response_e1 <- round(n2_prev * end_1_p2)
      trt_response_e2 <- round(n2_prev * end_2_p2)
      control_group_pop_samples <- list(end_points_1 = rbeta(nrep, control_response_e1,
                                                             n1_prev - control_response_e1),
                                        end_points_2 = rbeta(nrep, control_response_e2,
                                                             n1_prev - control_response_e2))
      trt_group_pop_samples <- list(end_points_1 = rbeta(nrep, trt_response_e1,
                                                         n2_prev - trt_response_e1),
                                    end_points_2 = rbeta(nrep, trt_response_e2,
                                                         n2_prev - trt_response_e2))

      ## 2. Generate data using the population samples
      control_group_e1_samples <- sapply(control_group_pop_samples[['end_points_1']], function(x) {
        return(rbinom(1, sample_size, prob = x))
      })
      control_group_e2_samples <- sapply(control_group_pop_samples[['end_points_2']], function(x) {
        return(rbinom(1, sample_size, prob = x))
      })
      trt_group_e1_samples <- sapply(trt_group_pop_samples[['end_points_1']], function(x) {
        return(rbinom(1, sample_size, prob = x))
      })
      trt_group_e2_samples <- sapply(trt_group_pop_samples[['end_points_2']], function(x) {
        return(rbinom(1, sample_size, prob = x))
      })

      ## 3. Calculate the confidence intervals for the treatment difference
      confidence_interval <- function(sample_dat_1, sample_dat_2, interval) {
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
      }

      # For endpoints 1
      twenty_percentile_e1 <- confidence_interval(control_group_e1_samples,
                                                  trt_group_e1_samples, interval = 1 - FGR)
      nighty_percentile_e1 <- confidence_interval(control_group_e1_samples,
                                                  trt_group_e1_samples, 0.9)

      # For endpoints 2
      twenty_percentile_e2 <- confidence_interval(control_group_e2_samples,
                                                  trt_group_e2_samples, interval = 1 - FGR)
      nighty_percentile_e2 <- confidence_interval(control_group_e2_samples,
                                                  trt_group_e2_samples, 0.9)

      ## 4. Make a decision based on the decision criteria
      # For endpoints 1
      if (nighty_percentile_e1[2] <= TV_1 & twenty_percentile_e1[1] <= LRV_1) {
        e1_decision <- "stop"
      } else if (nighty_percentile_e1[2] > TV_1 & twenty_percentile_e1[1] > LRV_1) {
        e1_decision <- "go"
      } else {
        e1_decision <- "consider"
      }
      # For endpoints 2
      if (nighty_percentile_e2[2] <= TV_2 & twenty_percentile_e2[1] <= LRV_2) {
        e2_decision <- "stop"
      } else if (nighty_percentile_e2[2] > TV_2 & twenty_percentile_e2[1] > LRV_2) {
        e2_decision <- "go"
      } else {
        e2_decision <- "consider"
      }
      ## Different decision criteria
      if (strategy == "conservative") {
        decision_vecs <- c(decision_vecs,
                           ifelse(e1_decision == "go" & e2_decision == "go", "go",
                                  ifelse(e1_decision == "stop" | e2_decision == "stop",
                                         "stop", "consider")))
      } else {
        decision_vecs <- c(decision_vecs,
                           ifelse(e1_decision == "stop" | e2_decision == "stop", "stop",
                                  ifelse(e1_decision == "go" | e2_decision == "go", "go",
                                         "consider")))
      }

      setTxtProgressBar(pb, idx)
      rm(control_group_e1_samples, control_group_e2_samples,
         trt_group_e1_samples, trt_group_e2_samples)
    } else {
      ## 1. Generate samples for the population parameters
      control_response_e1 <- round(n1_prev * end_1_p1)
      control_response_e2 <- round(n1_prev * end_2_p1)
      trt_response_e1 <- round(n2_prev * end_1_p2)
      trt_response_e2 <- round(n2_prev * end_2_p2)
      control_group_pop_samples <- list(end_points_1 = rbeta(nrep, control_response_e1,
                                                             n1_prev - control_response_e1),
                                        end_points_2 = rbeta(nrep, control_response_e2,
                                                             n1_prev - control_response_e2))
      trt_group_pop_samples <- list(end_points_1 = rbeta(nrep, trt_response_e1,
                                                         n2_prev - trt_response_e1),
                                    end_points_2 = rbeta(nrep, trt_response_e2,
                                                         n2_prev - trt_response_e2))
      ## 2. Generate data using the population samples - correlated
      # Build correaltion matrix
      cor_matrix <- matrix(c(1, correlation,
                                     correlation, 1), 2, 2)
      # Generate samples for each group separately
      # control_group_samples <- lapply(1:nrep, function(x) {
      #   tmp_sim <- simstudy::genCorGen(sample_size, nvars = 2, params1 = c(control_group_pop_samples[['end_points_1']][x],
      #                                                 control_group_pop_samples[['end_points_2']][x]),
      #             dist = "binary", corMatrix = cor_matrix, wide = TRUE)
      #   return(colMeans(tmp_sim[, c("V1", "V2")]))
      # })
      #
      # trt_group_samples <- lapply(1:nrep, function(x) {
      #   tmp_sim <- simstudy::genCorGen(sample_size, nvars = 2, params1 = c(trt_group_pop_samples[['end_points_1']][x],
      #                                                            trt_group_pop_samples[['end_points_2']][x]),
      #                        dist = "binary", corMatrix = cor_matrix, wide = TRUE)
      #   return(colMeans(tmp_sim[, c("V1", "V2")]))
      # })
      #
      copula_simulation <- function(cor_matrix, distribution = NULL, sample_size,
                                    p1, p2,
                                    rate1, rate2) {
        ## 1. Generate a Gaussian copula with the defined correlation
        # gaussian_cop <- copula::normalCopula(param = cor_params, dim = 2, dispstr = "un")
        gaussian_cop <- MASS::mvrnorm(sample_size, mu = c(0, 0), Sigma = cor_matrix)

        ## 2. Simulate correlated uniform variables
        # u <- copula::rCopula(sample_size, gaussian_cop)
        u <- pnorm(gaussian_cop)

        ## 3. Transform uniform variables to binary
        if (distribution == "binary") {
          binary_var1 <- as.numeric(u[, 1] < p1)
          binary_var2 <- as.numeric(u[, 2] < p2)

          ## Combine into a data frame
          simulated_data <- data.frame(binary_var1, binary_var2)
          return(simulated_data)
        } else if (distribution == "poisson") {
          poisson_var1 <- qpois(u[, 1], lambda = rate1)
          poisson_var2 <- qpois(u[, 2], lambda = rate2)

          ## Combine into a data frame
          simulated_data <- data.frame(poisson_var1, poisson_var2)
          return(simulated_data)
        }
      }
      control_group_samples <- lapply(1:nrep, function(x) {
        tmp_sim <- copula_simulation(cor_matrix = cor_matrix, distribution = "binary",
                                     sample_size = sample_size,
                                     p1 = control_group_pop_samples[['end_points_1']][x],
                                     p2 = control_group_pop_samples[['end_points_2']][x])
        return(colMeans(tmp_sim))
      })

      trt_group_samples <- lapply(1:nrep, function(x) {
        tmp_sim <- copula_simulation(cor_matrix = cor_matrix, distribution = "binary",
                                     sample_size = sample_size,
                                     p1 = trt_group_pop_samples[['end_points_1']][x],
                                     p2 = trt_group_pop_samples[['end_points_2']][x])
        return(colMeans(tmp_sim))
      })

      # Reframe the data structure
      control_group_e1_samples <- do.call(c, lapply(control_group_samples, function(x) {
        return(x[1])
      }))
      control_group_e2_samples <- do.call(c, lapply(control_group_samples, function(x) {
        return(x[2])
      }))
      trt_group_e1_samples <- do.call(c, lapply(trt_group_samples, function(x) {
        return(x[1])
      }))
      trt_group_e2_samples <- do.call(c, lapply(trt_group_samples, function(x) {
        return(x[2])
      }))

      ## 3. Calculate the confidence intervals for the treatment difference
      confidence_interval <- function(sample_dat_1, sample_dat_2, interval) {
        mean_est_1 <- sample_dat_1
        mean_est_2 <- sample_dat_2
        ## CLT to get the confidence intervals
        vec_sd <- sd(mean_est_2 - mean_est_1)
        n <- length(mean_est_2)
        vec_mean <- mean(mean_est_2 - mean_est_1)
        # Error according to normal distribution
        error <- qnorm(interval) * vec_sd / sqrt(n)
        # Confidence interval as a vector
        result <- c("lower" = vec_mean - error, "upper" = vec_mean + error)
        return(result)
      }
      # For endpoints 1
      twenty_percentile_e1 <- confidence_interval(control_group_e1_samples,
                                                  trt_group_e1_samples, interval = 1 - FGR)
      nighty_percentile_e1 <- confidence_interval(control_group_e1_samples,
                                                  trt_group_e1_samples, 0.9)

      # For endpoints 2
      twenty_percentile_e2 <- confidence_interval(control_group_e2_samples,
                                                  trt_group_e2_samples, interval = 1 - FGR)
      nighty_percentile_e2 <- confidence_interval(control_group_e2_samples,
                                                  trt_group_e2_samples, 0.9)

      ## 4. Make a decision based on the decision criteria
      # For endpoints 1
      if (nighty_percentile_e1[2] <= TV_1 & twenty_percentile_e1[1] <= LRV_1) {
        e1_decision <- "stop"
      } else if (nighty_percentile_e1[2] > TV_1 & twenty_percentile_e1[1] > LRV_1) {
        e1_decision <- "go"
      } else {
        e1_decision <- "consider"
      }
      # For endpoints 2
      if (nighty_percentile_e2[2] <= TV_2 & twenty_percentile_e2[1] <= LRV_2) {
        e2_decision <- "stop"
      } else if (nighty_percentile_e2[2] > TV_2 & twenty_percentile_e2[1] > LRV_2) {
        e2_decision <- "go"
      } else {
        e2_decision <- "consider"
      }
      ## Different decision criteria
      if (strategy == "conservative") {
        decision_vecs <- c(decision_vecs,
                           ifelse(e1_decision == "go" & e2_decision == "go", "go",
                                  ifelse(e1_decision == "stop" | e2_decision == "stop",
                                         "stop", "consider")))
      } else {
        decision_vecs <- c(decision_vecs,
                           ifelse(e1_decision == "stop" | e2_decision == "stop", "stop",
                                  ifelse(e1_decision == "go" | e2_decision == "go", "go",
                                         "consider")))
      }

      setTxtProgressBar(pb, idx)
      rm(control_group_samples, trt_group_samples, control_group_e1_samples, control_group_e2_samples,
         trt_group_e1_samples, trt_group_e2_samples)
    }
  }
  ## Return the results with decision proportions
  res <- list(sample_size, decision_vecs)
  names(res) <- c("sample_size", "decision")
  return(res)
}




#' Three outcome approach decision making for count type of endpoints
#' The limits are determined by the CLT only
#' @param end_1_shape1 Numeric. Shape parameter to generate Poisson rate for endpoints one from control group in the previous studies
#' @param end_1_rate1 Numeric. Rate parameter to generate Poisson rate for endpoints one from control group in the previous studies
#' @param end_2_shape1 Numeric. Shape parameter to generate Poisson rate for endpoints two from control group in the previous studies
#' @param end_2_rate1 Numeric. Rate parameter to generate Poisson rate for endpoints two from control group in the previous studies
#' @param end_1_shape2 Numeric. Shape parameter to generate Poisson rate for endpoints one from treatment group in the previous studies
#' @param end_1_rate2 Numeric. Rate parameter to generate Poisson rate for endpoints one from treatment group in the previous studies
#' @param end_2_shape2 Numeric. Shape parameter to generate Poisson rate for endpoints two from treatment group in the previous studies
#' @param end_2_rate2 Numeric. Rate parameter to generate Poisson rate for endpoints two from treatment group in the previous studies
#' @param nrep Numeric. Replicates to evaluate the decision once
#' @param NREP Numeric. Replicates to calculate the decision proportions
#' @param sample_size Numeric. Target Sample size in the new studies
#' @param independent Logical. TRUE if two endpoints are independent or FALSE if they are correlated
#' @param correlation Numeric. correlation strength between two endpoints if they are correlated
#' @param strategy Character. Two values: conservative or aggressive for decision criteria
#' @param LRV_1 Numeric. Lower reference value (treatment difference) for endpoint one
#' @param TV_1 Numeric. Target value (treatment difference) for endpoint one
#' @param LRV_2 Numeric. Lower reference value (treatment difference) for endpoint two
#' @param TV_2 Numeric. Target value (treatment difference) for endpoint two
#' @param FGR Numeric. False go risk level
#' @param FSR Numeric. False stop risk level
#' @param seed Numeric. Random seed
#' @examples
#' library(TOAST)
#' ## Example code for two binary primary endpoint
#' count_examples_ind <- three_outcomes_count_coprim(end_1_shape1 = 7, end_1_rate1 = 8.8,
#'                                                   end_2_shape1 = 3, end_2_rate1 = 5,
#'                                                   end_1_shape2 = 4, end_1_rate2 = 2,
#'                                                  end_2_shape2 = 5, end_2_rate2 = 6,
#'                         nrep = 300, NREP = 200, sample_size = 100, strategy = "conservative",
#'                         LRV_1 = 1.0, LRV_2 = 0.1, independent = T,
#'                         TV_1 = 1.204, TV_2 = 0.233, FGR = 0.2, FSR = 0.1, seed = 2024)
#' count_examples_cor <- three_outcomes_count_coprim(end_1_shape1 = 7, end_1_rate1 = 8.8,
#'                                                   end_2_shape1 = 3, end_2_rate1 = 5,
#'                                                   end_1_shape2 = 4, end_1_rate2 = 2,
#'                                                  end_2_shape2 = 5, end_2_rate2 = 6,
#'                         nrep = 300, NREP = 200, sample_size = 100, strategy = "conservative",
#'                         LRV_1 = 1.0, LRV_2 = 0.1, independent = F, correlation = 0.3,
#'                         TV_1 = 1.204, TV_2 = 0.233, FGR = 0.2, FSR = 0.1, seed = 2024)
#'
#' @export
three_outcomes_count_coprim <- function(end_1_shape1, end_1_rate1, end_2_shape1, end_2_rate1,
                                        end_1_shape2, end_1_rate2, end_2_shape2, end_2_rate2,
                                         nrep, NREP, sample_size,
                                         independent = T, correlation = NULL,
                                        strategy = "conservative",
                                        LRV_1, LRV_2,
                                        TV_1, TV_2,
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
    ## Different simulation schemes for correlated or independent endpoints pair
    if (independent == T) {
      ## 1. Generate samples for the population parameters
      control_group_pop_samples <- list(end_points_1 = rgamma(n = nrep,
                                                              shape = end_1_shape1,
                                                              rate = end_1_rate1),
                                        end_points_2 = rgamma(n = nrep,
                                                              shape = end_2_shape1,
                                                              rate = end_2_rate1))
      trt_group_pop_samples <- list(end_points_1 = rgamma(n = nrep,
                                                          shape = end_1_shape2,
                                                          rate = end_1_rate2),
                                    end_points_2 = rgamma(n = nrep,
                                                          shape = end_2_shape2,
                                                          rate = end_2_rate2))

      ## 2. Generate data using the population samples
      control_group_e1_samples <- sapply(control_group_pop_samples[['end_points_1']], function(x) {
        return(rpois(n = sample_size, lambda = x))
      })
      control_group_e2_samples <- sapply(control_group_pop_samples[['end_points_2']], function(x) {
        return(rpois(n = sample_size, lambda = x))
      })
      trt_group_e1_samples <- sapply(trt_group_pop_samples[['end_points_1']], function(x) {
        return(rpois(n = sample_size, lambda = x))
      })
      trt_group_e2_samples <- sapply(trt_group_pop_samples[['end_points_2']], function(x) {
        return(rpois(n = sample_size, lambda = x))
      })

      ## 3. Calculate the confidence intervals for the treatment difference
      confidence_interval <- function(sample_dat_1, sample_dat_2, interval) {
        mean_est_1 <- colMeans(sample_dat_1)
        mean_est_2 <- colMeans(sample_dat_2)
        ## CLT to get the confidence intervals
        vec_sd <- sd(mean_est_2 - mean_est_1)
        n <- length(mean_est_2)
        vec_mean <- mean(mean_est_2 - mean_est_1)
        # Error according to normal distribution - should use one sided critical value
        error <- qnorm(interval) * vec_sd / sqrt(n)
        # Confidence interval as a vector
        result <- c("lower" = vec_mean - error, "upper" = vec_mean + error)
        return(result)
      }

      # For endpoints 1
      twenty_percentile_e1 <- confidence_interval(control_group_e1_samples,
                                                  trt_group_e1_samples, interval = 1 - FGR)
      nighty_percentile_e1 <- confidence_interval(control_group_e1_samples,
                                                  trt_group_e1_samples, 0.9)

      # For endpoints 2
      twenty_percentile_e2 <- confidence_interval(control_group_e2_samples,
                                                  trt_group_e2_samples, interval = 1 - FGR)
      nighty_percentile_e2 <- confidence_interval(control_group_e2_samples,
                                                  trt_group_e2_samples, 0.9)

      ## 4. Make a decision based on the decision criteria
      # For endpoints 1
      if (nighty_percentile_e1[2] <= TV_1 & twenty_percentile_e1[1] <= LRV_1) {
        e1_decision <- "stop"
      } else if (nighty_percentile_e1[2] > TV_1 & twenty_percentile_e1[1] > LRV_1) {
        e1_decision <- "go"
      } else {
        e1_decision <- "consider"
      }
      # For endpoints 2
      if (nighty_percentile_e2[2] <= TV_2 & twenty_percentile_e2[1] <= LRV_2) {
        e2_decision <- "stop"
      } else if (nighty_percentile_e2[2] > TV_2 & twenty_percentile_e2[1] > LRV_2) {
        e2_decision <- "go"
      } else {
        e2_decision <- "consider"
      }
      ## Different decision criteria
      if (strategy == "conservative") {
        decision_vecs <- c(decision_vecs,
                           ifelse(e1_decision == "go" & e2_decision == "go", "go",
                                  ifelse(e1_decision == "stop" | e2_decision == "stop",
                                         "stop", "consider")))
      } else {
        decision_vecs <- c(decision_vecs,
                           ifelse(e1_decision == "stop" | e2_decision == "stop", "stop",
                                  ifelse(e1_decision == "go" | e2_decision == "go", "go",
                                         "consider")))
      }

      setTxtProgressBar(pb, idx)
      rm(control_group_e1_samples, control_group_e2_samples,
         trt_group_e1_samples, trt_group_e2_samples)
    } else {
      ## 1. Generate samples for the population parameters
      control_group_pop_samples <- list(end_points_1 = rgamma(n = nrep,
                                                              shape = end_1_shape1,
                                                              rate = end_1_rate1),
                                        end_points_2 = rgamma(n = nrep,
                                                              shape = end_2_shape1,
                                                              rate = end_2_rate1))
      trt_group_pop_samples <- list(end_points_1 = rgamma(n = nrep,
                                                          shape = end_1_shape2,
                                                          rate = end_1_rate2),
                                    end_points_2 = rgamma(n = nrep,
                                                          shape = end_2_shape2,
                                                          rate = end_2_rate2))

      ## 2. Generate data using the population samples - correlated
      # Build correaltion matrix
      cor_matrix <- matrix(c(1, correlation,
                             correlation, 1), 2, 2)

      copula_simulation <- function(cor_matrix, distribution = NULL, sample_size,
                                    p1 = NULL, p2 = NULL,
                                    rate1 = NULL, rate2 = NULL) {
        ## 1. Generate a Gaussian copula with the defined correlation
        # gaussian_cop <- copula::normalCopula(param = cor_params, dim = 2, dispstr = "un")
        gaussian_cop <- MASS::mvrnorm(sample_size, mu = c(0, 0), Sigma = cor_matrix)

        ## 2. Simulate correlated uniform variables
        # u <- copula::rCopula(sample_size, gaussian_cop)
        u <- pnorm(gaussian_cop)

        ## 3. Transform uniform variables to binary
        if (distribution == "binary") {
          if (is.null(p1) || is.null(p2)) stop("p1 and p2 must be provided for binary distribution")
          binary_var1 <- as.numeric(u[, 1] < p1)
          binary_var2 <- as.numeric(u[, 2] < p2)

          ## Combine into a data frame
          simulated_data <- data.frame(binary_var1, binary_var2)
          return(simulated_data)
        } else if (distribution == "poisson") {
          if (is.null(rate1) || is.null(rate2)) stop("rate1 and rate2 must be provided for Poisson distribution")
          poisson_var1 <- qpois(u[, 1], lambda = rate1)
          poisson_var2 <- qpois(u[, 2], lambda = rate2)

          ## Combine into a data frame
          simulated_data <- data.frame(poisson_var1, poisson_var2)
          return(simulated_data)
        } else {
          stop("Unsupported distribution type")
        }
      }

      ## Correlatin here needs some adjustment
      # Generate samples for each group separately
      control_group_samples <- lapply(1:nrep, function(x) {
        tmp_sim <- copula_simulation(cor_matrix = cor_matrix, distribution = "poisson",
                                     sample_size = sample_size,
                                     rate1 = control_group_pop_samples[['end_points_1']][x],
                                     rate2 = control_group_pop_samples[['end_points_2']][x])
        return(colMeans(tmp_sim))
      })

      trt_group_samples <- lapply(1:nrep, function(x) {
        tmp_sim <- copula_simulation(cor_matrix = cor_matrix, distribution = "poisson",
                                     sample_size = sample_size,
                                     rate1 = trt_group_pop_samples[['end_points_1']][x],
                                     rate2 = trt_group_pop_samples[['end_points_2']][x])
        return(colMeans(tmp_sim))
      })
      # # Generate samples for each group separately
      # control_group_samples <- lapply(1:nrep, function(x) {
      #   tmp_sim <- simstudy::genCorGen(sample_size, nvars = 2, params1 = c(control_group_pop_samples[['end_points_1']][x],
      #                                                                      control_group_pop_samples[['end_points_2']][x]),
      #                                  dist = "poisson", corMatrix = cor_matrix, wide = TRUE)
      #   return(tmp_sim[, c("V1", "V2")])
      # })
      #
      # trt_group_samples <- lapply(1:nrep, function(x) {
      #   tmp_sim <- simstudy::genCorGen(sample_size, nvars = 2, params1 = c(trt_group_pop_samples[['end_points_1']][x],
      #                                                                      trt_group_pop_samples[['end_points_2']][x]),
      #                                  dist = "poisson", corMatrix = cor_matrix, wide = TRUE)
      #   return(tmp_sim[, c("V1", "V2")])
      # })

      # Reframe the data structure
      control_group_e1_samples <- do.call(cbind, lapply(control_group_samples, function(x) {
        return(x[1])
      }))
      control_group_e2_samples <- do.call(cbind, lapply(control_group_samples, function(x) {
        return(x[2])
      }))
      trt_group_e1_samples <- do.call(cbind, lapply(trt_group_samples, function(x) {
        return(x[1])
      }))
      trt_group_e2_samples <- do.call(cbind, lapply(trt_group_samples, function(x) {
        return(x[2])
      }))

      ## 3. Calculate the confidence intervals for the treatment difference
      confidence_interval <- function(sample_dat_1, sample_dat_2, interval) {
        mean_est_1 <- colMeans(sample_dat_1)
        mean_est_2 <- colMeans(sample_dat_2)
        ## CLT to get the confidence intervals
        vec_sd <- sd(mean_est_2 - mean_est_1)
        n <- length(mean_est_2)
        vec_mean <- mean(mean_est_2 - mean_est_1)
        # Error according to normal distribution
        error <- qnorm(interval) * vec_sd / sqrt(n)
        # Confidence interval as a vector
        result <- c("lower" = vec_mean - error, "upper" = vec_mean + error)
        return(result)
      }

      # For endpoints 1
      twenty_percentile_e1 <- confidence_interval(control_group_e1_samples,
                                                  trt_group_e1_samples, interval = 1 - FGR)
      nighty_percentile_e1 <- confidence_interval(control_group_e1_samples,
                                                  trt_group_e1_samples, 0.9)

      # For endpoints 2
      twenty_percentile_e2 <- confidence_interval(control_group_e2_samples,
                                                  trt_group_e2_samples, interval = 1 - FGR)
      nighty_percentile_e2 <- confidence_interval(control_group_e2_samples,
                                                  trt_group_e2_samples, 0.9)

      ## 4. Make a decision based on the decision criteria
      # For endpoints 1
      if (nighty_percentile_e1[2] <= TV_1 & twenty_percentile_e1[1] <= LRV_1) {
        e1_decision <- "stop"
      } else if (nighty_percentile_e1[2] > TV_1 & twenty_percentile_e1[1] > LRV_1) {
        e1_decision <- "go"
      } else {
        e1_decision <- "consider"
      }
      # For endpoints 2
      if (nighty_percentile_e2[2] <= TV_2 & twenty_percentile_e2[1] <= LRV_2) {
        e2_decision <- "stop"
      } else if (nighty_percentile_e2[2] > TV_2 & twenty_percentile_e2[1] > LRV_2) {
        e2_decision <- "go"
      } else {
        e2_decision <- "consider"
      }
      ## Different decision criteria
      if (strategy == "conservative") {
        decision_vecs <- c(decision_vecs,
                           ifelse(e1_decision == "go" & e2_decision == "go", "go",
                                  ifelse(e1_decision == "stop" | e2_decision == "stop",
                                         "stop", "consider")))
      } else {
        decision_vecs <- c(decision_vecs,
                           ifelse(e1_decision == "stop" | e2_decision == "stop", "stop",
                                  ifelse(e1_decision == "go" | e2_decision == "go", "go",
                                         "consider")))
      }

      setTxtProgressBar(pb, idx)
      rm(control_group_samples, trt_group_samples, control_group_e1_samples, control_group_e2_samples,
         trt_group_e1_samples, trt_group_e2_samples)
    }
  }
  ## Return the results with decision proportions
  res <- list(sample_size, decision_vecs)
  names(res) <- c("sample_size", "decision")
  return(res)
}



#' Three outcome approach decision making for continous + binary co-primary endpoints
#' The limits are determined by the CLT only
#' @param n1_prev Numeric. Sample size from control group in the previous studies
#' @param n2_prev Numeric. Sample size from treatment group in the previous studies
#' @param end_1_m1 Numeric. Sample mean for endpoints one from control group in the previous studies
#' @param end_1_std1 Numeric. Std for endpoints one from control group in the previous studies
#' @param end_1_m2 Numeric. Sample mean for endpoints one from treatment group in the previous studies
#' @param end_1_std2 Numeric. Std for endpoints one from treatment group in the previous studies
#' @param end_2_p1 Numeric. Response rate from control group in the previous studies for endpoints two
#' @param end_2_p2 Numeric. Response rate from treatment group in the previous studies for endpoints two
#' @param nrep Numeric. Replicates to evaluate the decision once
#' @param NREP Numeric. Replicates to calculate the decision proportions
#' @param sample_size Numeric. Target Sample size in the new studies
#' @param independent Logical. TRUE if two endpoints are independent or FALSE if they are correlated
#' @param correlation Numeric. correlation strength between two endpoints if they are correlated
#' @param strategy Character. Two values: conservative or aggressive for decision criteria
#' @param LRV_1 Numeric. Lower reference value (treatment difference) for endpoint one
#' @param TV_1 Numeric. Target value (treatment difference) for endpoint one
#' @param LRV_2 Numeric. Lower reference value (treatment difference) for endpoint two
#' @param TV_2 Numeric. Target value (treatment difference) for endpoint two
#' @param FGR Numeric. False go risk level
#' @param FSR Numeric. False stop risk level
#' @param seed Numeric. Random seed
#' @examples
#' library(TOAST)
#' ## Example code for two binary primary endpoint
#' cont_bin_examples_ind <- three_outcomes_cont_binary_coprimary(n1_prev = 50, n2_prev = 50,
#'                                               end_1_m1 = 1.2, end_1_std1 = 0.3,
#'                                               end_1_m2 = 1.45, end_1_std2 = 0.3,
#'                                               end_2_p1 = 0.3, end_2_p2 = 0.5,
#'                         nrep = 300, NREP = 200, sample_size = 100, strategy = "conservative",
#'                         LRV_1 = 0.2, LRV_2 = 0.1, independent = T,
#'                         TV_1 = 0.25, TV_2 = 0.2, FGR = 0.2, FSR = 0.1, seed = 2024)
#' cont_bin_examples_cor <- three_outcomes_cont_binary_coprimary(n1_prev = 50, n2_prev = 50,
#'                                               end_1_m1 = 1.2, end_1_std1 = 0.3,
#'                                               end_1_m2 = 1.45, end_1_std2 = 0.3,
#'                                               end_2_p1 = 0.3, end_2_p2 = 0.5,
#'                         nrep = 300, NREP = 200, sample_size = 100, strategy = "conservative",
#'                         LRV_1 = 0.2, LRV_2 = 0.1, independent = F, correlation = 0.3,
#'                         TV_1 = 0.25, TV_2 = 0.2, FGR = 0.2, FSR = 0.1, seed = 2024)
#'
#' @export
three_outcomes_cont_binary_coprimary <- function(n1_prev, n2_prev,
                                         end_1_m1, end_1_std1,
                                         end_1_m2, end_1_std2,
                                         end_2_p1, end_2_p2,
                                         nrep, NREP, sample_size,
                                         independent = T, correlation = NULL,
                                         strategy = "conservative",
                                         LRV_1, LRV_2,
                                         TV_1, TV_2,
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
    ## Different simulation schemes for correlated or independent endpoints pair
    if (independent == T) {
      ## 1. Generate samples for the population parameters
      control_response_e2 <- round(n1_prev * end_2_p1)
      trt_response_e2 <- round(n2_prev * end_2_p2)
      control_group_pop_samples <- list(end_points_1 = rnorm(nrep, mean = end_1_m1, sd = end_1_std1),
                                        end_points_2 = rbeta(nrep, control_response_e2,
                                                             n1_prev - control_response_e2))
      trt_group_pop_samples <- list(end_points_1 = rnorm(nrep, mean = end_1_m2, sd = end_1_std2),
                                    end_points_2 = rbeta(nrep, trt_response_e2,
                                                         n2_prev - trt_response_e2))

      ## 2. Generate data using the population samples
      control_group_e1_samples <- sapply(control_group_pop_samples[['end_points_1']], function(x) {
        return(rnorm(sample_size, mean = x, sd = end_1_std1 * sqrt(n1_prev)))
      })
      control_group_e2_samples <- sapply(control_group_pop_samples[['end_points_2']], function(x) {
        return(rbinom(1, sample_size, prob = x))
      })
      trt_group_e1_samples <- sapply(trt_group_pop_samples[['end_points_1']], function(x) {
        return(rnorm(sample_size, mean = x, sd = end_1_std2 * sqrt(n2_prev)))
      })
      trt_group_e2_samples <- sapply(trt_group_pop_samples[['end_points_2']], function(x) {
        return(rbinom(1, sample_size, prob = x))
      })

      ## 3. Calculate the confidence intervals for the treatment difference
      confidence_interval <- function(sample_dat_1, sample_dat_2, interval, type) {
        if (type == "continuous") {
          mean_est_1 <- colMeans(sample_dat_1)
          mean_est_2 <- colMeans(sample_dat_2)
          ## CLT to get the confidence intervals
          vec_sd <- sd(mean_est_2 - mean_est_1)
          n <- length(mean_est_2)
          vec_mean <- mean(mean_est_2 - mean_est_1)
          # Error according to normal distribution - should use one sided critical value
          error <- qt(interval, df = n - 1) * vec_sd / sqrt(n)
          # Confidence interval as a vector
          result <- c("lower" = vec_mean - error, "upper" = vec_mean + error)
          return(result)
        } else {
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
        }
      }

      # For endpoints 1
      twenty_percentile_e1 <- confidence_interval(control_group_e1_samples,
                                                  trt_group_e1_samples, interval = 1 - FGR, type = "continuous")
      nighty_percentile_e1 <- confidence_interval(control_group_e1_samples,
                                                  trt_group_e1_samples, interval = 1 - FSR, type = "continuous")

      # For endpoints 2
      twenty_percentile_e2 <- confidence_interval(control_group_e2_samples,
                                                  trt_group_e2_samples, interval = 1 - FGR, type = "binary")
      nighty_percentile_e2 <- confidence_interval(control_group_e2_samples,
                                                  trt_group_e2_samples, interval = 1 - FSR, type = "binary")

      ## 4. Make a decision based on the decision criteria
      # For endpoints 1
      if (nighty_percentile_e1[2] <= TV_1 & twenty_percentile_e1[1] <= LRV_1) {
        e1_decision <- "stop"
      } else if (nighty_percentile_e1[2] > TV_1 & twenty_percentile_e1[1] > LRV_1) {
        e1_decision <- "go"
      } else {
        e1_decision <- "consider"
      }
      # For endpoints 2
      if (nighty_percentile_e2[2] <= TV_2 & twenty_percentile_e2[1] <= LRV_2) {
        e2_decision <- "stop"
      } else if (nighty_percentile_e2[2] > TV_2 & twenty_percentile_e2[1] > LRV_2) {
        e2_decision <- "go"
      } else {
        e2_decision <- "consider"
      }

      ## Different decision criteria
      if (strategy == "conservative") {
        decision_vecs <- c(decision_vecs,
                           ifelse(e1_decision == "go" & e2_decision == "go", "go",
                                  ifelse(e1_decision == "stop" | e2_decision == "stop",
                                         "stop", "consider")))
      } else {
        decision_vecs <- c(decision_vecs,
                           ifelse(e1_decision == "stop" | e2_decision == "stop", "stop",
                                  ifelse(e1_decision == "go" | e2_decision == "go", "go",
                                         "consider")))
      }

      setTxtProgressBar(pb, idx)
      rm(control_group_e1_samples, control_group_e2_samples,
         trt_group_e1_samples, trt_group_e2_samples)
    } else {
      ## 1. Generate samples for the population parameters
      control_response_e2 <- round(n1_prev * end_2_p1)
      trt_response_e2 <- round(n2_prev * end_2_p2)
      control_group_pop_samples <- list(end_points_1 = rnorm(nrep, mean = end_1_m1, sd = end_1_std1),
                                        end_points_2 = rbeta(nrep, control_response_e2,
                                                             n1_prev - control_response_e2))
      trt_group_pop_samples <- list(end_points_1 = rnorm(nrep, mean = end_1_m2, sd = end_1_std2),
                                    end_points_2 = rbeta(nrep, trt_response_e2,
                                                         n2_prev - trt_response_e2))

      ## 2. Generate data using the population samples - correlated
      # Build correaltion matrix
      cor_matrix <- matrix(c(1, correlation,
                             correlation, 1), 2, 2)

      copula_simulation <- function(cor_matrix, distribution = NULL, sample_size,
                                    mean1, sd1, p2) {
        ## 1. Generate a Gaussian copula with the defined correlation
        gaussian_cop <- MASS::mvrnorm(sample_size, mu = c(0, 0), Sigma = cor_matrix)

        ## 2. Simulate correlated uniform variables
        u <- pnorm(gaussian_cop)

        ## 3. Transform uniform variables to desired distributions
        continuous_var <- qnorm(u[, 1], mean = mean1, sd = sd1)
        binary_var <- as.numeric(u[, 2] < p2)

        ## 4. Combine into a data frame and return the results
        simulated_data <- data.frame(continuous_var, binary_var)
        return(simulated_data)
      }

      control_group_samples <- lapply(1:nrep, function(x) {
        tmp_sim <- copula_simulation(cor_matrix = cor_matrix,
                                     sample_size = sample_size,
                                     mean1 = control_group_pop_samples[['end_points_1']][x],
                                     sd1 = end_1_std1 * sqrt(n1_prev),
                                     p2 = control_group_pop_samples[['end_points_2']][x])
        return(colMeans(tmp_sim))
      })

      trt_group_samples <- lapply(1:nrep, function(x) {
        tmp_sim <- copula_simulation(cor_matrix = cor_matrix,
                                     sample_size = sample_size,
                                     mean1 = trt_group_pop_samples[['end_points_1']][x],
                                     sd1 = end_1_std2 * sqrt(n1_prev),
                                     p2 = trt_group_pop_samples[['end_points_2']][x])
        return(colMeans(tmp_sim))
      })

      # Reframe the data structure
      control_group_e1_samples <- do.call(c, lapply(control_group_samples, function(x) {
        return(x[1])
      }))
      control_group_e2_samples <- do.call(c, lapply(control_group_samples, function(x) {
        return(x[2])
      }))
      trt_group_e1_samples <- do.call(c, lapply(trt_group_samples, function(x) {
        return(x[1])
      }))
      trt_group_e2_samples <- do.call(c, lapply(trt_group_samples, function(x) {
        return(x[2])
      }))

      ## 3. Calculate the confidence intervals for the treatment difference
      confidence_interval <- function(sample_dat_1, sample_dat_2, interval) {
        mean_est_1 <- sample_dat_1
        mean_est_2 <- sample_dat_2
        ## CLT to get the confidence intervals
        vec_sd <- sd(mean_est_2 - mean_est_1)
        n <- length(mean_est_2)
        vec_mean <- mean(mean_est_2 - mean_est_1)
        # Error according to normal distribution
        error <- qnorm(interval) * vec_sd / sqrt(n)
        # Confidence interval as a vector
        result <- c("lower" = vec_mean - error, "upper" = vec_mean + error)
        return(result)
      }
      # For endpoints 1
      twenty_percentile_e1 <- confidence_interval(control_group_e1_samples,
                                                  trt_group_e1_samples, interval = 1 - FGR)
      nighty_percentile_e1 <- confidence_interval(control_group_e1_samples,
                                                  trt_group_e1_samples, interval = 1 - FSR)

      # For endpoints 2
      twenty_percentile_e2 <- confidence_interval(control_group_e2_samples,
                                                  trt_group_e2_samples, interval = 1 - FGR)
      nighty_percentile_e2 <- confidence_interval(control_group_e2_samples,
                                                  trt_group_e2_samples, interval = 1 - FSR)

      ## 4. Make a decision based on the decision criteria
      # For endpoints 1
      if (nighty_percentile_e1[2] <= TV_1 & twenty_percentile_e1[1] <= LRV_1) {
        e1_decision <- "stop"
      } else if (nighty_percentile_e1[2] > TV_1 & twenty_percentile_e1[1] > LRV_1) {
        e1_decision <- "go"
      } else {
        e1_decision <- "consider"
      }
      # For endpoints 2
      if (nighty_percentile_e2[2] <= TV_2 & twenty_percentile_e2[1] <= LRV_2) {
        e2_decision <- "stop"
      } else if (nighty_percentile_e2[2] > TV_2 & twenty_percentile_e2[1] > LRV_2) {
        e2_decision <- "go"
      } else {
        e2_decision <- "consider"
      }

      ## Different decision criteria
      if (strategy == "conservative") {
        decision_vecs <- c(decision_vecs,
                           ifelse(e1_decision == "go" & e2_decision == "go", "go",
                                  ifelse(e1_decision == "stop" | e2_decision == "stop",
                                         "stop", "consider")))
      } else {
        decision_vecs <- c(decision_vecs,
                           ifelse(e1_decision == "stop" | e2_decision == "stop", "stop",
                                  ifelse(e1_decision == "go" | e2_decision == "go", "go",
                                         "consider")))
      }

      setTxtProgressBar(pb, idx)
      rm(control_group_samples, trt_group_samples, control_group_e1_samples, control_group_e2_samples,
         trt_group_e1_samples, trt_group_e2_samples)
    }
  }
  ## Return the results with decision proportions
  res <- list(sample_size, decision_vecs)
  names(res) <- c("sample_size", "decision")
  return(res)
}




#' Three outcome approach decision making for continous + count co-primary endpoints
#' @param n1_prev Numeric. Sample size from control group in the previous studies
#' @param n2_prev Numeric. Sample size from treatment group in the previous studies
#' @param end_1_m1 Numeric. Sample mean for endpoints one from control group in the previous studies
#' @param end_1_std1 Numeric. Std for endpoints one from control group in the previous studies
#' @param end_1_m2 Numeric. Sample mean for endpoints one from treatment group in the previous studies
#' @param end_1_std2 Numeric. Std for endpoints one from treatment group in the previous studies
#' @param end_2_shape1 Numeric. Shape parameter to generate Poisson rate for endpoints two from control group in the previous studies
#' @param end_2_rate1 Numeric. Rate parameter to generate Poisson rate for endpoints two from control group in the previous studies
#' @param end_2_shape2 Numeric. Shape parameter to generate Poisson rate for endpoints two from treatment group in the previous studies
#' @param end_2_rate2 Numeric. Rate parameter to generate Poisson rate for endpoints two from treatment group in the previous studies
#' @param nrep Numeric. Replicates to evaluate the decision once
#' @param NREP Numeric. Replicates to calculate the decision proportions
#' @param sample_size Numeric. Target Sample size in the new studies
#' @param independent Logical. TRUE if two endpoints are independent or FALSE if they are correlated
#' @param correlation Numeric. correlation strength between two endpoints if they are correlated
#' @param strategy Character. Two values: conservative or aggressive for decision criteria
#' @param LRV_1 Numeric. Lower reference value (treatment difference) for endpoint one
#' @param TV_1 Numeric. Target value (treatment difference) for endpoint one
#' @param LRV_2 Numeric. Lower reference value (treatment difference) for endpoint two
#' @param TV_2 Numeric. Target value (treatment difference) for endpoint two
#' @param FGR Numeric. False go risk level
#' @param FSR Numeric. False stop risk level
#' @param seed Numeric. Random seed
#' @examples
#' library(TOAST)
#' ## Example code for two binary primary endpoint
#' cont_count_examples_ind <- three_outcomes_cont_count_coprimary(n1_prev = 50, n2_prev = 50,
#'                                               end_1_m1 = 1.2, end_1_std1 = 0.3,
#'                                               end_1_m2 = 1.45, end_1_std2 = 0.3,
#'                                               end_2_shape1 = 3, end_2_rate1 = 5,
#'                                               end_2_shape2 = 5, end_2_rate2 = 6,
#'                         nrep = 300, NREP = 200, sample_size = 100, strategy = "conservative",
#'                         LRV_1 = 0.2, LRV_2 = 0.1, independent = T,
#'                         TV_1 = 0.25, TV_2 = 0.233, FGR = 0.2, FSR = 0.1, seed = 2024)
#' cont_count_examples_cor <- three_outcomes_cont_count_coprimary(n1_prev = 50, n2_prev = 50,
#'                                               end_1_m1 = 1.2, end_1_std1 = 0.3,
#'                                               end_1_m2 = 1.45, end_1_std2 = 0.3,
#'                                               end_2_shape1 = 3, end_2_rate1 = 5,
#'                                               end_2_shape2 = 5, end_2_rate2 = 6,
#'                         nrep = 300, NREP = 200, sample_size = 100, strategy = "conservative",
#'                         LRV_1 = 0.2, LRV_2 = 0.1, independent = F, correlation = 0.3,
#'                         TV_1 = 0.25, TV_2 = 0.233, FGR = 0.2, FSR = 0.1, seed = 2024)
#'
#' @export
three_outcomes_cont_count_coprimary <- function(n1_prev, n2_prev,
                                                 end_1_m1, end_1_std1,
                                                 end_1_m2, end_1_std2,
                                                 end_2_shape1, end_2_rate1,
                                                 end_2_shape2, end_2_rate2,
                                                 nrep, NREP, sample_size,
                                                 independent = T, correlation = NULL,
                                                 strategy = "conservative",
                                                 LRV_1, LRV_2,
                                                 TV_1, TV_2,
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
    ## Different simulation schemes for correlated or independent endpoints pair
    if (independent == T) {
      ## 1. Generate samples for the population parameters
      control_group_pop_samples <- list(end_points_1 = rnorm(nrep, mean = end_1_m1, sd = end_1_std1),
                                        end_points_2 = rgamma(n = nrep,
                                                              shape = end_2_shape1,
                                                              rate = end_2_rate1))
      trt_group_pop_samples <- list(end_points_1 = rnorm(nrep, mean = end_1_m2, sd = end_1_std2),
                                    end_points_2 = rgamma(n = nrep,
                                                          shape = end_2_shape2,
                                                          rate = end_2_rate2))

      ## 2. Generate data using the population samples
      control_group_e1_samples <- sapply(control_group_pop_samples[['end_points_1']], function(x) {
        return(rnorm(sample_size, mean = x, sd = end_1_std1 * sqrt(n1_prev)))
      })
      control_group_e2_samples <- sapply(control_group_pop_samples[['end_points_2']], function(x) {
        return(rpois(n = sample_size, lambda = x))
      })
      trt_group_e1_samples <- sapply(trt_group_pop_samples[['end_points_1']], function(x) {
        return(rnorm(sample_size, mean = x, sd = end_1_std2 * sqrt(n2_prev)))
      })
      trt_group_e2_samples <- sapply(trt_group_pop_samples[['end_points_2']], function(x) {
        return(rpois(n = sample_size, lambda = x))
      })

      ## 3. Calculate the confidence intervals for the treatment difference
      confidence_interval <- function(sample_dat_1, sample_dat_2, interval, type) {
        if (type == "continuous") {
          mean_est_1 <- colMeans(sample_dat_1)
          mean_est_2 <- colMeans(sample_dat_2)
          ## CLT to get the confidence intervals
          vec_sd <- sd(mean_est_2 - mean_est_1)
          n <- length(mean_est_2)
          vec_mean <- mean(mean_est_2 - mean_est_1)
          # Error according to normal distribution - should use one sided critical value
          error <- qt(interval, df = n - 1) * vec_sd / sqrt(n)
          # Confidence interval as a vector
          result <- c("lower" = vec_mean - error, "upper" = vec_mean + error)
          return(result)
        } else {
          mean_est_1 <- colMeans(sample_dat_1)
          mean_est_2 <- colMeans(sample_dat_2)
          ## CLT to get the confidence intervals
          vec_sd <- sd(mean_est_2 - mean_est_1)
          n <- length(mean_est_2)
          vec_mean <- mean(mean_est_2 - mean_est_1)
          # Error according to normal distribution
          error <- qnorm(interval) * vec_sd / sqrt(n)
          # Confidence interval as a vector
          result <- c("lower" = vec_mean - error, "upper" = vec_mean + error)
          return(result)
        }
      }

      # For endpoints 1
      twenty_percentile_e1 <- confidence_interval(control_group_e1_samples,
                                                  trt_group_e1_samples, interval = 1 - FGR, type = "continuous")
      nighty_percentile_e1 <- confidence_interval(control_group_e1_samples,
                                                  trt_group_e1_samples, interval = 1 - FSR, type = "continuous")

      # For endpoints 2
      twenty_percentile_e2 <- confidence_interval(control_group_e2_samples,
                                                  trt_group_e2_samples, interval = 1 - FGR, type = "count")
      nighty_percentile_e2 <- confidence_interval(control_group_e2_samples,
                                                  trt_group_e2_samples, interval = 1 - FSR, type = "count")

      ## 4. Make a decision based on the decision criteria
      # For endpoints 1
      if (nighty_percentile_e1[2] <= TV_1 & twenty_percentile_e1[1] <= LRV_1) {
        e1_decision <- "stop"
      } else if (nighty_percentile_e1[2] > TV_1 & twenty_percentile_e1[1] > LRV_1) {
        e1_decision <- "go"
      } else {
        e1_decision <- "consider"
      }
      # For endpoints 2
      if (nighty_percentile_e2[2] <= TV_2 & twenty_percentile_e2[1] <= LRV_2) {
        e2_decision <- "stop"
      } else if (nighty_percentile_e2[2] > TV_2 & twenty_percentile_e2[1] > LRV_2) {
        e2_decision <- "go"
      } else {
        e2_decision <- "consider"
      }

      ## Different decision criteria
      if (strategy == "conservative") {
        decision_vecs <- c(decision_vecs,
                           ifelse(e1_decision == "go" & e2_decision == "go", "go",
                                  ifelse(e1_decision == "stop" | e2_decision == "stop",
                                         "stop", "consider")))
      } else {
        decision_vecs <- c(decision_vecs,
                           ifelse(e1_decision == "stop" | e2_decision == "stop", "stop",
                                  ifelse(e1_decision == "go" | e2_decision == "go", "go",
                                         "consider")))
      }

      setTxtProgressBar(pb, idx)
      rm(control_group_e1_samples, control_group_e2_samples,
         trt_group_e1_samples, trt_group_e2_samples)
    } else {
      ## 1. Generate samples for the population parameters
      control_group_pop_samples <- list(end_points_1 = rnorm(nrep, mean = end_1_m1, sd = end_1_std1),
                                        end_points_2 = rgamma(n = nrep,
                                                              shape = end_2_shape1,
                                                              rate = end_2_rate1))
      trt_group_pop_samples <- list(end_points_1 = rnorm(nrep, mean = end_1_m2, sd = end_1_std2),
                                    end_points_2 = rgamma(n = nrep,
                                                          shape = end_2_shape2,
                                                          rate = end_2_rate2))

      ## 2. Generate data using the population samples - correlated
      # Build correaltion matrix
      cor_matrix <- matrix(c(1, correlation,
                             correlation, 1), 2, 2)

      copula_simulation <- function(cor_matrix, distribution = NULL, sample_size,
                                    mean1, sd1, rate2) {
        ## 1. Generate a Gaussian copula with the defined correlation
        gaussian_cop <- MASS::mvrnorm(sample_size, mu = c(0, 0), Sigma = cor_matrix)

        ## 2. Simulate correlated uniform variables
        u <- pnorm(gaussian_cop)

        ## 3. Transform uniform variables to desired distributions
        continuous_var <- qnorm(u[, 1], mean = mean1, sd = sd1)
        poisson_var <- qpois(u[, 2], lambda = rate2)

        ## 4. Combine into a data frame and return the results
        simulated_data <- data.frame(continuous_var, poisson_var)
        return(simulated_data)
      }

      control_group_samples <- lapply(1:nrep, function(x) {
        tmp_sim <- copula_simulation(cor_matrix = cor_matrix,
                                     sample_size = sample_size,
                                     mean1 = control_group_pop_samples[['end_points_1']][x],
                                     sd1 = end_1_std1 * sqrt(n1_prev),
                                     rate2 = control_group_pop_samples[['end_points_2']][x])
        return(colMeans(tmp_sim))
      })

      trt_group_samples <- lapply(1:nrep, function(x) {
        tmp_sim <- copula_simulation(cor_matrix = cor_matrix,
                                     sample_size = sample_size,
                                     mean1 = trt_group_pop_samples[['end_points_1']][x],
                                     sd1 = end_1_std2 * sqrt(n1_prev),
                                     rate2 = trt_group_pop_samples[['end_points_2']][x])
        return(colMeans(tmp_sim))
      })

      # Reframe the data structure
      control_group_e1_samples <- do.call(c, lapply(control_group_samples, function(x) {
        return(x[1])
      }))
      control_group_e2_samples <- do.call(c, lapply(control_group_samples, function(x) {
        return(x[2])
      }))
      trt_group_e1_samples <- do.call(c, lapply(trt_group_samples, function(x) {
        return(x[1])
      }))
      trt_group_e2_samples <- do.call(c, lapply(trt_group_samples, function(x) {
        return(x[2])
      }))

      ## 3. Calculate the confidence intervals for the treatment difference
      confidence_interval <- function(sample_dat_1, sample_dat_2, interval) {
        mean_est_1 <- sample_dat_1
        mean_est_2 <- sample_dat_2
        ## CLT to get the confidence intervals
        vec_sd <- sd(mean_est_2 - mean_est_1)
        n <- length(mean_est_2)
        vec_mean <- mean(mean_est_2 - mean_est_1)
        # Error according to normal distribution
        error <- qnorm(interval) * vec_sd / sqrt(n)
        # Confidence interval as a vector
        result <- c("lower" = vec_mean - error, "upper" = vec_mean + error)
        return(result)
      }
      # For endpoints 1
      twenty_percentile_e1 <- confidence_interval(control_group_e1_samples,
                                                  trt_group_e1_samples, interval = 1 - FGR)
      nighty_percentile_e1 <- confidence_interval(control_group_e1_samples,
                                                  trt_group_e1_samples, interval = 1 - FSR)

      # For endpoints 2
      twenty_percentile_e2 <- confidence_interval(control_group_e2_samples,
                                                  trt_group_e2_samples, interval = 1 - FGR)
      nighty_percentile_e2 <- confidence_interval(control_group_e2_samples,
                                                  trt_group_e2_samples, interval = 1 - FSR)

      ## 4. Make a decision based on the decision criteria
      # For endpoints 1
      if (nighty_percentile_e1[2] <= TV_1 & twenty_percentile_e1[1] <= LRV_1) {
        e1_decision <- "stop"
      } else if (nighty_percentile_e1[2] > TV_1 & twenty_percentile_e1[1] > LRV_1) {
        e1_decision <- "go"
      } else {
        e1_decision <- "consider"
      }
      # For endpoints 2
      if (nighty_percentile_e2[2] <= TV_2 & twenty_percentile_e2[1] <= LRV_2) {
        e2_decision <- "stop"
      } else if (nighty_percentile_e2[2] > TV_2 & twenty_percentile_e2[1] > LRV_2) {
        e2_decision <- "go"
      } else {
        e2_decision <- "consider"
      }

      ## Different decision criteria
      if (strategy == "conservative") {
        decision_vecs <- c(decision_vecs,
                           ifelse(e1_decision == "go" & e2_decision == "go", "go",
                                  ifelse(e1_decision == "stop" | e2_decision == "stop",
                                         "stop", "consider")))
      } else {
        decision_vecs <- c(decision_vecs,
                           ifelse(e1_decision == "stop" | e2_decision == "stop", "stop",
                                  ifelse(e1_decision == "go" | e2_decision == "go", "go",
                                         "consider")))
      }

      setTxtProgressBar(pb, idx)
      rm(control_group_samples, trt_group_samples, control_group_e1_samples, control_group_e2_samples,
         trt_group_e1_samples, trt_group_e2_samples)
    }
  }
  ## Return the results with decision proportions
  res <- list(sample_size, decision_vecs)
  names(res) <- c("sample_size", "decision")
  return(res)
}





#' Three outcome approach decision making for continous + count co-primary endpoints
#' @param n1_prev Numeric. Sample size from control group in the previous studies
#' @param n2_prev Numeric. Sample size from treatment group in the previous studies
#' @param end_1_p1 Numeric. Response rate from control group in the previous studies for endpoints one
#' @param end_1_p2 Numeric. Response rate from treatment group in the previous studies for endpoints one
#' @param end_2_shape1 Numeric. Shape parameter to generate Poisson rate for endpoints two from control group in the previous studies
#' @param end_2_rate1 Numeric. Rate parameter to generate Poisson rate for endpoints two from control group in the previous studies
#' @param end_2_shape2 Numeric. Shape parameter to generate Poisson rate for endpoints two from treatment group in the previous studies
#' @param end_2_rate2 Numeric. Rate parameter to generate Poisson rate for endpoints two from treatment group in the previous studies
#' @param nrep Numeric. Replicates to evaluate the decision once
#' @param NREP Numeric. Replicates to calculate the decision proportions
#' @param sample_size Numeric. Target Sample size in the new studies
#' @param independent Logical. TRUE if two endpoints are independent or FALSE if they are correlated
#' @param correlation Numeric. correlation strength between two endpoints if they are correlated
#' @param strategy Character. Two values: conservative or aggressive for decision criteria
#' @param LRV_1 Numeric. Lower reference value (treatment difference) for endpoint one
#' @param TV_1 Numeric. Target value (treatment difference) for endpoint one
#' @param LRV_2 Numeric. Lower reference value (treatment difference) for endpoint two
#' @param TV_2 Numeric. Target value (treatment difference) for endpoint two
#' @param FGR Numeric. False go risk level
#' @param FSR Numeric. False stop risk level
#' @param seed Numeric. Random seed
#' @examples
#' library(TOAST)
#' ## Example code for two binary primary endpoint
#' bin_count_examples_ind <- three_outcomes_bin_count_coprimary(n1_prev = 50, n2_prev = 50,
#'                                               end_1_p1 = 0.3, end_1_p2 = 0.5,
#'                                               end_2_shape1 = 3, end_2_rate1 = 5,
#'                                               end_2_shape2 = 5, end_2_rate2 = 6,
#'                         nrep = 300, NREP = 200, sample_size = 100, strategy = "conservative",
#'                         LRV_1 = 0.1, LRV_2 = 0.1, independent = T,
#'                         TV_1 = 0.2, TV_2 = 0.233, FGR = 0.2, FSR = 0.1, seed = 2024)
#' bin_count_examples_cor <- three_outcomes_bin_count_coprimary(n1_prev = 50, n2_prev = 50,
#'                                               end_1_p1 = 0.3, end_1_p2 = 0.5,
#'                                               end_2_shape1 = 3, end_2_rate1 = 5,
#'                                               end_2_shape2 = 5, end_2_rate2 = 6,
#'                         nrep = 300, NREP = 200, sample_size = 100, strategy = "conservative",
#'                         LRV_1 = 0.1, LRV_2 = 0.1, independent = F, correlation = 0.3,
#'                         TV_1 = 0.2, TV_2 = 0.233, FGR = 0.2, FSR = 0.1, seed = 2024)
#'
#' @export
three_outcomes_bin_count_coprimary <- function(n1_prev, n2_prev,
                                                end_1_p1, end_1_p2,
                                                end_2_shape1, end_2_rate1,
                                                end_2_shape2, end_2_rate2,
                                                nrep, NREP, sample_size,
                                                independent = T, correlation = NULL,
                                                strategy = "conservative",
                                                LRV_1, LRV_2,
                                                TV_1, TV_2,
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
    ## Different simulation schemes for correlated or independent endpoints pair
    if (independent == T) {
      ## 1. Generate samples for the population parameters
      control_response_e1 <- round(n1_prev * end_1_p1)
      trt_response_e1 <- round(n2_prev * end_1_p2)
      control_group_pop_samples <- list(end_points_1 = rbeta(nrep, control_response_e1,
                                                             n1_prev - control_response_e1),
                                        end_points_2 = rgamma(n = nrep,
                                                              shape = end_2_shape1,
                                                              rate = end_2_rate1))
      trt_group_pop_samples <- list(end_points_1 = rbeta(nrep, trt_response_e1,
                                                         n2_prev - trt_response_e1),
                                    end_points_2 = rgamma(n = nrep,
                                                          shape = end_2_shape2,
                                                          rate = end_2_rate2))

      ## 2. Generate data using the population samples
      control_group_e1_samples <- sapply(control_group_pop_samples[['end_points_1']], function(x) {
        return(rbinom(1, sample_size, prob = x))
      })
      control_group_e2_samples <- sapply(control_group_pop_samples[['end_points_2']], function(x) {
        return(rpois(n = sample_size, lambda = x))
      })
      trt_group_e1_samples <- sapply(trt_group_pop_samples[['end_points_1']], function(x) {
        return(rbinom(1, sample_size, prob = x))
      })
      trt_group_e2_samples <- sapply(trt_group_pop_samples[['end_points_2']], function(x) {
        return(rpois(n = sample_size, lambda = x))
      })

      ## 3. Calculate the confidence intervals for the treatment difference
      confidence_interval <- function(sample_dat_1, sample_dat_2, interval, type) {
        if (type == "binary") {
          mean_est_1 <- sample_dat_1/sample_size
          mean_est_2 <- sample_dat_2/sample_size
        } else {
          mean_est_1 <- colMeans(sample_dat_1)
          mean_est_2 <- colMeans(sample_dat_2)
        }
        ## CLT to get the confidence intervals
        vec_sd <- sd(mean_est_2 - mean_est_1)
        n <- length(mean_est_2)
        vec_mean <- mean(mean_est_2 - mean_est_1)
        # Error according to normal distribution
        error <- qnorm(interval) * vec_sd / sqrt(n)
        # Confidence interval as a vector
        result <- c("lower" = vec_mean - error, "upper" = vec_mean + error)
        return(result)
      }

      # For endpoints 1
      twenty_percentile_e1 <- confidence_interval(control_group_e1_samples,
                                                  trt_group_e1_samples, interval = 1 - FGR, type = "binary")
      nighty_percentile_e1 <- confidence_interval(control_group_e1_samples,
                                                  trt_group_e1_samples, interval = 1 - FSR, type = "binary")

      # For endpoints 2
      twenty_percentile_e2 <- confidence_interval(control_group_e2_samples,
                                                  trt_group_e2_samples, interval = 1 - FGR, type = "count")
      nighty_percentile_e2 <- confidence_interval(control_group_e2_samples,
                                                  trt_group_e2_samples, interval = 1 - FSR, type = "count")

      ## 4. Make a decision based on the decision criteria
      # For endpoints 1
      if (nighty_percentile_e1[2] <= TV_1 & twenty_percentile_e1[1] <= LRV_1) {
        e1_decision <- "stop"
      } else if (nighty_percentile_e1[2] > TV_1 & twenty_percentile_e1[1] > LRV_1) {
        e1_decision <- "go"
      } else {
        e1_decision <- "consider"
      }
      # For endpoints 2
      if (nighty_percentile_e2[2] <= TV_2 & twenty_percentile_e2[1] <= LRV_2) {
        e2_decision <- "stop"
      } else if (nighty_percentile_e2[2] > TV_2 & twenty_percentile_e2[1] > LRV_2) {
        e2_decision <- "go"
      } else {
        e2_decision <- "consider"
      }

      ## Different decision criteria
      if (strategy == "conservative") {
        decision_vecs <- c(decision_vecs,
                           ifelse(e1_decision == "go" & e2_decision == "go", "go",
                                  ifelse(e1_decision == "stop" | e2_decision == "stop",
                                         "stop", "consider")))
      } else {
        decision_vecs <- c(decision_vecs,
                           ifelse(e1_decision == "stop" | e2_decision == "stop", "stop",
                                  ifelse(e1_decision == "go" | e2_decision == "go", "go",
                                         "consider")))
      }

      setTxtProgressBar(pb, idx)
      rm(control_group_e1_samples, control_group_e2_samples,
         trt_group_e1_samples, trt_group_e2_samples)
    } else {
      ## 1. Generate samples for the population parameters
      control_response_e1 <- round(n1_prev * end_1_p1)
      trt_response_e1 <- round(n2_prev * end_1_p2)
      control_group_pop_samples <- list(end_points_1 = rbeta(nrep, control_response_e1,
                                                             n1_prev - control_response_e1),
                                        end_points_2 = rgamma(n = nrep,
                                                              shape = end_2_shape1,
                                                              rate = end_2_rate1))
      trt_group_pop_samples <- list(end_points_1 = rbeta(nrep, trt_response_e1,
                                                         n2_prev - trt_response_e1),
                                    end_points_2 = rgamma(n = nrep,
                                                          shape = end_2_shape2,
                                                          rate = end_2_rate2))

      ## 2. Generate data using the population samples - correlated
      # Build correaltion matrix
      cor_matrix <- matrix(c(1, correlation,
                             correlation, 1), 2, 2)

      copula_simulation <- function(cor_matrix, distribution = NULL, sample_size,
                                    p1, rate2) {
        ## 1. Generate a Gaussian copula with the defined correlation
        gaussian_cop <- MASS::mvrnorm(sample_size, mu = c(0, 0), Sigma = cor_matrix)

        ## 2. Simulate correlated uniform variables
        u <- pnorm(gaussian_cop)

        ## 3. Transform uniform variables to desired distributions
        binary_var <- as.numeric(u[, 1] < p1)
        poisson_var <- qpois(u[, 2], lambda = rate2)

        ## 4. Combine into a data frame and return the results
        simulated_data <- data.frame(binary_var, poisson_var)
        return(simulated_data)
      }

      control_group_samples <- lapply(1:nrep, function(x) {
        tmp_sim <- copula_simulation(cor_matrix = cor_matrix,
                                     sample_size = sample_size,
                                     p1 = control_group_pop_samples[['end_points_1']][x],
                                     rate2 = control_group_pop_samples[['end_points_2']][x])
        return(colMeans(tmp_sim))
      })

      trt_group_samples <- lapply(1:nrep, function(x) {
        tmp_sim <- copula_simulation(cor_matrix = cor_matrix,
                                     sample_size = sample_size,
                                     p1 = trt_group_pop_samples[['end_points_1']][x],
                                     rate2 = trt_group_pop_samples[['end_points_2']][x])
        return(colMeans(tmp_sim))
      })

      # Reframe the data structure
      control_group_e1_samples <- do.call(c, lapply(control_group_samples, function(x) {
        return(x[1])
      }))
      control_group_e2_samples <- do.call(c, lapply(control_group_samples, function(x) {
        return(x[2])
      }))
      trt_group_e1_samples <- do.call(c, lapply(trt_group_samples, function(x) {
        return(x[1])
      }))
      trt_group_e2_samples <- do.call(c, lapply(trt_group_samples, function(x) {
        return(x[2])
      }))

      ## 3. Calculate the confidence intervals for the treatment difference
      confidence_interval <- function(sample_dat_1, sample_dat_2, interval) {
        mean_est_1 <- sample_dat_1
        mean_est_2 <- sample_dat_2
        ## CLT to get the confidence intervals
        vec_sd <- sd(mean_est_2 - mean_est_1)
        n <- length(mean_est_2)
        vec_mean <- mean(mean_est_2 - mean_est_1)
        # Error according to normal distribution
        error <- qnorm(interval) * vec_sd / sqrt(n)
        # Confidence interval as a vector
        result <- c("lower" = vec_mean - error, "upper" = vec_mean + error)
        return(result)
      }
      # For endpoints 1
      twenty_percentile_e1 <- confidence_interval(control_group_e1_samples,
                                                  trt_group_e1_samples, interval = 1 - FGR)
      nighty_percentile_e1 <- confidence_interval(control_group_e1_samples,
                                                  trt_group_e1_samples, interval = 1 - FSR)

      # For endpoints 2
      twenty_percentile_e2 <- confidence_interval(control_group_e2_samples,
                                                  trt_group_e2_samples, interval = 1 - FGR)
      nighty_percentile_e2 <- confidence_interval(control_group_e2_samples,
                                                  trt_group_e2_samples, interval = 1 - FSR)

      ## 4. Make a decision based on the decision criteria
      # For endpoints 1
      if (nighty_percentile_e1[2] <= TV_1 & twenty_percentile_e1[1] <= LRV_1) {
        e1_decision <- "stop"
      } else if (nighty_percentile_e1[2] > TV_1 & twenty_percentile_e1[1] > LRV_1) {
        e1_decision <- "go"
      } else {
        e1_decision <- "consider"
      }
      # For endpoints 2
      if (nighty_percentile_e2[2] <= TV_2 & twenty_percentile_e2[1] <= LRV_2) {
        e2_decision <- "stop"
      } else if (nighty_percentile_e2[2] > TV_2 & twenty_percentile_e2[1] > LRV_2) {
        e2_decision <- "go"
      } else {
        e2_decision <- "consider"
      }

      ## Different decision criteria
      if (strategy == "conservative") {
        decision_vecs <- c(decision_vecs,
                           ifelse(e1_decision == "go" & e2_decision == "go", "go",
                                  ifelse(e1_decision == "stop" | e2_decision == "stop",
                                         "stop", "consider")))
      } else {
        decision_vecs <- c(decision_vecs,
                           ifelse(e1_decision == "stop" | e2_decision == "stop", "stop",
                                  ifelse(e1_decision == "go" | e2_decision == "go", "go",
                                         "consider")))
      }

      setTxtProgressBar(pb, idx)
      rm(control_group_samples, trt_group_samples, control_group_e1_samples, control_group_e2_samples,
         trt_group_e1_samples, trt_group_e2_samples)
    }
  }
  ## Return the results with decision proportions
  res <- list(sample_size, decision_vecs)
  names(res) <- c("sample_size", "decision")
  return(res)
}
