#' @import ggplot2
#' @import shiny
## Server function to execute
server <- function(input, output, session) {
  #####
  ## Single Evaluation
  #####
  ## Reactive programming for single parameter
  run_threeoutcome <- shiny::eventReactive(input$simulate, {
    endpoints_type <- input$dist_type
    if (endpoints_type == "normal") {
      three_outcomes_cont(n1_prev = input$n1prev, n2_prev = input$n2prev,
                          mean_1 = input$mean1, mean_2 = input$mean2,
                          std_1 = input$std1, std_2 = input$std2,
                          Method = "empirical",
                          nrep = 300, NREP = 300, sample_size = input$sample_size,
                          LRV = input$LRV, TV = input$TV, seed = 2024)
    } else if (endpoints_type == "binary") {
      three_outcomes_binary(n1_prev = input$n1prev, n2_prev = input$n2prev,
                            p_1 = input$p1, p_2 = input$p2,
                            Method = "empirical",
                            nrep = 300, NREP = 300, sample_size = input$sample_size,
                            LRV = input$LRV, TV = input$TV, seed = 2024)
    } else if (endpoints_type == "count") {
      three_outcomes_count(shape_1 = input$shape1, shape_2 = input$shape2,
                           rate_1 = input$rate1, rate_2 = input$rate2,
                           Method = "empirical",
                           nrep = 300, NREP = 300, sample_size = input$sample_size,
                           LRV = input$LRV, TV = input$TV, seed = 2024)
    }
  })

  ## Plot the line plot to show the results
  output$lineplot <- shiny::renderPlot({
    power <- c(table(run_threeoutcome()$decision)["go"],
               table(run_threeoutcome()$decision)["consider"],
               table(run_threeoutcome()$decision)["stop"]) / 300
    power[is.na(power)] <- 0
    tmp_dat <- data.frame(Outcome = c("Go", "Consider", "Stop"),
                          power = as.numeric(power))
    tmp_dat$Outcome <- factor(tmp_dat$Outcome,
                              levels = c("Go", "Consider", "Stop"))
    ggplot2::ggplot(tmp_dat, ggplot2::aes(x = "", y = power, fill = Outcome)) +
      ggplot2::geom_bar(stat = "identity", width = 1) +
      ggplot2::geom_text(ggplot2::aes(label = paste0(round(power, 2) * 100, "%")),
                position = ggplot2::position_stack(vjust = 0.5)) + ggplot2::coord_polar("y", start = 0) +
      ggplot2::scale_fill_manual(values = c("#52b788", "#ffe97f", "#ff758f")) + ggplot2::theme_bw() +
      ggplot2::ggtitle("Proportions for Three Outcomes") + ggplot2::xlab("") + ggplot2::ylab("") +
      ggplot2::theme(legend.text = ggplot2::element_text(size = 15),
            legend.title = ggplot2::element_text(size = 15),
            title = ggplot2::element_text(size = 17))
  }, res = 96)

  #####
  ## Range Evaluation
  #####
  ## Design the output UI according to the type of endpoints users choose
  output$param_ui <- shiny::renderUI({
    req(input$range_dist_type)
    dist <- input$range_dist_type
    ## Different UIs
    if (dist == "normal") {
      shiny::tagList(
        shiny::selectInput("param_to_range", "Parameter to Set Range",
                    choices = c("Population SD (Shared between Control and Treatment)" = "std",
                                "LRV/TV Ratio" = "lrvtv_ratio")),
        shiny::numericInput("range_min", "Minimum Value for Selected Parameter", value = 0.1),
        shiny::numericInput("range_max", "Maximum Value for Selected Parameter", value = 1)
      )
    } else if (dist == "binary") {
      shiny::tagList(
        shiny::selectInput("param_to_range", "Parameter to Set Range",
                    choices = c("LRV/TV Ratio" = "lrvtv_ratio")),
        shiny::numericInput("range_min", "Minimum Value for Selected Parameter", value = 0.1),
        shiny::numericInput("range_max", "Maximum Value for Selected Parameter", value = 0.99)
      )
    } else if (dist == "count") {
      shiny::tagList(
        shiny::selectInput("param_to_range", "Parameter to Set Range",
                    choices = c("LRV/TV Ratio" = "lrvtv_ratio")),
        shiny::numericInput("range_min", "Minimum Value for Selected Parameter", value = 0.1),
        shiny::numericInput("range_max", "Maximum Value for Selected Parameter", value = 0.99)
      )
    }
  })

  ## Reactive programming for range parameter
  run_range_threeoutcome <- shiny::eventReactive(input$range_simulate, {
    dist_type <- input$range_dist_type
    param_to_range <- input$param_to_range
    range_min <- input$range_min
    range_max <- input$range_max
    range_sample_size <- input$range_sample_size
    ## The range values to run
    param_values <- seq(range_min, range_max, length.out = 5)

    result_list <- lapply(param_values, function(val) {
      if (dist_type == "normal") {
        if (param_to_range == "std") {
          args <- list(n1_prev = 50, n2_prev = 50, mean_1 = 1.2, mean_2 = 1.4, std_1 = val, std_2 = val)
          do.call(three_outcomes_cont, c(args, list(Method = "empirical", nrep = 300, NREP = 300,
                                                    sample_size = range_sample_size, LRV = 0.1, TV = 0.2, seed = 2024)))
        } else if (param_to_range == "lrvtv_ratio") {
          ## Fix TV for now and seek alternative solutions
          args <- list(n1_prev = 50, n2_prev = 50, mean_1 = 1.2, mean_2 = 1.4, LRV = 0.2 * val)
          do.call(three_outcomes_cont, c(args, list(Method = "empirical",
                                                    std_1 = 0.3, std_2 = 0.3,
                                                    nrep = 300, NREP = 300,
                                                    sample_size = range_sample_size,
                                                    TV = 0.2, seed = 2024)))
        }
      } else if (dist_type == "binary") {
        args <- list(n1_prev = 50, n2_prev = 50, LRV = 0.2 * val)
        do.call(three_outcomes_binary, c(args, list(Method = "empirical",
                                                    p_1 = 0.3, p_2 = 0.5,
                                                    nrep = 300, NREP = 300,
                                                    sample_size = range_sample_size,
                                                    TV = 0.2, seed = 2024)))
      } else if (dist_type == "count") {
        args <- list(LRV = 0.2 * val)
        do.call(three_outcomes_count, c(args, list(Method = "empirical",
                                                   shape_1 = 3, shape_2 = 4,
                                                   rate_1 = 5, rate_2 = 5,
                                                   nrep = 300, NREP = 300,
                                                   sample_size = range_sample_size,
                                                   TV = 0.2, seed = 2024)))
      }
    })
    result_list
  })

  ## Reactive expression for x-axis label
  x_axis_label <- shiny::reactive({
    if (input$param_to_range == "std") {
      "Population SD"
    } else if (input$param_to_range == "lrvtv_ratio") {
      "LRV/TV Ratio"
    }
  })

  ## Plot the results for the range evaluation
  output$rangeLineplot <- shiny::renderPlot({
    result_list <- run_range_threeoutcome()

    param_values <- seq(input$range_min, input$range_max, length.out = 5)
    power_list <- lapply(result_list, function(res) {
      # Initialize power vector with zeros for each possible outcome
      power <- c(go = 0, consider = 0, stop = 0)
      # Compute the power values only for the outcomes that are present
      outcome_counts <- table(res$decision)
      power[names(outcome_counts)] <- outcome_counts / 300
      power
    })

    tmp_dat <- do.call(rbind, lapply(1:length(param_values), function(i) {
      power_values <- power_list[[i]]
      power_values[is.na(power_values)] <- 0
      data.frame(Parameter = param_values[i], Outcome = c("Go", "Consider", "Stop"), power = power_values)
    }))
    tmp_dat$Outcome <- factor(tmp_dat$Outcome, levels = c("Go", "Consider", "Stop"))

    ggplot2::ggplot(tmp_dat, ggplot2::aes(x = Parameter, y = power, color = Outcome)) +
      ggplot2::geom_line(ggplot2::aes(linetype = Outcome), linewidth = 1) +
      ggplot2::geom_point(size = 3) +
      ggplot2::theme_bw() + ggplot2::scale_color_manual(values = c("#52b788", "#ffe97f", "#ff758f")) +
      ggplot2::ggtitle("Proportions for Three Outcomes across Parameter Range") + ggplot2::xlab(x_axis_label()) + ggplot2::ylab("Proportions") +
      ggplot2::theme(legend.text = ggplot2::element_text(size = 15),
            axis.title = ggplot2::element_text(size = 15),
            axis.text = ggplot2::element_text(size = 15),
            legend.title = ggplot2::element_text(size = 15),
            title = ggplot2::element_text(size = 17))
    # plotly::plot_ly(tmp_dat, x = ~Parameter, y = ~power, color = ~Outcome,
    #                 type = "scatter", mode = "lines+markers")
  }, res = 96)
}
