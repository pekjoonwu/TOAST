#' @import shiny
#' @import shinydashboard
#' @import shinythemes
## UI design
ui <- shinydashboard::dashboardPage(
  skin = "blue",  # Change the skin color
  shinydashboard::dashboardHeader(title = "TOAST"),

  shinydashboard::dashboardSidebar(
    shinydashboard::sidebarMenu(
      shinydashboard::menuItem("One Time Evaluation", tabName = "one_time_eval", icon = shiny::icon("dashboard")),
      shinydashboard::menuItem("Evaluation of Parameter Ranges", tabName = "param_range_eval", icon = shiny::icon("th"))
    )
  ),

  shinydashboard::dashboardBody(
    ## shinythemes::themeSelector(),
    theme = shinythemes::shinytheme("spacelab"),
    shinydashboard::tabItems(
      shinydashboard::tabItem(tabName = "one_time_eval",
              shiny::fluidRow(
                shinydashboard::box(
                  title = "Evaluation Parameters", status = "primary", solidHeader = TRUE, collapsible = TRUE,
                  shiny::selectInput("dist_type", "Select Endpoints Type:",
                              choices = c("Continuous" = "normal", "Binary" = "binary", "Count" = "count")),
                  shiny::conditionalPanel(
                    condition = "input.dist_type == 'normal'",
                    shiny::numericInput("n1prev", "Previous Sample Size (Control):", value = 50),
                    shiny::numericInput("n2prev", "Previous Sample Size (Treatment):", value = 50),
                    shiny::numericInput("mean1", "Population Mean (Control):", value = 1.2),
                    shiny::numericInput("std1", "Population SD (Control):", value = 0.3, min = 0),
                    shiny::numericInput("mean2", "Population Mean (Treatment):", value = 1.4),
                    shiny::numericInput("std2", "Population SD (Treatment):", value = 0.3, min = 0)
                  ),
                  shiny::conditionalPanel(
                    condition = "input.dist_type == 'binary'",
                    shiny::numericInput("n1prev", "Previous Sample Size (Control):", value = 50),
                    shiny::numericInput("n2prev", "Previous Sample Size (Treatment):", value = 50),
                    shiny::numericInput("p1", "Success Probability (Control):", value = 0.3, min = 0, max = 1, step = 0.01),
                    shiny::numericInput("p2", "Success Probability (Treatment):", value = 0.5, min = 0, max = 1, step = 0.01)
                  ),
                  shiny::conditionalPanel(
                    condition = "input.dist_type == 'count'",
                    shiny::numericInput("shape1", "Shape Parameter of Poisson Rate (Control):", value = 3, min = 0.1),
                    shiny::numericInput("rate1", "Rate Parameter of Poisson Rate (Control):", value = 5, min = 0.1),
                    shiny::numericInput("shape2", "Shape Parameter of Poisson Rate (Treatment):", value = 4, min = 0.1),
                    shiny::numericInput("rate2", "Rate Parameter of Poisson Rate (Treatment):", value = 5, min = 0.1)
                  ),
                  shiny::numericInput("sample_size", "Sample Size:", value = 100, min = 10),
                  shiny::numericInput("LRV", "Lower Reference Value:", value = 0.1, min = 0),
                  shiny::numericInput("TV", "Target Value:", value = 0.2, min = 0),
                  shiny::actionButton("simulate", "Simulate!")
                ),
                shinydashboard::box(
                  title = "Results", status = "primary", solidHeader = TRUE, collapsible = TRUE,
                  shiny::plotOutput("lineplot", width = "800px", height = "600px")
                )
              )
      ),

      shinydashboard::tabItem(tabName = "param_range_eval",
              shiny::fluidRow(
                shinydashboard::box(
                  title = "Parameter Range Evaluation", status = "primary", solidHeader = TRUE, collapsible = TRUE,
                  shiny::selectInput("range_dist_type", "Select Endpoints Type:",
                              choices = c("Continuous" = "normal", "Binary" = "binary", "Count" = "count")),
                  shiny::uiOutput("param_ui"),
                  shiny::numericInput("range_sample_size", "Sample Size for Evaluation", value = 100),
                  shiny::actionButton("range_simulate", "Simulate for Evaluation!")
                ),
                shinydashboard::box(
                  title = "Range Evaluation Results", status = "primary", solidHeader = TRUE, collapsible = TRUE,
                  shiny::plotOutput("rangeLineplot", width = "800px", height = "600px")
                )
              )
      )
    )
  )
)

