#' Run the TOAST Shiny Application
#' @importFrom shiny shinyApp runApp
#' @export
#' @examples
#'
#'\dontrun{
#'  ## Will Load an Interactive Session
#' shinyOutput <- TOAST_shiny()
#'}

TOAST_shiny <- function() {
  out <- shiny::runApp(shiny::shinyApp(ui, server))
}
