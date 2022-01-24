#' Run mlemur
#'
#' This function initializes mlemur.
#'
#' @export
mlemur <- function(options=list("display.mode" = "normal", "launch.browser" = T)) {
  shiny::shinyApp(ui=mlemurUI, server=mlemurServer, options=options)
}
# mlemurGUI <- function() {
#   shiny::shinyApp(ui=mlemurUI, server=mlemurServer, options=list("display.mode" = "normal", "launch.browser" = T))
# }