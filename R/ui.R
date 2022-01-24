mlemurUI <- function(request) {
  fluidPage(
    shinyjs::useShinyjs(),
    shinyFeedback::useShinyFeedback(),
    rclipboard::rclipboardSetup(),
    tags$head(
      tags$link(rel= "stylesheet", type = "text/css", href = "www/bootstrap.css"),
      tags$script(src = "www/javascript.js"),
      tags$link(rel="shortcut icon", href="www/mlemur.ico")
    ),
    a(name="top"),
    navbarPage(fluid = FALSE,
               id = "navigation",
               title = div(img(src = ("www/mlemur.svg"), style="margin-top:-11px; padding-right:10px; padding-bottom:0px;", class="unselectable", draggable="false", dragstart="false;", height = 45)),
               windowTitle="mlemur: MLE Mutation Rate Calculator",
               footer = column(12,
                               hr(),
                               HTML("<center>mlemur: MLE Mutation Rate Calculator (beta) | 2022 | GPL-2</center>"),
                               br(),
                               br()
               ),
               collapsible = TRUE,
               #### Mutation Rate Tab ####
               tabPanel(
                 "Rate",
                 titlePanel("Mutation rate"),
                 hr(),
                 fluidRow(
                   column(
                     7,
                     settingsPlatingUI("SettingsRate", c(0, 1, 1)),
                     hr(),
                     countsPlatingUI("CountsRate", "SettingsRate", FALSE, usePreset = 1),
                     fluidRow(
                       column(
                         4,
                         shinyFeedback::loadingButton(
                           "calculate",
                           label = HTML("<i class = \"fas fa-calculator\"></i> Calculate"),
                           style = "width:100%;",
                           loadingLabel = "Calculating\U2026"
                         ),
                         br(),
                         br()
                       ),
                       column(
                         4,
                         actionButton("erase",
                                      label = "Clear all cells",
                                      width = "100%",
                                      icon = icon("trash", lib = "font-awesome")
                         ),
                         br(),
                         br()
                       ),
                       column(
                         4,
                         actionButton("sample",
                                      label = "Load sample",
                                      width = "100%",
                                      icon = icon("sync", lib = "font-awesome")
                         ),
                         br(),
                         br()
                       )
                     )
                   ),
                   column(
                     5,
                     h3("Results"),
                     htmlOutput(outputId = "errorBarRate"),
                     tags$style(type = "text/css", "#errorBarRate {white-space: pre-wrap; line-height:18px;}"),
                     shinyjs::hidden(
                       div(
                         id = "advanced",
                         reactable::reactableOutput("tableRate"),
                         br(),
                         fluidRow(
                           column(8,
                                  offset = 2,
                                  div(uiOutput("clip"), style =
                                        "text-align:center")
                           )
                         ),
                         br(),
                         hr(),
                         fluidRow(
                           column(8,
                                  offset = 2,
                                  textInput(
                                    inputId = "datasetName",
                                    label = "Dataset name:",
                                    value = "Strain 1",
                                    width = "100%"
                                  ),
                                  actionButton("appendToReport",
                                               label = "Add to report",
                                               width = "100%",
                                               icon = icon("plus-square", lib = "font-awesome")
                                  ),
                                  br(),
                                  br(),
                                  div(downloadButton("dl", "Export report to XLSX"),
                                      style = "text-align:center")
                           )
                         ),
                         br(),
                         hr(),
                         abbrevsUI("abbrevsrate", FALSE)
                       )
                     )
                   )
                 )
               ),
               #### p-value Tab ####
               tabPanel(
                 title = HTML("<em>P</em> value"),
                 titlePanel(div(HTML(
                   "<em>P</em> value"
                 ))),
                 hr(),
                 fluidRow(
                   column(
                     7,
                     settingsPlatingUI("SettingsPval", c(0, 1, 2)),
                     hr(),
                     fluidRow(
                       column(
                         6,
                         h4(HTML("<b>Strain 1</b>")),
                         countsPlatingUI("CountsStrain1", "SettingsPval", TRUE, usePreset = 1),
                       ),
                       column(
                         6,
                         h4(HTML("<b>Strain 2</b>")),
                         countsPlatingUI("CountsStrain2", "SettingsPval", TRUE, usePreset = 2),
                       )
                     ),
                     br(),
                     fluidRow(
                       column(
                         4,
                         shinyFeedback::loadingButton(
                           "calculate2",
                           label = HTML("<i class = \"fas fa-calculator\"></i> Calculate"),
                           style = "width:100%;",
                           loadingLabel = "Calculating\U2026"
                         ),
                         br(),
                         br()
                       ),
                       column(
                         4,
                         actionButton("erase2",
                                      label = "Clear all cells",
                                      width = "100%",
                                      icon = icon("trash", lib = "font-awesome")
                         ),
                         br(),
                         br()
                       ),
                       column(
                         4,
                         actionButton("sample2",
                                      label = "Load sample",
                                      width = "100%",
                                      icon = icon("sync", lib = "font-awesome")
                         ),
                         br(),
                         br()
                       )
                     )
                   ),
                   column(5,
                          h3("Results"),
                          htmlOutput(outputId = "errorBarPvalue"),
                          tags$style(type = "text/css", "#errorBarPvalue {white-space: pre-wrap; line-height:18px;}"),
                          shinyjs::hidden(
                            div(
                              id = "advanced2",
                              reactable::reactableOutput("tablePvalue"),
                              hr(),
                              htmlOutput("pvalueinfo"),
                              textOutput("effsizeinfo"),
                              br(),
                              hr(),
                              abbrevsUI("abbrevspval", TRUE)
                            )
                          ))
                 )
               ),
               #### Pvalue Correction Tab ####
               tabPanel(
                 title = HTML("<em>P</em> correction"),
                 titlePanel(div(HTML(
                   "<em>P</em> value correction"
                 ))),
                 hr(),
                 fluidRow(
                   column(
                     7,
                     fluidRow(
                       column(
                         6,
                         shinyWidgets::awesomeRadio(
                           inputId = "correctionMethod",
                           label = "Select method of correction:",
                           choices = c("Bonferroni" = "bonferroni", "Bonferroni-Holm" = "holm", "Benjamini-Hochberg" = "BH"),
                           selected = "bonferroni"
                         ),
                         h6()
                       ),
                       column(
                         6,
                         textAreaInput(
                           inputId = "enteredPvalues",
                           label = HTML(paste("Paste <em>P</em> values here:", infoTooltip("Inputs must be within 0 and 1. Scientific notation can be used as needed. Paste one value under the other. Do not separate the values with a comma, semicolon, or any other character."))),
                           value = "0.064\n0.00007\n0.002",
                           rows = 10,
                           width = "100%",
                           resize = "vertical"
                         ),
                         br()
                       )
                     ),
                     fluidRow(
                       column(
                         4,
                         actionButton(
                           "calculate3",
                           label = "Calculate",
                           width = "100%",
                           class = "btn btn-primary",
                           icon = icon("calculator", lib = "font-awesome")
                         ),
                         br(),
                         br()
                       ),
                       column(
                         4,
                         actionButton("erase3",
                                      label = "Clear all cells",
                                      width = "100%",
                                      icon = icon("trash", lib = "font-awesome")
                         ),
                         br(),
                         br()
                       ),
                       column(
                         4,
                         actionButton("sample3",
                                      label = "Load sample",
                                      width = "100%",
                                      icon = icon("sync", lib = "font-awesome")
                         ),
                         br(),
                         br()
                       )
                     )
                   ),
                   column(5,
                          h3("Results"),
                          htmlOutput(outputId = "errorBarCorrector"),
                          tags$style(type = "text/css", "#errorBarCorrector {white-space:pre-wrap; line-height:18px;}"),
                          shinyjs::hidden(
                            div(
                              id = "advanced3",
                              reactable::reactableOutput("tableCorrector"),
                              br(),
                              fluidRow(
                                column(8,
                                       offset = 2,
                                       div(uiOutput("clip2"), style =
                                             "text-align:center")
                                )
                              )
                            )
                          ))
                 )
               ),
               #### Batch Tab ####
               tabPanel("Batch",
                        titlePanel(div(HTML(
                          "Batch mutation rate and <em>P</em> value"
                        ))),
                        hr(),
                        "XLS(X) files only! For instructions, check Help, or see example.",
                        a(href="example.xlsx", HTML("<i class=\"fa fa-download\" aria-hidden=\"true\"></i> Download example with comments"), download=NA, target="_blank"),
                        a(href="template.xlsx", HTML("<i class=\"fa fa-download\" aria-hidden=\"true\"></i> Download empty template to fill"), download=NA, target="_blank"),
                        fluidRow(
                          column(4,
                                 offset = 4,
                                 h3("Step 1: Upload file"),
                                 fileInput(inputId = "userData", width = "100%", label = "Select XLS/XLSX file", accept = c("application/vnd.ms-excel", "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet", ".xls", ".xlsx"), multiple = FALSE),
                                 htmlOutput("errorBarFile"),
                                 shinyjs::hidden(
                                   div(
                                     id = "onDataUploaded",
                                     h3("Step 2: Load data"),
                                     shinyFeedback::loadingButton(
                                       "dataCheck",
                                       class = "btn btn-default",
                                       label = "Load & Check data",
                                       style = "width:100%",
                                       loadingLabel = "Checking\U2026"
                                     )
                                   )
                                 )
                          )
                        ),
                        br(),
                        htmlOutput(outputId = "batchInfo"),
                        tags$style(type = "text/css", "#batchInfo {white-space: pre-wrap; line-height:18px;}"),
                        br(),
                        shinyjs::hidden(
                          div(
                            id = "onDataChecked",
                            tabsetPanel(
                              tabPanel(title = "Experiment parameters",
                                       id = "selectorPlating",
                                       reactable::reactableOutput("platingTable")
                              ),
                              tabPanel(title = "Counts on selective medium",
                                       id = "selectorSel",
                                       reactable::reactableOutput("selectiveTable")
                              ),
                              tabPanel(title = "Counts on non-selective medium",
                                       id = "selectorNsel",
                                       reactable::reactableOutput("nonselectiveTable")
                              )
                            ),
                            shinyjs::hidden(
                              div(
                                id = "ifDataCorrect",
                                fluidRow(
                                  column(4,
                                         offset = 4,
                                         h3("Step 3: Select options"),
                                         uiOutput("BatchOptions", inline = FALSE),
                                         h3("Step 4: Calculate"),
                                         shinyFeedback::loadingButton(
                                           "calculate4",
                                           label = HTML("<i class = \"fas fa-calculator\"></i> Calculate"),
                                           style = "width:100%;",
                                           loadingLabel = "Calculating\U2026"
                                         ),
                                         br()
                                  )
                                )
                              )
                            ),
                            shinyjs::hidden(
                              div(
                                id = "onOutputGenerated",
                                br(),
                                textOutput("infoBatch"),
                                HTML("<h4>Mutation Rates</h4>"),
                                reactable::reactableOutput("batchRateResults"),
                                br(),
                                shinyjs::hidden(
                                  div(
                                    id = "BatchPvalue",
                                    HTML("<h4><i>P</i> values</h4>"),
                                    reactable::reactableOutput("BatchPvalueResults"),
                                    br()
                                  )
                                ),
                                fluidRow(
                                  column(4,
                                         offset = 4,
                                         h3("Step 5: Download results"),
                                         div(downloadButton("dlBatch", "Export the results to XLSX"), style =
                                               "text-align:center; width:100%")
                                  )
                                )
                              )
                            )
                          )
                        ),
                        br()
            #             #### Help Tab ####
            #             tabPanel("Help",
            #                      titlePanel("Help"),
            #                      hr(),
            #                      withMathJax(),
            #                      tags$div(HTML("<script type='text/x-mathjax-config'>
            # MathJax.Hub.Config({
            # tex2jax: {inlineMath: [['$','$']]}
            # });
            # </script>
            # ")),
            # div(
            #   includeHTML(system.file("www", "help.txt", package = "mlemur"))
            # )
            #             )
               )
    )
  )
}