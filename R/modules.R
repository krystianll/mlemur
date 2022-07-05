settingsPlatingUI <- function(id, defaultSettings) {
  ns <- NS(id)
  tagList(
    fluidRow(
      column(
        width = 12,
        shinyWidgets::awesomeRadio(
          inputId = ns("setSelective"),
          label = HTML(paste("Plating efficiency:", infoTooltip("Choose whether you want to calculate the fraction of culture plated on selective medium (plating efficiency) using appropriate volumes, or already have a pre-calculated value."))),
          choices = c("Calculate using provided volumes and dilutions" = 1, "Use a pre-calculated value" = 2),
          selected = defaultSettings[2],
          inline = TRUE
        ),
        shinyWidgets::awesomeRadio(
          inputId = ns("setNonselective"),
          label = HTML(paste("Average number of cells in the culture:", infoTooltip("Choose whether you want to calculate average number of cells in culture (culture size) using colony counts and appropriate volumes, or already have a pre-calculated value."))),
          choices = c("Calculate using the colony counts and volumes" = 1, "Use a pre-calculated value" = 2),
          selected = defaultSettings[3],
          inline = TRUE
        ),
        shinyWidgets::awesomeCheckbox(
          inputId = ns("setModel"),
          label = HTML(paste("Use advanced experiment parameters", infoTooltip("A more detailed description of available mutant models can be found in Help."))),
          value = FALSE
        )
      )
    )
  )
}

settingsPlating <- function(input, output, session) {
  rv <- reactiveValues()
  observe({
    rv$value <- list("setModel" = input$setModel,
                     "setNonselective" = input$setNonselective,
                     "setSelective" = input$setSelective
    )
  })
  return(reactive(rv$value))
}

countsPlatingUI <- function(id, currentTab, stack_cols=FALSE, usePreset) {
  ns <- NS(id)
  number <- ifelse(stack_cols, 12, 6)
  if (usePreset==1) {useValues <- c("500", "100", "1", "1\n2\n3", "100", "1000000", "100\n200\n300", "0.9", "0.2", "1e+9", "0", "0", "0", "0", "0", "0")}
  else if (usePreset) {useValues <- c("500", "150", "1", "7\n8\n9", "100", "1000000", "90\n100\n110", "1", "0.3", "5e+8", "0", "0", "0", "0", "0", "0")}
  tagList(
    shinyFeedback::useShinyFeedback(),
    conditionalPanel(
      condition="input.setModel",
      ns = NS(currentTab),
      fluidRow(
        column(
          width = 12,
          shinyWidgets::awesomeRadio(
            inputId = ns("useLagFitness"),
            label = HTML(paste("Specify phenotypic lag or fitness:", infoTooltip("Choose whether you want to input phenotypic lag or mutant fitness."))),
            choices = c("Don't use" = 0, "Specify phenotypic lag" = 1, "Specify mutant fitnes" = 2),
            selected = 0,
            inline = TRUE
          )
        ),
        column(
          width = number,
          conditionalPanel(
            condition="input.useLagFitness==1",
            ns = NS(id),
            textInput(
              inputId = ns("Lag"),
              label = HTML(paste("Phenotypic lag (generations):", infoTooltip("Please provide a positive number."))),
              value = useValues[12],
              width = "100%"
            )
          ),
          conditionalPanel(
            condition="input.useLagFitness==2",
            ns = NS(id),
            textInput(
              inputId = ns("Fitness"),
              label = HTML(paste("Mutant relative fitness:", infoTooltip("Please provide a positive number. Fitness &lt;&nbsp;1 means that mutants grow <u>slower</u> than non-mutants, while &gt;&nbsp;1 means that mutants grow faster. Fitness equal to 1 is equivalent to L-C model.<br><b>Note: underflow may be encountered when when both plating efficiency and fitness are much smaller than 0.1.</b>"))),
              value = useValues[8],
              width = "100%"
            )
          )
        ),
        column(
          width = 12,
          textInput(
            inputId = ns("Death"),
            label = HTML(paste("Death rate:", infoTooltip("Please provide a non-negative number that is a fraction of the birth rate (&ge;&nbsp;0 but &lt;&nbsp;1)."))),
            value = useValues[14],
            width = "100%"
          )
        ),
        column(
          width = number,
          textInput(
            inputId = ns("Residual"),
            label = HTML(paste("Residual mutations:", infoTooltip("Please provide a non-negative number. Floating-point numbers will be rounded down to nearest integer. Residual mutations apply only to the portion of culture that was plated."))),
            value = useValues[13],
            width = "100%"
          )
        ),
        column(
          width = number,
          textInput(
            inputId = ns("Inoculum"),
            label = HTML(paste("Size of the inoculum (number of cells):", infoTooltip("Please provide a non-negative number. Floating-point numbers will be rounded down to nearest integer. Inoculum size cannot exceed the total number of cells in grown culture."))),
            value = useValues[16],
            width = "100%"
          )
        ),
        column(
          width = 12,
          div(
            style="display:inline-block; vertical-align:top;",
            shinyWidgets::awesomeRadio(
              inputId = ns("setCV"),
              label = HTML(paste("Coefficient of variation:", infoTooltip("Use if you want to account for the variation of the final number of cells in each culture. If using custom value, please provide a positive number. Scientific notation can be used, e.g. 2e-2. Values smaller than 1e-6 will be set to 0, and values bigger than 10 will be set to 10. Leave 0 if you do not want to use the coefficient of variation. When only one colony count or only the final number of cells is provided, it will default to 0 unless supplied with another value."))),
              choices = c("Calculate using the colony counts" = 0, "Use custom value:" = 1),
              selected = 1,
              inline = TRUE
            )
          ),
          div(
            style = "display:inline-block; vertical-align:bottom;",
            textInput(
              inputId = ns("CV"),
              label = NULL,
              value = useValues[11],
              width = "100%",
            )
          )
        )
      ),
      hr()
    ),
    conditionalPanel(
      condition="input.setSelective==1 || input.setNonselective==1",
      ns = NS(currentTab),
      textInput(
        inputId = ns("VolumeTotal"),
        label = HTML(paste("Total culture volume (&mu;l):", infoTooltip("Please provide a positive number."))),
        value = useValues[1],
        width = "100%"
      )
    ),
    fluidRow(
      column(
        width = number,
        h4("Selective medium:"),
        conditionalPanel(
          condition="input.setSelective==1",
          ns = NS(currentTab),
          textInput(
            inputId = ns("VolumeSelective"),
            label = HTML(paste("Volume plated on selective medium (&mu;l):", infoTooltip("Please provide a positive number. Please note that this value, divided by the dilution factor, cannot exceed total culture volume."))),
            value = useValues[2],
            width = "100%"
          ),
          textInput(
            inputId = ns("DilutionSelective"),
            label = HTML(paste("Dilution factor:", infoTooltip("Please provide a positive number &gt;&nbsp;1 if culture was diluted prior plating. If culture was concentrated, use a number &lt;&nbsp;1. If culture was neither diluted nor concentrated, type 1."))),"",
            value = useValues[3],
            width = "100%"
          )
        ),
        conditionalPanel(
          condition="input.setSelective==2",
          ns = NS(currentTab),
          textInput(
            inputId = ns("PlatingEfficiency"),
            label = HTML(paste("Plating efficiency:", infoTooltip("Please provide a positive number so that 0&nbsp;&lt;&nbsp;&epsilon;&nbsp;&le;&nbsp;1."))),
            value = useValues[9],
            width = "100%"
          ),
        ),
        textAreaInput(
          inputId = ns("CountsSelective"),
          label = HTML(paste("Colony counts:", infoTooltip("Zeros and positive numbers are accepted. At least two numbers must be provided, and at least one must be bigger than 0. Paste one number under the other. Do not separate the numbers with a comma, semicolon, or any other character."))),
          value = useValues[4],
          rows = 10,
          width = "100%",
          resize = "vertical"
        ),
        br()
      ),
      column(
        width = number,
        h4("Non-selective medium:"),
        conditionalPanel(
          condition = "input.setNonselective==1",
          ns = NS(currentTab),
          textInput(
            inputId = ns("VolumeNonselective"),
            label = HTML(paste("Volume plated on non-selective medium (&mu;l):", infoTooltip("Please provide a positive number. Please note that this value, divided by the dilution factor, cannot exceed total culture volume."))),
            value = useValues[5],
            width = "100%"
          ),
          textInput(
            inputId = ns("DilutionNonselective"),
            label = HTML(paste("Dilution factor:", infoTooltip("Please provide a positive number &gt;&nbsp;1 if culture was diluted prior plating. If culture was concentrated, use a number &lt;&nbsp;1. If culture was neither diluted nor concentrated, type 1."))),"",
            value = useValues[6],
            width = "100%"
          ),
          textAreaInput(
            inputId = ns("CountsNonselective"),
            label = HTML(paste("Colony counts:", infoTooltip("Zeros and positive numbers are accepted. At least one number must be provided if you have chosen not to account for the variation of the culture size (L-C or M-K models), and at least two if chosen otherwise (B<sup>0</sup> model). At least one colony count must be bigger than 0. Paste one number under the other. Do not separate the numbers with a comma, semicolon, or any other character."))),
            value = useValues[7],
            rows = 10,
            width = "100%",
            resize = "vertical"
          )
        ),
        conditionalPanel(
          condition = "input.setNonselective==2",
          ns = NS(currentTab),
          textInput(
            inputId = ns("MeanCells"),
            label = HTML(paste("Average number of cells in culture:", infoTooltip("Please provide a positive number. Scientific notation can be used, e.g. 4e9."))),
            value = useValues[10],
            width = "100%"
          )
        )
      )
    )
  )
}

countsPlating <- function(input, output, session, userSettings, stack_cols) {
  ReactValue <- reactiveValues()
  
  observeEvent(input$useLagFitness,{
    if (!stack_cols) {
      if (input$useLagFitness!=0) {
        shinyjs::runjs("document.getElementById('CountsRate-Death-label').parentElement.parentElement.setAttribute('class','col-sm-6')")
      } else {
        shinyjs::runjs("document.getElementById('CountsRate-Death-label').parentElement.parentElement.setAttribute('class','col-sm-12')")
      }
    }
  })
  
  observe({
    if (input$setCV==0) {
      shinyjs::disable(id = "CV")
    } else {
      shinyjs::enable(id = "CV")
    }
    
    ReactValue$setNonselective <- userSettings()$setNonselective
    ReactValue$setSelective <- userSettings()$setSelective
    ReactValue$setModel <- userSettings()$setModel
    ReactValue$setCV <- input$setCV
    ReactValue$useLagFitness <- input$useLagFitness
    
    if (ReactValue$setNonselective==1 | ReactValue$setSelective==1) {
      ReactValue$VolumeTotal <- numerise(input$VolumeTotal)
    } else {
      ReactValue$VolumeTotal <- NA
    }
    if (ReactValue$setSelective==1) {
      ReactValue$VolumeSelective <- numerise(input$VolumeSelective)
      ReactValue$DilutionSelective <- numerise(input$DilutionSelective)
      ReactValue$PlatingEfficiency <- NA
    } else if (ReactValue$setSelective==2) {
      ReactValue$VolumeSelective <- NA
      ReactValue$DilutionSelective <- NA
      ReactValue$PlatingEfficiency <- numerise(input$PlatingEfficiency)
    }
    ReactValue$CountsSelective <- numerise(input$CountsSelective)
    if (ReactValue$setModel) {
      if (ReactValue$useLagFitness == 1) {
        ReactValue$Lag <- numerise(input$Lag)
      } else {
        ReactValue$Lag <- NA
      }
      if (ReactValue$useLagFitness == 2) {
        ReactValue$Fitness <- numerise(input$Fitness)
      } else {
        ReactValue$Fitness <- NA
      }
      ReactValue$Death <- numerise(input$Death)
      ReactValue$Residual <- numerise(input$Residual)
      ReactValue$Inoculum <- numerise(input$Inoculum)
      if (ReactValue$setCV==1) {
        ReactValue$CV <- numerise(input$CV)
      }
    } else {
      ReactValue$Fitness <- NA
      ReactValue$Lag <- NA
      ReactValue$Residual <- NA
      ReactValue$Death <- NA
      ReactValue$Inoculum <- NA
      ReactValue$CV <- NA
      shinyFeedback::hideFeedback(inputId = "CV")
    }
    if (ReactValue$setNonselective==1) {
      ReactValue$VolumeNonselective <- numerise(input$VolumeNonselective)
      ReactValue$DilutionNonselective <- numerise(input$DilutionNonselective)
      ReactValue$CountsNonselective <- numerise(input$CountsNonselective)
      ReactValue$MeanCells <- NA
    } else {
      ReactValue$VolumeNonselective <- NA
      ReactValue$DilutionNonselective <- NA
      ReactValue$CountsNonselective <- NA
      ReactValue$MeanCells <- numerise(input$MeanCells)
    }
    
    ReactValue$errorsDetected <- FALSE
    
    if (ReactValue$setNonselective == 1 || ReactValue$setSelective == 1) {
      if (ValueValidator(ReactValue$VolumeTotal) != "") {
        ReactValue$errorsDetected <- TRUE
        textInputError(inputId = "VolumeTotal", text = paste(ValueValidator(ReactValue$VolumeTotal)))
      } else {
        shinyFeedback::hideFeedback(inputId = "VolumeTotal")
      }
    }
    
    if (ReactValue$setModel) {
      
      if (ReactValue$useLagFitness == 1) {
        if (NonNegValueValidator(ReactValue$Lag) != "") {
          ReactValue$errorsDetected <- TRUE
          textInputError(inputId = "Lag", text = paste(NonNegValueValidator(ReactValue$Lag)))
        } else {
          shinyFeedback::hideFeedback(inputId = "Lag")
        }
      }
      
      if (ReactValue$useLagFitness == 2) {
        if (ValueValidator(ReactValue$Fitness) != "") {
          ReactValue$errorsDetected <- TRUE
          textInputError(inputId = "Fitness", text = paste(ValueValidator(ReactValue$Fitness)))
        } else {
          shinyFeedback::hideFeedback(inputId = "Fitness")
        }
      }
      
      if (NonNegValueValidator(ReactValue$Residual) != "") {
        ReactValue$errorsDetected <- TRUE
        textInputError(inputId = "Residual", text = paste(NonNegValueValidator(ReactValue$Residual)))
      } else {
        shinyFeedback::hideFeedback(inputId = "Residual")
      }
      
      if (NonNegValueValidator(ReactValue$Inoculum) != "") {
        ReactValue$errorsDetected <- TRUE
        textInputError(inputId = "Inoculum", text = paste(NonNegValueValidator(ReactValue$Inoculum)))
      } else {
        shinyFeedback::hideFeedback(inputId = "Inoculum")
      }
      
      if (DeathValidator(ReactValue$Death) != "") {
        ReactValue$errorsDetected <- TRUE
        textInputError(inputId = "Death", text = paste(DeathValidator(ReactValue$Death)))
      } else {
        shinyFeedback::hideFeedback(inputId = "Death")
      }
      
      if (ReactValue$setCV == 1) {
        if (NonNegValueValidator(ReactValue$CV) != "") {
          ReactValue$errorsDetected <- TRUE
          textInputError(inputId = "CV", text = paste(NonNegValueValidator(ReactValue$CV)))
        } else {
          shinyFeedback::hideFeedback(inputId = "CV")
        }
      }
    }
    
    if (ReactValue$setSelective == 1) {
      if (ValueValidator(ReactValue$VolumeSelective) != "") {
        ReactValue$errorsDetected <- TRUE
        textInputError(inputId = "VolumeSelective", text = paste(ValueValidator(ReactValue$VolumeSelective)))
      } else {
        shinyFeedback::hideFeedback(inputId = "VolumeSelective")
      }
      
      if (ValueValidator(ReactValue$DilutionSelective) != "") {
        ReactValue$errorsDetected <- TRUE
        textInputError(inputId = "DilutionSelective", text = paste(ValueValidator(ReactValue$DilutionSelective)))
      } else {
        shinyFeedback::hideFeedback(inputId = "DilutionSelective")
      }
      
    } else if (ReactValue$setSelective == 2) {
      if (PlatEffValidator(ReactValue$PlatingEfficiency) != "") {
        ReactValue$errorsDetected <- TRUE
        textInputError(inputId = "PlatingEfficiency", text = paste(PlatEffValidator(ReactValue$PlatingEfficiency)))
      } else {
        shinyFeedback::hideFeedback(inputId = "PlatingEfficiency")
      }
    }
    
    if (SelectiveValidator(ReactValue$CountsSelective) != "") {
      ReactValue$errorsDetected <- TRUE
      textInputError(inputId = "CountsSelective", text = paste(SelectiveValidator(ReactValue$CountsSelective)))
    } else {
      shinyFeedback::hideFeedback(inputId = "CountsSelective")
    }
    
    if (ReactValue$setNonselective==1) {
      if (ValueValidator(ReactValue$VolumeNonselective) != "") {
        ReactValue$errorsDetected <- TRUE
        textInputError(inputId = "VolumeNonselective", text = paste(ValueValidator(ReactValue$VolumeNonselective)))
      } else {
        shinyFeedback::hideFeedback(inputId = "VolumeNonselective")
      }
      
      if (ValueValidator(ReactValue$DilutionNonselective) != "") {
        ReactValue$errorsDetected <- TRUE
        textInputError(inputId = "DilutionNonselective", text = paste(ValueValidator(ReactValue$DilutionNonselective)))
      } else {
        shinyFeedback::hideFeedback(inputId = "DilutionNonselective")
      }
      
      if (NonselectiveValidator(ReactValue$CountsNonselective) != "") {
        ReactValue$errorsDetected <- TRUE
        textInputError(inputId = "CountsNonselective", text = paste(NonselectiveValidator(ReactValue$CountsNonselective)))
      } else {
        shinyFeedback::hideFeedback(inputId = "CountsNonselective")
      }
      
    } else if (ReactValue$setNonselective==2) {
      if (ValueValidator(ReactValue$MeanCells) != "") {
        ReactValue$errorsDetected <- TRUE
        textInputError(inputId = "MeanCells", text = paste(ValueValidator(ReactValue$MeanCells)))
      } else {
        shinyFeedback::hideFeedback(inputId = "MeanCells")
      }
    }
    
    if (ReactValue$setSelective==1) {
      if (paste(ValueValidator(ReactValue$VolumeTotal), ValueValidator(ReactValue$VolumeSelective), ValueValidator(ReactValue$DilutionSelective), sep = "") == "") {
        if (PlatingValidator(ReactValue$VolumeSelective, ReactValue$DilutionSelective, ReactValue$VolumeTotal) != "") {
          ReactValue$errorsDetected <- TRUE
          textInputError(inputId = "VolumeSelective", text = paste(PlatingValidator(ReactValue$VolumeSelective, ReactValue$DilutionSelective, ReactValue$VolumeTotal)))
          textInputError(inputId = "DilutionSelective", text = "")
          textInputError(inputId = "VolumeTotal", text = "")
        } else {
          shinyFeedback::hideFeedback(inputId = "VolumeSelective")
          shinyFeedback::hideFeedback(inputId = "DilutionSelective")
          shinyFeedback::hideFeedback(inputId = "VolumeTotal")
        }
      }
    }
    
    if (ReactValue$setNonselective==1) {
      if (paste(ValueValidator(ReactValue$VolumeTotal), ValueValidator(ReactValue$VolumeNonselective), ValueValidator(ReactValue$DilutionNonselective), sep = "") == "") {
        if (PlatingValidator(ReactValue$VolumeNonselective, ReactValue$DilutionNonselective, ReactValue$VolumeTotal) != "") {
          ReactValue$errorsDetected <- TRUE
          textInputError(inputId = "VolumeNonselective", text = paste(PlatingValidator(ReactValue$VolumeNonselective, ReactValue$DilutionNonselective, ReactValue$VolumeTotal)))
          textInputError(inputId = "DilutionNonselective", text = "")
          textInputError(inputId = "VolumeTotal", text = "")
        } else {
          shinyFeedback::hideFeedback(inputId = "VolumeNonselective")
          shinyFeedback::hideFeedback(inputId = "DilutionNonselective")
          shinyFeedback::hideFeedback(inputId = "VolumeTotal")
        }
      }
    }
    
    if (ReactValue$setModel) {
      if (ReactValue$setNonselective == 2) {
        if (paste(NonNegValueValidator(ReactValue$Inoculum), ValueValidator(ReactValue$MeanCells), sep = "") == "") {
          if (ReactValue$Inoculum >= ReactValue$MeanCells) {
            ReactValue$errorsDetected <- TRUE
            textInputError(inputId = "Inoculum", text = paste("Size of the inoculum must be smaller than mean culture size."))
            textInputError(inputId = "MeanCells", text = "")
          } else {
            shinyFeedback::hideFeedback(inputId = "Inoculum")
            shinyFeedback::hideFeedback(inputId = "MeanCells")
          }
        }
      } else if (ReactValue$setNonselective == 1) {
        if (paste(NonNegValueValidator(ReactValue$Inoculum), NonselectiveValidator(ReactValue$CountsNonselective), ValueValidator(ReactValue$VolumeTotal), ValueValidator(ReactValue$VolumeNonselective), ValueValidator(ReactValue$DilutionNonselective), sep = "") == "") {
          if (ReactValue$Inoculum >= (mean(ReactValue$CountsNonselective) * ReactValue$DilutionNonselective * ReactValue$VolumeNonselective / ReactValue$VolumeTotal)) {
            ReactValue$errorsDetected <- TRUE
            textInputError(inputId = "Inoculum", text = paste("Size of the inoculum must be smaller than mean culture size."))
            textInputError(inputId = "CountsNonselective", text = "")
            textInputError(inputId = "DilutionNonselective", text = "")
            textInputError(inputId = "VolumeNonselective", text = "")
            textInputError(inputId = "VolumeTotal", text = "")
          } else {
            shinyFeedback::hideFeedback(inputId = "Inoculum")
            shinyFeedback::hideFeedback(inputId = "CountsNonselective")
            shinyFeedback::hideFeedback(inputId = "DilutionNonselective")
            shinyFeedback::hideFeedback(inputId = "VolumeNonselective")
            shinyFeedback::hideFeedback(inputId = "VolumeTotal")
          }
        }
      }
    }
  })
  return(reactive(list("VolumeTotal" = ReactValue$VolumeTotal, "Fitness" = ReactValue$Fitness, "VolumeSelective" = ReactValue$VolumeSelective, "DilutionSelective" = ReactValue$DilutionSelective,
                       "PlatingEfficiency" = ReactValue$PlatingEfficiency, "CountsSelective" = ReactValue$CountsSelective, "VolumeNonselective" = ReactValue$VolumeNonselective,
                       "DilutionNonselective" = ReactValue$DilutionNonselective, "CountsNonselective" = ReactValue$CountsNonselective, "MeanCells" = ReactValue$MeanCells, "CV" = ReactValue$CV,
                       "Lag" = ReactValue$Lag, "Residual" = ReactValue$Residual, "Death" = ReactValue$Death, "Inoculum" = ReactValue$Inoculum,
                       "model" = ReactValue$setModel, "errors" = ReactValue$errorsDetected, "setSel" = ReactValue$setSelective, "setNsel" = ReactValue$setNonselective, "setCV" = ReactValue$setCV)))
}

countsPlatingUpdate <- function(input, output, session, usePreset) {
  inputsVec <- c("VolumeTotal", "VolumeSelective", "DilutionSelective", "CountsSelective", "VolumeNonselective", "DilutionNonselective",
                 "CountsNonselective", "Fitness", "PlatingEfficiency", "MeanCells", "CV", "Lag", "Residual", "Death", "Inoculum")
  if (usePreset==0) {useValues <- rep("", 16)}
  else if (usePreset==1) {useValues <- c("500", "100", "1", "1\n2\n3", "100", "1000000", "100\n200\n300", "0.9", "0.2", "1e+9", "0", "0", "0", "0", "0")}
  else if (usePreset==2) {useValues <- c("500", "150", "1", "7\n8\n9", "100", "1000000", "90\n100\n110", "1", "0.3", "5e+8", "0", "0", "0", "0", "0")}
  else {return(NULL)}
  for (i in 1:16) {
    updateTextInput(session, inputId = inputsVec[i], value = useValues[i])
  }
}

abbrevsUI <- function(id, alpha=TRUE) {
  ns <- NS(id)
  tagList(
    h5("Abbreviations:"),
    h6(HTML("&epsilon; &mdash; plating efficiency")),
    h6(HTML("N<sub>t</sub> &mdash; total number of cells in the culture")),
    h6(HTML("CV &mdash; coefficient of variation of N<sub>t</sub>")),
    h6(HTML("&rho; &mdash; relative fitness of the mutant cells")),
    h6(HTML("&lambda; &mdash; mean phenotypic lag (generations)")),
    h6(HTML("d &mdash; death rate (as a fraction of growth rate)")),
    h6(HTML("m<sub>p</sub> &mdash; number of residual mutations")),
    h6(HTML("&phi; &mdash; size of inoculum relative to final culture size")),
    h6(HTML("m &mdash; number of mutations per culture")),
    h6(HTML("m<sup>95&percnt;&ndash;</sup> &mdash; lower limit of 95&percnt; CI for m")),
    h6(HTML("m<sup>95&percnt;&plus;</sup> &mdash; upper limit of 95&percnt; CI for m")),
    h6(HTML("&mu; &mdash; mutation rate per cell per generation")),
    h6(HTML("&mu;<sup>95&percnt;&ndash;</sup> &mdash; lower limit of 95&percnt; CI for &mu;")),
    h6(HTML("&mu;<sup>95&percnt;&plus;</sup> &mdash; upper limit of 95&percnt; CI for &mu;")),
    if (alpha==TRUE) {h6(HTML("&alpha; &mdash; significance level"))}
  )
}