mlemurServer <- function(input, output, session) {
  rv <- reactiveValues()
  
  #### Rate Validate Data ####
  rv$SettingsRateUserInput <- callModule(settingsPlating, id = "SettingsRate")
  rv$CountsRateUserInput <- callModule(countsPlating, id = "CountsRate", userSettings=reactive(rv$SettingsRateUserInput()))
  
  #### Rate Calculate ####
  observeEvent(input$calculate, {
    userData <- rv$CountsRateUserInput()
    if (userData$errors == TRUE) {
      output$errorBarRate <- renderText("There are problems with your data that need to be resolved. If you have trouble, check Help.")
      shinyjs::hide("advanced")
      shinyFeedback::resetLoadingButton("calculate")
    } else {
      outputData <- calc.rate.int(userData)[1:15]
      results <- data.frame("Calculations" = outputData,
                            row.names = c("&epsilon;", "N<sub>t</sub>", "CV", "&rho;", "&Lambda;", "&delta;<sub>1</sub>", "&delta;<sub>2</sub>", "m<sub>p</sub>", "&phi;",
                                          "m", "m<sup>95&percnt;&ndash;</sup>", "m<sup>95&percnt;&plus;</sup>", "&mu;", "&mu;<sup>95&percnt;&ndash;</sup>", "&mu;<sup>95&percnt;&plus;</sup>"),
                            check.names = FALSE)
      output$tableRate <- reactable::renderReactable({
        reactable::reactable(results,
                             rownames = TRUE,
                             sortable = FALSE,
                             pagination = FALSE,
                             defaultColDef = reactable::colDef(html = TRUE))
      })
      rv$RatetoClipboard <- outputData
      output$errorBarRate <- renderText("")
      shinyjs::show("advanced")
      shinyjs::enable("calculate")
      shinyjs::runjs("$('#calculate').html(\"<i class = 'fas fa-calculator'></i> Calculate\")")
    }
    gc()
  })
  
  #### Rate Copy to Clipboard ####
  output$clip <- renderUI({
    rclipboard::rclipButton(
      inputId = "clipbtn",
      label = "Copy to Clipboard",
      clipText = paste(as.character(rv$RatetoClipboard), collapse = "\n"),
      modal = FALSE,
      icon = icon("clipboard"),
      width = "100%"
    )
  })
  
  #### Rate Add to Report ####
  rv$Ratereport <- data.frame(
    "Parameters" = c("eff", "Nt", "CV", "Fitness", "Lag", "WT Death", "Mut Death", "Residual m", "Rel size of inoc", "m", "m_low", "m_up", "mu", "mu_low", "mu_up")
  )
  
  observeEvent(input$appendToReport, {
    addedSet <- data.frame(c(rv$RatetoClipboard), check.names = FALSE)
    names(addedSet) <- c(input$datasetName)
    rv$Ratereport <- cbind(rv$Ratereport, addedSet)
    shinyFeedback::showToast(type = "success", message = "Appended succesfully", .options = list(positionClass = "toast-bottom-right", progressBar = FALSE, timeOut = 1500))
  })
  
  #### Rate Download Report ####
  output$dl <- downloadHandler(
    filename = function() {
      paste("Report-", Sys.Date(), ".xlsx", sep = "")
    },
    content = function(file) {
      writexl::write_xlsx(rv$Ratereport, path = file)
    }
  )
  
  #### Rate Erase ####
  observeEvent(input$erase, {
    callModule(countsPlatingUpdate, "CountsRate", usePreset = 0)
  })
  
  #### Rate Load Sample ####
  observeEvent(input$sample, {
    callModule(countsPlatingUpdate, "CountsRate", usePreset = 1)
  })
  
  #### P value Validate Data ####
  rv$SettingsPvalUserInput <- callModule(settingsPlating, id = "SettingsPval")
  rv$CountsStrain1UserInput <- callModule(countsPlating, id = "CountsStrain1", userSettings=reactive(rv$SettingsPvalUserInput()))
  rv$CountsStrain2UserInput <- callModule(countsPlating, id = "CountsStrain2", userSettings=reactive(rv$SettingsPvalUserInput()))
  
  #### P value Calculate ####
  observeEvent(input$calculate2, {
    userDataStrain1 <- rv$CountsStrain1UserInput()
    userDataStrain2 <- rv$CountsStrain2UserInput()
    if (userDataStrain1$errors == TRUE || userDataStrain2$errors == TRUE) {
      output$errorBarPvalue <- renderText("There are problems with your data that need to be resolved. If you have trouble, check Help.")
      shinyjs::hide("advanced2")
      shinyFeedback::resetLoadingButton("calculate2")
    } else {
      Strain1 <- calc.rate.int(userDataStrain1)
      Strain2 <- calc.rate.int(userDataStrain2)
      results2 <- data.frame(
        "Strain 1" = Strain1[1:15],
        "Strain 2" = Strain2[1:15],
        row.names = c("&epsilon;", "N<sub>t</sub>", "CV", "&rho;", "&Lambda;", "&delta;<sub>1</sub>", "&delta;<sub>2</sub>", "m<sub>p</sub>", "&phi;",
                      "m", "m<sup>95&percnt;&ndash;</sup>", "m<sup>95&percnt;&plus;</sup>", "&mu;", "&mu;<sup>95&percnt;&ndash;</sup>", "&mu;<sup>95&percnt;&plus;</sup>"),
        check.names = FALSE
      )
      output$tablePvalue <- reactable::renderReactable({
        reactable::reactable(
          results2,
          rownames = TRUE,
          sortable = FALSE,
          wrap = FALSE,
          pagination = FALSE,
          defaultColDef = reactable::colDef(html = TRUE))
      })
      statistics <- calc.pval.int(userDataStrain1, userDataStrain2, as.numeric(Strain1[16]), as.numeric(Strain2[16]))
      pvalue <- statistics[1]
      output$pvalueinfo <- renderText({
        if (is.na(pvalue) == TRUE | is.nan(pvalue) == TRUE) {
          Outp <- paste("impossible to calculate using provided datasets.")
        }
        else if (pvalue >= 0.05) {
          Outp <- paste( signif(pvalue, 4), "&mdash; the difference is not significant at &alpha;&equals;0.05.")
        }
        else {
          Outp <- paste( signif(pvalue, 4), "&mdash; the difference is significant at &alpha;&equals;0.05.")
        }
        paste("<em>", "P", "</em> value is ", Outp)
      })
      output$errorBarPvalue <- renderText("")
      shinyjs::show("advanced2")
      shinyjs::enable("calculate2")
      shinyjs::runjs("$('#calculate2').html(\"<i class = 'fas fa-calculator'></i> Calculate\")")
    }
    gc()
  })
  
  #### P value Erase ####
  observeEvent(input$erase2, {
    callModule(countsPlatingUpdate, "CountsStrain1", usePreset = 0)
    callModule(countsPlatingUpdate, "CountsStrain2", usePreset = 0)
  })
  
  #### P value Load Sample ####
  observeEvent(input$sample2, {
    callModule(countsPlatingUpdate, "CountsStrain1", usePreset = 1)
    callModule(countsPlatingUpdate, "CountsStrain2", usePreset = 2)
  })
  
  #### Corrector Data Validation ####
  observe({
    rv$enteredPvalues <- numerise(input$enteredPvalues)
    rv$selectedMethod <- input$correctionMethod
    
    rv$errorCorrector <- 0
    
    if (PvalueValidator(rv$enteredPvalues, rv$selectedMethod) != "") {
      rv$errorCorrector <- 1
      textInputError(inputId = "enteredPvalues", text = paste(PvalueValidator(rv$enteredPvalues, rv$selectedMethod)))
    } else {
      shinyFeedback::hideFeedback(inputId = "enteredPvalues")
    }
    
  })
  
  #### Corrector Error Bar ####
  observeEvent(input$calculate3, {
    if (rv$errorCorrector == 1) {
      output$errorBarCorrector <- renderText("There are problems with your data that need to be resolved. If you have trouble, check Help.")
      shinyjs::hide("advanced3")
      #### Corrector Calculator ####
    } else {
      output$errorBarCorrector <- renderText("")
      correctedPvalues <- p.adjust(p=rv$enteredPvalues, method=rv$selectedMethod)
      results3 <-
        data.frame(
          "Entered <em>P</em> values" = as.character(rv$enteredPvalues),
          "Corrected <em>P</em> values" = as.character(correctedPvalues),
          check.names = FALSE
        )
      rv$oldPvalues <- c("Entered P values", as.character(rv$enteredPvalues))
      rv$newPvalues <- c("Corrected P values", as.character(correctedPvalues))
      output$tableCorrector <- reactable::renderReactable({
        reactable::reactable(
          results3,
          sortable = FALSE,
          wrap = TRUE,
          defaultColDef = reactable::colDef(html = TRUE)
        )
      })
      shinyjs::show("advanced3")
    }
  })
  
  output$clip2 <- renderUI({
    rclipboard::rclipButton(
      inputId = "clipbtn",
      label = "Copy to Clipboard",
      clipText = paste(rv$oldPvalues, rv$newPvalues, sep = "\t", collapse = "\n"),
      modal = FALSE,
      icon = icon("clipboard")
    )
  })
  
  #### Corrector Eraser ####
  observeEvent(input$erase3, {
    updateTextAreaInput(session, inputId = "enteredPvalues", value = "")
  })
  
  #### Corrector Sample ####
  observeEvent(input$sample3, {
    updateTextAreaInput(session, inputId = "enteredPvalues", value = "0.064\n0.00007\n0.002")
  })
  
  #### BatchCalc file upload ####
  observeEvent(input$userData, {
    output$batchInfo <- renderText("")
    shinyjs::hide("onDataChecked")
    shinyjs::hide("ifDataCorrect")
    shinyjs::hide("onOutputGenerated")
    shinyjs::runjs("$('a[data-value=\"Experiment parameters\"]').click()")
    
    userInput <- tryCatch(readxl::read_excel(input$userData$datapath, sheet = 1, col_types = "text"), error = function(err) {NULL})
    if (is.null(userInput)) {
      output$errorBarFile <- renderText("<p class=\"text-danger\">Invalid file. Please make sure you uploaded the correct file.</p>")
      shinyjs::hide("onDataUploaded")
      return(NULL)
    }
    sheet_number <- length(readxl::excel_sheets(input$userData$datapath))
    if (sheet_number != 1) {
      text1 <- "This file has more than 1 sheets. Only first sheet was loaded.\n"
    } else {text1 <- ""}
    userInput <- userInput[!apply(is.na(userInput), 1, all),!apply(is.na(userInput), 2, all)]
    text2 <- ""
    if (is.na(any(as.character(userInput[[1]])=="CountsSelective", na.rm = TRUE))) {
      output$errorBarFile <- renderText("<p class=\"text-danger\">Invalid file. Please make sure you uploaded the correct file.</p>")
      shinyjs::hide("onDataUploaded")
      return(NULL)
    } else if (!any(as.character(userInput[[1]])=="CountsSelective", na.rm = TRUE)) {
      text2 <- "CountsSelective was not found in the first column, therefore mlemur assumed that you provided only colony counts on selective medium."
      userInput <- as.data.frame(cbind(rep("CountsSelective", nrow(userInput)), as.matrix(userInput)))
    }
    for (i in 2:nrow(userInput)) {
      if (is.na(userInput[[(i-1),1]])==TRUE) {userInput[[(i-1),1]] <- "NoName"}
      if (is.na(userInput[[(i),1]])==TRUE) {userInput[[i,1]] <- userInput[[(i-1),1]]}
    }
    rv$DataSelective <- subset(userInput, subset = userInput[[1]]=='CountsSelective')
    rv$DataNonselective <- subset(userInput, subset = userInput[[1]]=='CountsNonselective')
    rv$DataPlating <- subset(userInput, subset = userInput[[1]]!='CountsSelective' & userInput[[1]]!='CountsNonselective')
    output$errorBarFile <- renderText(paste(text1, text2, sep = ""))
    shinyjs::show("onDataUploaded")
  })
  
  #### BatchCalc data checker ####
  observeEvent(input$dataCheck, {
    
    validation <- comboValidator(rv$DataPlating, rv$DataSelective, rv$DataNonselective)
    
    if (is.null(validation)==FALSE) {
      
      DataPlating <- validation[[1]]
      DataSelective <- validation[[2]]
      DataNonselective <- validation[[3]]
      
      rv$DataPlatingForCalc <- validation[[4]]
      rv$DataSelectiveForCalc <- validation[[5]]
      rv$DataNonselectiveForCalc <- validation[[6]]
      
      warningList <- validation[[7]]
      errorList <- validation[[8]]
      
      adjustTabNames(validation[[9]], validation[[10]], validation[[11]])
      
      rv$CVFromCounts <- validation[[12]]
      rv$PvaluePossible <- validation[[13]]
      
      if (warningList != "" & errorList != "") {
        output$batchInfo <- renderText({paste("<p class=\"text-danger\">There are problems with your data that need to be resolved. Look for red triangles for more details:</p>", errorList, "\n",
                                              "<p>Additionally, mlemur has generated the following non-error informations and warnings:</p>", warningList, sep = "")})
        shinyjs::hide("ifDataCorrect")
      } else if (errorList != "") {
        output$batchInfo <- renderText({paste("<p class=\"text-danger\">There are problems with your data that need to be resolved. Look for red triangles for more details:</p>", errorList, sep = "")})
        shinyjs::hide("ifDataCorrect")
      } else {
        if (warningList != "") {
          output$batchInfo <- renderText({paste("<p class=\"text-success\">Your data appears to be OK.</p>", "\n",
                                                "<p>Additionally, mlemur has generated the following non-error informations and warnings:</p>", warningList, sep = "")})
        } else {
          output$batchInfo <- renderText({paste("<p class=\"text-success\">Your data appears to be OK.</p>")})
        }
        shinyjs::show("ifDataCorrect")
        
        output$BatchOptions <- renderUI({
          tagList(
            shinyWidgets::awesomeCheckboxGroup(
              inputId = "batchSetParameters",
              label = "Select parameters to use", 
              choices = c("Coefficient of variation", "Mutant Fitness", "Phenotypic Lag", "Residual mutations", "Wild-type Cell Death Rate", "Mutant Cell Death Rate", "Inoculum size"),
              selected = c("Coefficient of variation", "Mutant Fitness", "Phenotypic Lag", "Residual mutations", "Wild-type Cell Death Rate", "Mutant Cell Death Rate", "Inoculum size")
            ),
            shinyjs::disabled(div(class = "checkboxdiv",
                                  shinyWidgets::awesomeCheckbox("BatchCalculateRates",
                                                                HTML("Mutation rates are always calculated."),
                                                                status = "primary",
                                                                value = TRUE))),
            if (rv$PvaluePossible)
            {tagList(div(class = "checkboxdiv",
                         shinyWidgets::awesomeCheckbox("BatchCalculatePvalues",
                                                       HTML("Calculate <em>P</em> values"),
                                                       status = "primary")),
                     conditionalPanel(condition = "input.BatchCalculatePvalues",
                                      div(class = "checkboxdiv",
                                          shinyWidgets::awesomeCheckbox("batchApplyCorrection",
                                                                        HTML("Apply a <em>P</em> value correction"),
                                                                        status = "primary")),
                                      conditionalPanel(condition = "input.batchApplyCorrection",
                                                       div(class = "checkboxdiv",
                                                           style = "padding-left:27px; padding-bottom:5px;",
                                                           shinyWidgets::awesomeRadio("batchCorrectionMethod",
                                                                                      label = HTML(paste("Select method of correction:", infoTooltip("A more detailed description of available correction methods can be found in Help."))),
                                                                                      choices = c("Bonferroni" = "bonferroni", "Bonferroni-Holm" = "holm", "Benjamini-Hochberg" = "BH"),
                                                                                      selected = "bonferroni",
                                                                                      width = "100%") )
                                      )
                     )
            )}
          )
        })
      }
      
      output$platingTable <- reactable::renderReactable({defaultTable(DataPlating, "Input", TRUE)})
      output$selectiveTable <- reactable::renderReactable({defaultTable(DataSelective, "No")})
      output$nonselectiveTable <- reactable::renderReactable({defaultTable(DataNonselective, "No")})
      
      shinyjs::show("onDataChecked")
      
    } else {
      shinyjs::hide("ifDataCorrect")
      output$batchInfo <- renderText("<p class=\"text-danger\">There are no data in this file.</p>")
    }
    
    shinyFeedback::resetLoadingButton("dataCheck")
    
  })
  
  #### BatchCalc calculations ####
  observeEvent(input$calculate4, {
    t1 <- Sys.time()
    rv$PlatingList <- lapply(rv$DataPlatingForCalc, as.list)
    for (item in colnames(rv$DataPlatingForCalc)) {
      eval(parse(text = paste("names(rv$PlatingList$`", item, "`) <- rownames(rv$DataPlatingForCalc)", sep = "")))
      eval(parse(text = paste("rv$PlatingList$`", item, "`$CountsSelective <- rv$DataSelectiveForCalc$`", item, "`", sep = "")))
      eval(parse(text = paste("rv$PlatingList$`", item, "`$CountsNonselective <- rv$DataNonselectiveForCalc$`", item, "`", sep = "")))
      eval(parse(text = paste("rv$PlatingList$`", item, "`$model <- TRUE", sep = "")))
      eval(parse(text = paste("rv$PlatingList$`", item, "`$setCV <- !rv$CVFromCounts[which(colnames(rv$DataPlatingForCalc)=='", item, "')]", sep = "")))
    }
    rv$PlatingList <- removeParameters(rv$PlatingList, input$batchSetParameters)

    DataRate <- lapply(rv$PlatingList, calc.rate.int)
    NoOfMutations <- as.numeric(lapply(DataRate, function(x) {x[16]}))
    DataRate <- as.data.frame(lapply(DataRate, function(x) {x[-16]}), check.names = FALSE)
    
    rownames(DataRate) <- c("&epsilon;", "N<sub>t</sub>", "CV", "&rho;", "&Lambda;", "&delta;<sub>1</sub>", "&delta;<sub>2</sub>", "m<sub>p</sub>", "&phi;",
                            "m", "m<sup>95&percnt;&ndash;</sup>", "m<sup>95&percnt;&plus;</sup>", "&mu;", "&mu;<sup>95&percnt;&ndash;</sup>", "&mu;<sup>95&percnt;&plus;</sup>")
    output$batchRateResults <- reactable::renderReactable({defaultTable(DataRate, "Parameters")})
    rv$RateOutput <- data.frame("Parameters" = c("eff", "Nt", "CV", "Fitness", "Lag", "WTDeath", "MutantDeath", "Residual", "Inoculum", "m", "m_low", "m_up", "mu", "mu_low", "mu_up"), DataRate, check.names = FALSE)

    if (input$BatchCalculatePvalues) {
      
      table.length <- length(DataRate)
      xarg <- rep(1:(table.length-1), (table.length-1):1)
      yarg <- unlist(sapply(2:table.length, function(x) x:table.length))
      
      proxy <- mapply(function(x,y) {calc.pval.int(rv$PlatingList[[x]], rv$PlatingList[[y]], NoOfMutations[[x]], NoOfMutations[[y]])}, x = xarg, y = yarg)
      
      DataPvalue <- matrix(NA,nrow=table.length,ncol=table.length,dimnames = list(names(DataRate),names(DataRate)))
      DataPvalue[lower.tri(DataPvalue)] <- proxy
      
      if (input$batchApplyCorrection) {
        selectedMethod <- input$batchCorrectionMethod
        DataPvalue[] <- p.adjust(DataPvalue, method=selectedMethod)
      }
      
      DataPvalue[upper.tri(DataPvalue)] <- t(DataPvalue)[upper.tri(DataPvalue)] # do it after correction to avoid doubling; must use transpose to correctly fill second triangle
      
      output$BatchPvalueResults <- reactable::renderReactable({
        DataPvalueR <- matrix(NA,nrow=table.length,ncol=table.length,dimnames = list(names(DataRate),names(DataRate)))
        m <- mapply(function(x,y) {colourPvalue(DataPvalue[x,y])}, x = xarg, y = yarg)
        DataPvalueR[lower.tri(DataPvalueR)] <- m
        DataPvalueR[upper.tri(DataPvalueR)] <- t(DataPvalueR)[upper.tri(DataPvalueR)]
        rownames(DataPvalueR) <- sapply(rownames(DataPvalue), function(x) paste("<b>", x, "</b>", sep = ""))
        defaultTable(DataPvalueR, "Strain")
      })
      
      rv$PvalueOutput <- data.frame("Strain" = colnames(DataPvalue), as.data.frame(DataPvalue, check.names = FALSE))
      shinyjs::show("BatchPvalue")
      
    } else {
      shinyjs::hide("BatchPvalue")
      rv$PvalueOutput <- NULL
    }
    
    shinyjs::show("onOutputGenerated")
    
    shinyjs::enable("calculate4")
    shinyjs::runjs("$('#calculate4').html(\"<i class = 'fas fa-calculator'></i> Calculate\")")
    gc()
  })
  
  #### BatchCalc file export ####
  output$dlBatch <- downloadHandler(
    filename = function() {
      paste("Report-", Sys.Date(), ".xlsx", sep = "")
    },
    content = function(file) {
      FileList <- list("Rates" = rv$RateOutput)
      if (is.null(rv$PvalueOutput)==FALSE) {
        FileList$`P values` <- rv$PvalueOutput
      }
      writexl::write_xlsx(FileList, path = file)
    }
  )
  
  session$onSessionEnded(function() {
    stopApp()
  })
  
}