#
# This is the server logic of a Shiny web application. You can run the 
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#


# Set the size limit for file uploads. If unset, the maximum request size defaults to 5MB.
options(shiny.maxRequestSize = 30 * 1024^2)

pkgTest <- function(x) {
  if (!require(x, character.only = TRUE)) {
    install.packages(x, dep = TRUE, repos = "http://cran.rstudio.com/")
    if (!require(x, character.only = TRUE))
      stop("Package not found")
  }
}
pkgTest("shiny")


# Maximum capacity for sequencing (192 sequences per week)
maxSeqCapacity <- 192
# Maximum capacity for sequencing (6000 sequences per year)
maxAnnualCapacity <- 6000




ageGroups <- c(-1, 1, 6, 14, 65, 100)
selectedWeek <- NA


smallFont <- list(
  # family = "Courier New, monospace",
  size = 10,
  # color = "#7f7f7f"
  color = "black"
)




lastDayOfWeek <- function(weekNum, year) {
  # Option 1. UK week
  # first.day <- as.numeric(format(as.Date(paste(year, "-01-01", sep = "")), "%w"))
  # return(as.Date(paste(year, "-01-01", sep = "")) + weekNum * 7 - first.day)
  
  
  # as.Date(as.POSIXlt("2016 1 1", format = "%Y %U %u")) + 1
  
  # Option 2. ISO week
  return(as.Date(as.POSIXlt(paste(year, weekNum - 1, "6"), format = "%Y %U %u")) + 1)
  # return(strptime(paste0(year, "Sunday", weekNum), "%Y%A%U"))
}


getRegionData <- function(sData, region) {
  region <- as.character(region)
  sData$region <- as.character(sData$region)
  
  # sData$region <- replace(sData$region, which(Encoding(sData$region) == "UTF-8"), 
  #                         iconv(sData$region[which(Encoding(sData$region) == "UTF-8")], 
  #                               from = "UTF-8", to = "ASCII//TRANSLIT"))
  data <- sData
  
  # if (region != "All regions")
  #   data <- sData[which(substr(iconv(sData$region, from = "windows-1252", to = "ASCII//TRANSLIT"), 1, 4) == 
  #                         substr(iconv(region, from = "UTF-8", to = "ASCII//TRANSLIT"), 1, 4)), ]
  
  if (region != "All regions")
    data <- sData[which(substr(sData$region, 1, 4) == substr(region, 1, 4)), ]
  
  return(data)
}

getHistForRegion <- function(region, year, sData, datesType) {
  data <- getRegionData(sData, region)
  
  if (datesType == "Date of reception")
    dates <- data$Date.de.reception
  else
    dates <- data$Date.isolement
  
  # dates <- data$Date.de.reception
  
  dates <- dates[which(dates >= as.Date(paste(year, "-01-01", sep = "")) & 
                         dates < as.Date(paste(as.integer(year) + 1, "-01-01", sep = "")))]
  
  # Option 1. UK week
  # weeks <- as.integer(format(dates,
  #                            format = "%W" # UK week
  #                            # format = "%G-W%V" # ISO week
  #                            ))
  
  # Option 2. ISO week
  weeks <- data.table::isoweek(dates)
  
  # h <- hist(weeks, breaks = 53, plot = F)
  h <- hist(weeks, breaks = seq(-1, 53, 1), plot = F)
  # h <- hist(weeks, breaks = seq(min(weeks), max(weeks), 1), plot = FALSE)
  
  # print(range(weeks))
  
  return(h)
}

# Returns monthly histogram for provided data
getMonthlyHist <- function(data, datesType) {
  if (datesType == "Date of reception")
    dates <- data$Date.de.reception
  else
    dates <- data$Date.isolement
  
  months <- data.table::month(dates)
  h <- hist(months, breaks = seq(0, 12, 1), plot = F)
  
  return(h)
}

getDateBreaks <- function(dates, year) {
  if (length(dates) > 0) {
    # Convert each date to starting days of respective weeks
    dateBreaks <- as.Date(cut(dates, breaks = "week", start.on.monday = TRUE))
    return(dateBreaks)
  } else {
    print("No dateBreaks!")
    return(c())
  }
}

findAlerts <- function(dateBreaks, year, method) {
  require(surveillance)
  
  firstWeek <- as.Date(cut(as.Date(paste(year, "-01-01", sep = "")), breaks = "week", start.on.monday = TRUE))
  lastWeek <- as.Date(cut(as.Date(paste(year, "-12-31", sep = "")), breaks = "week", start.on.monday = TRUE))
  
  print(paste0("First week: ", firstWeek))
  print(paste0("Last week: ", lastWeek))
  
  # Exclude data from the following years
  dateBreaks <- dateBreaks[which(dateBreaks <= lastWeek)]
  
  
  if (length(dateBreaks) == 0) {
    print("length(dateBreaks) == 0")
    return(NULL)
  }
  
  # salmonellaDisProg <- linelist2sts(sData, dateCol = "Date.de.reception", aggregate.by = "1 week")
  
  salmonellaDisProg <- linelist2sts(data.frame(counts = dateBreaks), dateCol = "counts", aggregate.by = "1 week")
  
  # as.Date(salmonellaDisProg@epoch, origin = "1970-01-01")
  
  #Plot the result
  # plot(salmonellaDisProg, main = "Salmonella", 
  #      xaxis.tickFreq = list("%d" = atChange, "%m" = atChange),
  #      xaxis.labelFreq = list("%d" = at2ndChange),
  #      xaxis.labelFormat = "%d %b",
  #      xlab = "", las = 2, cex.axis = 0.5)
  
  
  
  Num_serotype <- sts2disProg(salmonellaDisProg)
  
  # print(firstWeek)
  # print(range(dateBreaks))
  # print(sort(unique(dateBreaks)))
  
  # start <- which(sort(unique(dateBreaks)) == firstWeek)
  start <- which(sort(unique(dateBreaks)) >= firstWeek)[1] # Closest dateBreak to firstWeek, which is after firstWeek
  
  # end <- which(sort(unique(dateBreaks)) == lastWeek)
  end <- which(sort(unique(dateBreaks)) == max(dateBreaks))
  # end <- length(Num_serotype$observed) # Should be the same after we excluded data from the following years
  
  # print(paste0("start: ", start))
  # print(paste0("end: ", end))
  # print(length(start:end))
  # print(lastWeek)
  
  
  if (is.na(start) | is.na(end)) {
    print("start or end is NA")
    return(NULL)
  }
  
  cntrl  <- list(range = start:end,
                 m = 1,
                 w = 3, # windows size, i.e. number of weeks to include before and after the current week
                 b = 4, # how many years back in time to include when forming the base counts.
                 alpha = 0.01, # An approximate (two-sided) (1 − α) prediction interval is calculated.
                 reweight = TRUE # Boolean specifying whether to perform reweight step
                 # , name = serotypeName#as.vector(ListeSerocomp[i,2])
  )
  
  # sts.cdsc <- algo.farrington(Num_serotype, control = cntrl)
  # sts.cdc <- algo.cdc(Num_serotype, control = cntrl)
  # sts.rki3 <- algo.rki3(Num_serotype, control = cntrl)
  
  status <- tryCatch({
    if (method == "farrington")
      sts <- algo.farrington(Num_serotype, control = cntrl)
    else if (method == "bayes")
      sts <- algo.bayes(Num_serotype, control = cntrl)
    else if (method == "cdc")
      sts <- algo.cdc(Num_serotype, control = cntrl)
    else if (method == "rki1")
      sts <- algo.rki1(Num_serotype, control = cntrl)
    else if (method == "rki2")
      sts <- algo.rki2(Num_serotype, control = cntrl)
    else if (method == "rki3")
      sts <- algo.rki3(Num_serotype, control = cntrl)
    0
  }, warning = function(w) {
    warning(w)
    NULL
  }, error = function(e) {
    print(e)
    NULL
    # }, finally = {
    #   cleanup-code
  })
  
  if (is.null(status))
    # return(c())
    return(NULL)
  
  # plot(sts, startyear = year)
  
  start <- start - 1
  
  # print(sts)
  # print(sort(unique(dateBreaks))[which(sts$alarm == 1) + start])
  
  return(sort(unique(dateBreaks))[which(sts$alarm == 1) + start])
}


# Update old regions with new ones in vector
useNewRegions <- function(vec) {
  vec <- replace(vec, which(substr(vec, 1, 2) == "Rh"), "Auvergne-Rhone-Alpes")
  vec <- replace(vec, which(substr(vec, 1, 4) == "La R"), "La Reunion")
  vec <- replace(vec, which(substr(vec, 1, 9) == "Nouvelle "), "Nouvelle Caledonie")
  vec <- replace(vec, which(substr(vec, 1, 4) == "Poly"), "Polynesie Francaise")
  vec <- replace(vec, which(substr(vec, 1, 4) == "Prov"), "Provence-Alpes-Cote d'Azur")
  vec <- replace(vec, which(substr(vec, 1, 4) == "Pays"), "Pays de la Loire")
  
  vec <- replace(vec, which(substr(vec, 1, 4) == "Bour"), "Bourgogne-Franche-Comte")
  vec <- replace(vec, which(substr(vec, 1, 4) == "Fran"), "Bourgogne-Franche-Comte")
  
  vec <- replace(vec, which(substr(vec, 1, 4) == "Aqui"), "Nouvelle-Aquitaine")
  vec <- replace(vec, which(substr(vec, 1, 4) == "Limo"), "Nouvelle-Aquitaine")
  vec <- replace(vec, which(substr(vec, 1, 4) == "Poit"), "Nouvelle-Aquitaine")
  
  vec <- replace(vec, which(substr(vec, 1, 4) == "Haut"), "Normandie")
  vec <- replace(vec, which(substr(vec, 1, 4) == "Bass"), "Normandie")
  
  vec <- replace(vec, which(substr(vec, 1, 4) == "Midi"), "Occitanie")
  vec <- replace(vec, which(substr(vec, 1, 4) == "Lang"), "Occitanie")
  
  vec <- replace(vec, which(substr(vec, 1, 4) == "Alsa"), "Grand Est")
  vec <- replace(vec, which(substr(vec, 1, 4) == "Cham"), "Grand Est")
  vec <- replace(vec, which(substr(vec, 1, 4) == "Lorr"), "Grand Est")
  
  vec <- replace(vec, which(substr(vec, 1, 4) == "Nord"), "Hauts-de-France")
  vec <- replace(vec, which(substr(vec, 1, 4) == "Pica"), "Hauts-de-France")
  
  vec <- replace(vec, which(substr(vec, 1, 4) == "Auve"), "Auvergne-Rhone-Alpes")
  
  return(vec)
}





##########################################################################

library(shiny)
library(V8) # Necessary for shinyjs

# Define server logic
shinyServer(function(input, output, session) {
  
  shinyjs::disable("button")
  
  loadedData <- reactive({
    # input$file1 will be NULL initially. After the user selects
    # and uploads a file, it will be a data frame with 'name',
    # 'size', 'type', and 'datapath' columns. The 'datapath'
    # column will contain the local filenames where the data can be found.
    
    inFile <- input$file
    
    if (is.null(inFile))
      return(NULL)
    
    # print(inFile$datapath)
    
    
    salmData <- read.csv(inFile$datapath, sep = ";", encoding = "macintosh", stringsAsFactors = FALSE)
    
    # Try to read with comma as separator if semi-colon does not work
    if (ncol(salmData) == 1)
      salmData <- read.csv(inFile$datapath, sep = ",", encoding = "macintosh", stringsAsFactors = FALSE)
    
    # Use Serotype_FINAL2 if it is in the file
    if ("Serotype_FINAL2" %in% colnames(salmData))
      salmData$Serotype_FINAL <- salmData$Serotype_FINAL2
    
    
    # fileEncoding = "macintosh"
    
    
    # print(inFile$datapath)
    
    # library(readxl)
    # salmData <- read_xls(inFile$datapath, sheet = 1)
    # colnames(salmData) <- gsub('-', '.', colnames(salmData))
    # colnames(salmData) <- gsub(' ', '.', colnames(salmData))
    
    print(paste0("Loaded data: ", nrow(salmData), " rows"))
    
    salmData$Date.isolement <- as.Date(salmData$Date.isolement, "%d/%m/%Y")
    salmData$Date.de.reception <- as.Date(salmData$Date.de.reception, "%d/%m/%Y")
    
    # Update CT types: "CT_1" -> "1"
    str <- as.character(salmData$CRISPOL.type)
    salmData$CRISPOL.type <- ifelse(startsWith(str, "CT_"), substr(str, 4, length(str)), str)
    
    sData <- salmData
    sData$region <- as.character(sData$region)
    
    sData$region <- enc2native(sData$region)
    # sData$region <- replace(sData$region, which(Encoding(sData$region) == "UTF-8"), 
    #                         iconv(sData$region[which(Encoding(sData$region) == "UTF-8")], 
    #                               from = "UTF-8", to = "ASCII//TRANSLIT"))
    sData$region <- useNewRegions(as.vector(sData$region))
    
    
    sData$Serotype_FINAL <- enc2native(sData$Serotype_FINAL)
    
    # Remove accents
    sData$Voyage <- enc2native(sData$Voyage)
    sData$Epidemie <- enc2native(sData$Epidemie)
    # sData$Voyage <- iconv(sData$Voyage, from = "UTF-8", to = "ASCII//TRANSLIT")
    # sData$Epidemie <- iconv(sData$Epidemie, from = "UTF-8", to = "ASCII//TRANSLIT")
    
    sData$Voyage <- replace(sData$Voyage, which(startsWith(sData$Voyage, "Non pr")), "Not specified")
    sData$Epidemie <- replace(sData$Epidemie, which(startsWith(sData$Epidemie, "Isol")), "Isolated")
    
    updateSelectInput(session, inputId = "region", choices = c("All regions", sort(unique(sData$region))))
    updateSelectInput(session, inputId = "serotype", choices = c("All serotypes", sort(unique(sData$Serotype_FINAL))))
    
    updateSelectInput(session, inputId = "ageGroup", choices = c("All age groups", sapply(1:(length(ageGroups) - 1), function(x) sprintf("%i to %i", ageGroups[x] + 1, ageGroups[x + 1]))))
    
    updateSelectInput(session, inputId = "department", choices = c("All departments", sort(unique(sData$Departement))))
    
    # TODO: update subtypes menu
    # selSerotype <- getSelectedSerotype()
    # if (selSerotype != "All serotypes") {
    #   tmp <- sData[which(sData$Serotype_FINAL == selSerotype), ]
    #   
    #   if (selSerotype == "Typhimurium" | selSerotype == "4,12:i:- (monophasique)" | selSerotype == "4,5,12:i:- (monophasique)"
    #       # | s$x == "Enteritidis"
    #   )
    #     updateSelectInput(session, inputId = "subtype", choices = c("All subtypes", sort(unique(tmp$CRISPOL.type))))
    #   else
    #     updateSelectInput(session, inputId = "subtype", choices = c("All subtypes", sort(unique(tmp$profil_MLST))))
    # }
    
    updateSelectInput(session, inputId = "subtype", choices = c("All subtypes", sort(unique(sData$profil_MLST))))
    
    # require(lubridate)
    # years <- year(sData$Date.de.reception)
    # print(sData$Date.de.reception)
    
    years <- substring(as.character(sData$Date.de.reception), 1, 4)
    
    # print(unique(years))
    
    years <- sort(unique(years))
    
    updateSelectInput(session, inputId = "year", choices = years, selected = max(years))
    
    updateSelectInput(session, inputId = "week", choices = c("All weeks", 1:53), selected = 2)
    
    selectedWeek <<- lastDayOfWeek(1, max(years)) + 1
    
    return(sData)
  })
  
  # Reactive conductor for selected region
  getSelectedRegion <- reactive({
    s <- event_data("plotly_click", source = "regPlot")
    if (length(s)) {
      if (input$region != s$x) {
        updateSelectInput(session, inputId = "region", selected = s$x)
        
        # data <- loadedData()
        # tmp <- data[which(data$region == s$x), ]
        # 
        # print(paste0("Available deps: ", sort(unique(tmp$Departement))))
        # updateSelectInput(session, inputId = "department", choices = c("All departments", sort(unique(tmp$Departement))))
        # updateSelectInput(session, inputId = "department", selected = "All departments")
      } else {
        updateSelectInput(session, inputId = "region", selected = "All regions")
        
        # updateSelectInput(session, inputId = "department", choices = c("All departments"))
        # updateSelectInput(session, inputId = "department", selected = "All departments")
      }
      
      
      # # Update department menu: TODO: cancel double-update
      # selRegion <- input$region
      # selDep <- getSelectedDepartment()
      # if (selRegion != "All regions") {
      #   tmp <- data[which(data$region == selRegion), ]
      #   print(paste0("!!! ", selDep))
      #   
      #   if (selDep != "All departments") {
      #     updateSelectInput(session, inputId = "department", choices = c("All departments", sort(unique(tmp$Departement))))
      #     # updateSelectInput(session, inputId = "department", selected = selDep)
      #   }
      # }
      
      js$resetRegClick()
    }
    
    return(input$region)
  })
  
  # Reactive conductor for selected serotype
  getSelectedSerotype <- reactive({
    click <- event_data("plotly_click", source = "serPlot")
    if (length(click)) {
      if (input$serotype != click$x) {
        updateSelectInput(session, inputId = "serotype", selected = click$x)
      } else {
        updateSelectInput(session, inputId = "serotype", selected = "All serotypes")
      }
      
      js$resetSerClick()
    }
    
    return(input$serotype)
  })
  
  getSelectedWeek <- function() {
    print(paste0("input week: ", input$week))
    
    # s <- event_data("plotly_click", source = "mainPlot")
    # if (length(s)) {
    #   selectedWeek <- as.Date(s$x, origin = "1970-01-01")
    #   weekNum <- data.table::isoweek(selectedWeek)
    #   updateSelectInput(session, inputId = "week", selected = weekNum)
    #   
    #   # js$resetWeekClick()
    # }
    # 
    # if (input$week == "All weeks") {
    #   js$resetWeekClick()
    #   return("All weeks")
    # }
    # 
    # # if (input$week != ...)
    # 
    # # selectedWeek <- as.Date(input$week, origin = "1970-01-01")
    # return(selectedWeek)
    
    
    s <- event_data("plotly_click", source = "mainPlot")
    if (length(s)) {
      print(paste0("s$x: ", s$x))
      if (input$week != data.table::isoweek(as.Date(s$x, origin = "1970-01-01"))) { # Uncomment to enable selection of all weeks by double-click
        selectedWeek <- as.Date(s$x, origin = "1970-01-01")
        weekNum <- data.table::isoweek(selectedWeek)
        # weekNum <- lubridate::isoweek(selectedWeek)
        updateSelectInput(session, inputId = "week", selected = weekNum)
      } else {
        updateSelectInput(session, inputId = "week", selected = "All weeks")
        selectedWeek <- "All weeks"
      }
      
      js$resetWeekClick()
    }
    
    if (input$week == "All weeks")
      return("All weeks")
    
    if (input$week == "")
      return(NULL)
    
    selectedWeek <- ISOweek::ISOweek2date(paste0(input$year, "-W", sprintf("%02d", as.integer(input$week)), "-1"))
    
    return(selectedWeek)
  }
  
  
  # Reactive conductor for selected serotype
  getSelectedAgeGroup <- reactive({
    click <- event_data("plotly_click", source = "agePlot")
    if (length(click)) {
      # print(paste0("clicked: ", click$x))
      up <- as.integer(strsplit(click$x, " ")[[1]][2])
      i <- which(ageGroups == up)
      selAgeGroup <- sprintf("%i to %i", ageGroups[i - 1] + 1, ageGroups[i])
      
      print(paste0("Selected age group: ", selAgeGroup))
      
      if (input$ageGroup != selAgeGroup)
        updateSelectInput(session, inputId = "ageGroup", selected = selAgeGroup)
      else
        updateSelectInput(session, inputId = "ageGroup", selected = "All age groups")
      
      js$resetAgeClick()
    }
    
    return(input$ageGroup)
  })
  
  
  # Reactive conductor for selected department
  getSelectedDepartment <- reactive({
    click <- event_data("plotly_click", source = "depPlot")
    if (length(click)) {
      print(paste0("Selected department: ", click$x))
      print(paste0("Current department: ", input$department))
      
      if (input$department != click$x) {
        updateSelectInput(session, inputId = "department", selected = click$x)
        print(paste0("New department: ", input$department))
      } else {
        updateSelectInput(session, inputId = "department", selected = "All departments")
        print(paste0("New department: ", input$department))
      }
      
      js$resetDepClick()
      
      
      
      if (input$region == "All regions") {
        data <- loadedData()
        tmp <- data[which(data$Departement == click$x), ]
        print(paste0("Change to region: ", input$department, "", click$x))
        print(paste0("Change to region: ", tmp$region[1]))
        
        updateSelectInput(session, inputId = "region", selected = tmp$region[1])
        
        js$resetRegClick()
      }
    }
    
    return(input$department)
  })
  
  getSelectedSubtype <- reactive({
    # click <- event_data("plotly_click", source = "agePlot")
    # if (length(click)) {
    #   # print(paste0("clicked: ", click$x))
    #   up <- as.integer(strsplit(click$x, " ")[[1]][2])
    #   i <- which(ageGroups == up)
    #   selAgeGroup <- sprintf("%i to %i", ageGroups[i - 1], ageGroups[i])
    #   
    #   if (input$ageGroup != selAgeGroup)
    #     updateSelectInput(session, inputId = "ageGroup", selected = click$x)
    #   else
    #     updateSelectInput(session, inputId = "ageGroup", selected = "All age groups")
    # }
    
    return(input$subtype)
  })
  
  # Returns data for selected week, region and serotype
  subsetData <- function(data, selectedRegion = getSelectedRegion(), 
                         selectedWeek = getSelectedWeek(), 
                         selectedSerotype = getSelectedSerotype(), 
                         selectedSubtype = getSelectedSubtype(),
                         selectedDepartment = getSelectedDepartment(),
                         selectedAgeGroup = getSelectedAgeGroup(),
                         selectedYear = input$year) {
    regData <- getRegionData(data, selectedRegion)
    
    if (!is.null(selectedWeek)) {
      if (as.character(selectedWeek) == "All weeks") {
        firstDay <- as.Date(paste(selectedYear, "-01-01", sep = ""))
        lastDay <- as.Date(paste(selectedYear, "-12-31", sep = ""))
          
        if (input$datesType == "Date of reception")
          regData <- regData[which(regData$Date.de.reception >= firstDay &
                                     regData$Date.de.reception <= lastDay), ]
        else
          regData <- regData[which(regData$Date.isolement >= firstDay &
                                     regData$Date.isolement <= lastDay), ]
      } else {
        if (input$datesType == "Date of reception")
          regData <- regData[which(regData$Date.de.reception - selectedWeek >= 0 &
                                   regData$Date.de.reception - selectedWeek < 7), ]
        else
          regData <- regData[which(regData$Date.isolement - selectedWeek >= 0 &
                                   regData$Date.isolement - selectedWeek < 7), ]
      }
    }
    
    # TODO: Error in if: missing value where TRUE/FALSE needed (29 Jan 2018, Serotype alerts)
    if (selectedSerotype != "All serotypes")
      regData <- regData[which(regData$Serotype_FINAL == selectedSerotype), ]
    
    if (selectedAgeGroup != "All age groups") {
      print(selectedAgeGroup)
      age <- as.integer(strsplit(selectedAgeGroup, " ")[[1]][1]) - 1
      age2 <- ageGroups[which(ageGroups == age) + 1]
      
      ages <- regData$Age
      ages[which(ages == "< 1")] <- rep(0, length(which(ages == "< 1")))
      ages <- as.integer(ages)
      
      print(ages)
      
      print(paste0("Age between ", age, " and ", age2))
      
      print(which(ages > age & ages <= age2))
      regData <- regData[which(ages > age & ages <= age2), ]
      
      
      print(paste0(nrow(regData), " rows"))
    }
    
    if (selectedDepartment != "All departments") {
      regData <- regData[which(regData$Departement == selectedDepartment), ]
    }
    
    if (selectedSubtype != "All subtypes") {
      if (selectedSerotype == "Typhimurium" 
          | selectedSerotype == "4,12:i:- (monophasique)" 
          | selectedSerotype == "4,5,12:i:- (monophasique)"
          # | selectedSerotype == "Enteritidis"
      )
        regData <- regData[which(regData$CRISPOL.type == selectedSubtype), ]
      else
        regData <- regData[which(regData$profil_MLST == selectedSubtype), ]
    }
    
    return(regData)
  }
  
  
  # Main (temporal) distribution
  output$temporalPlotly <- renderPlotly({
    data <- loadedData()
    
    # Update subtype menu
    selSerotype <- input$serotype
    if (selSerotype != "All serotypes") {
      tmp <- data[which(data$Serotype_FINAL == selSerotype), ]
      
      if (selSerotype == "Typhimurium" | selSerotype == "4,12:i:- (monophasique)" | selSerotype == "4,5,12:i:- (monophasique)"
          # | s$x == "Enteritidis"
      ) {
        storeSubtype <- input$subtype
        updateSelectInput(session, inputId = "subtype", choices = c("All subtypes", sort(unique(tmp$CRISPOL.type))))
        updateSelectInput(session, inputId = "subtype", selected = storeSubtype)
        print(paste0("new subtype: ", input$subtype))
        if (input$subtype == "")
          updateSelectInput(session, inputId = "subtype", selected = "All subtypes")
      } else {
        storeSubtype <- input$subtype
        updateSelectInput(session, inputId = "subtype", choices = c("All subtypes", sort(unique(tmp$profil_MLST))))
        updateSelectInput(session, inputId = "subtype", selected = storeSubtype)
        print(paste0("new subtype: ", input$subtype))
        if (input$subtype == "")
          updateSelectInput(session, inputId = "subtype", selected = "All subtypes")
      }
    } else {
      updateSelectInput(session, inputId = "subtype", selected = "All subtypes")
    }
    
    updateDepMenu <- reactive({
      if (input$region != "All regions") {
        tmp <- data[which(data$region == input$region), ]
        
        print(paste0("Available deps: ", sort(unique(tmp$Departement))))
        updateSelectInput(session, inputId = "department", choices = c("All departments", sort(unique(tmp$Departement))))
        # updateSelectInput(session, inputId = "department", selected = "All departments")
      } else {
        updateSelectInput(session, inputId = "department", choices = c("All departments"))
        # updateSelectInput(session, inputId = "department", selected = "All departments")
      }
    })
    # updateDepMenu()
    
    if (is.null(data)) {
      plotly_empty(type = "scatter", mode = "markers") %>%
        layout(title = "<- Please upload data to start", titlefont = list(
          # family = 'Courier New, monospace',
          size = 18,
          color = "red"#'#7f7f7f'
        )) 
    } else {
      data <- subsetData(data, selectedWeek = NULL)
      year <- input$year
      
      # print(nrow(data))
      
      if (input$datesType == "Date of reception")
        dates <- data$Date.de.reception
      else
        dates <- data$Date.isolement
      
      dateBreaks <- getDateBreaks(dates, year)
      
      # print(paste0("YEAR: ", year))
      # print(year == "")
      # print(is.null(year))
      # print(is.na(year))
      # print(dateBreaks)
      if (year == "")
        return(NULL)
      
      firstWeek <- as.Date(cut(as.Date(paste(year, "-01-01", sep = "")), breaks = "week", start.on.monday = TRUE))
      lastWeek <- as.Date(cut(as.Date(paste(year, "-12-31", sep = "")), breaks = "week", start.on.monday = TRUE))
      
      # print(sort(unique(dateBreaks)))
      
      # Exclude data from other years
      dateBrks <- dateBreaks[which(dateBreaks <= lastWeek & firstWeek <= dateBreaks)]
      
      print(sort(unique(dateBrks)))
      
      
      if (length(dateBrks) != 0) {
        h <- hist(dateBrks, breaks = "week", start.on.monday = FALSE, plot = F)
        
        # Highlight selected week
        colors <- rep('rgba(0,0,0, 0)', length(h$breaks) - 1) # Transparent
        if (as.character(getSelectedWeek()) != "All weeks")
          colors[which(as.Date(h$breaks[-length(h$breaks)] + 1, origin = "1970-01-01") == getSelectedWeek())] <- "red"
        else
          colors <- rep("red", length(h$breaks) - 1) # Select all 
        
        print(paste0("h: ", h))
        
        p <- plot_ly(x = as.Date(h$breaks[-length(h$breaks)] + 1, origin = "1970-01-01"), y = h$counts, type = 'bar', #line = list(color = '#45171D'),
                     hoverinfo = "text",
                     marker = list(line = list(color = colors, width = 1.5)),
                     text = paste(h$counts, ' samples', " on week ",
                                  data.table::isoweek(as.Date(h$breaks[-length(h$breaks)] + 1, origin = "1970-01-01")), " (",
                                  format(as.Date(h$breaks[-length(h$breaks)] + 1, origin = "1970-01-01"), format = "%d-%b"), " to ",
                                  format(as.Date(h$breaks[-length(h$breaks)] + 7, origin = "1970-01-01"), format = "%d-%b"), ")",
                                  sep = ""),
                     name = year,
                     source = "mainPlot")
        
        
        # Visualise alerts
        pal <- c("red", "orange", "purple", "pink", "salmon")
        pal2 <- c("rgba(0, 0, 55, .3)", 
                  "rgba(0, 0, 55, .3)", 
                  "rgba(0, 0, 55, .3)", 
                  "rgba(0, 0, 55, .3)", 
                  "rgba(0, 0, 55, .3)")
        
        i <- 0
        for (method in c("bayes", "rki1", "farrington")) {
          i <- i + 1
          alerts <- findAlerts(dateBreaks, year, method)

          if (!is.null(alerts)) {
            numCasesAlerts <- h$counts[which(h$breaks %in% as.integer(alerts - 1))]

            print(paste0(method, " alerts: ", paste(alerts, collapse = ", ")))
            # print(numCasesAlerts)
            print(as.Date(h$breaks + 1, origin = "1970-01-01"))

            p <- add_trace(p, x = alerts,
                           y = numCasesAlerts,
                           marker = list(color = pal[i], size = 14 - i * 3),
                           inherit = FALSE,
                           type = 'scatter', mode = 'markers',
                           # hoverinfo = "none",
                           name = paste0("Alerts (", method, ")") #yaxis = 'y2',
                           )

            no_alerts <- h$breaks[which(!(h$breaks %in% as.integer(alerts - 1)))]
            no_alerts <- as.Date(no_alerts + 1, origin = "1970-01-01")
            no_numCasesAlerts <- h$counts[which(!(h$breaks %in% as.integer(alerts - 1)))]

            p <- add_trace(p, x = no_alerts,
                           y = no_numCasesAlerts,
                           marker = list(color = pal2[i], size = 14 - i * 3),
                           inherit = FALSE,
                           type = 'scatter', mode = 'markers',
                           # hoverinfo = "none",
                           name = paste0("No alerts (", method, ")") #yaxis = 'y2',
                           )
          }
        }
      } else { # If there are no samples for selected year
        p <- plot_ly(x = c(), y = c(), type = 'bar',
                     hoverinfo = "none",
                     marker = list(line = list(color = colors, width = 1.5)),
                     name = year,
                     source = "mainPlot")
      }
      
      p <- layout(p, title = paste0("Weekly distribution of ", getSelectedSerotype()," samples" , 
                                    ifelse(getSelectedSubtype() == "All subtypes", "", paste0(" (subtype ", getSelectedSubtype(), ")")),
                                    "<br>",
                                    "in ", getSelectedRegion(), 
                                    ifelse(getSelectedDepartment() == "All departments", "", paste0(" (dep ", getSelectedDepartment(), ")")),
                                    ifelse(getSelectedAgeGroup() == "All age groups", "", paste0(" (age ", getSelectedAgeGroup(), ")")),
                                    " in ", year), titlefont = smallFont, 
                  showlegend = TRUE,
                  xaxis = list(title = "Week", showgrid = FALSE, zeroline = FALSE, showticklabels = TRUE, tickfont = smallFont, titlefont = smallFont),
                  yaxis = list(title = "Number of samples", showgrid = TRUE, zeroline = FALSE, showticklabels = TRUE, tickfont = smallFont, titlefont = smallFont))
      
      return(p)
    }
  })
  
  
  # Regional distribution
  output$regionalPlotly <- renderPlotly({
    selectedWeek <- getSelectedWeek()
    
    data <- loadedData()
    if (is.null(data))
      return(plotly_empty(type = "scatter", mode = "markers"))
    
    regData <- subsetData(data, selectedRegion = "All regions", selectedDepartment = "All departments")
    
    sel <- which(regData$region == "")
    regData$region[sel] <- rep("Unknown", length(sel))
    
    # df <- as.data.frame(table(droplevels(regData$region))) # Was used for old regions
    df <- as.data.frame(table(regData$region))
    
    if (nrow(df) == 0) {
      # Create empty data frame with the same column names
      df <- data.frame(Var1 = c(0), Freq = c(0))
      df <- df[-1, ]
      colors <- c()
    } else {
      # df$Var1 <- iconv(df$Var1, to = "UTF-8")
      
      otherRegions <- c("Corse", "Guadeloupe", 
                        "Guyane","La Reunion",  "Martinique", 
                        "Mayotte", "Monaco", 
                        "Nouvelle Caledonie", "Polynesie Francaise", 
                        "Saint Pierre et Miquelon", 
                        "Etranger",
                        "Unknown")
      
      cols <- rep("Mainland France", length(df$Var1))
      sel <- which(df$Var1 %in% otherRegions)
      cols[sel] <- rep("Other regions", length(sel))
      cols[which(df$Var1 == "Unknown")] <- "Unknown"
      
      
      selectedRegion <- getSelectedRegion()
      
      colors <- rep("lightblue", length(df$Var1))
      sel <- which(df$Var1 %in% otherRegions)
      colors[sel] <- rep("cyan", length(sel))
      colors[which(df$Var1 == "Unknown")] <- "lightgray"
      colors[which(df$Var1 == selectedRegion)] <- "red"
    }
    
    p <- plot_ly(df, 
                 #labels = ~Var1, values = ~Freq, type = 'pie',
                 x = ~Var1, y = ~Freq, type = "bar", # barmode = "stack", 
                 name = "Regions",
                 # color = cols,
                 marker = list(color = "lightblue", line = list(color = colors, width = 1.5)),
                 # inherit = FALSE, # deprecated
                 source = "regPlot") %>%
      layout(title = paste0("Regional distribution of ", getSelectedSerotype(), " samples",
                            ifelse(getSelectedSubtype() == "All subtypes", "", paste0(" (subtype ", getSelectedSubtype(), ")")),
                            ifelse(getSelectedAgeGroup() == "All age groups", "", paste0(" (age ", getSelectedAgeGroup(), ")")),
                            " <br> on week starting ", selectedWeek),
             showlegend = FALSE, titlefont = smallFont,
             xaxis = list(title = "", showgrid = FALSE, zeroline = FALSE, showticklabels = TRUE, tickfont = smallFont, titlefont = smallFont),
             yaxis = list(title = "Number of samples", showgrid = TRUE, zeroline = FALSE, showticklabels = TRUE, tickfont = smallFont, titlefont = smallFont))
    
    # Regional alerts
    if (input$detectReg) {
      for (selRegion in df$Var1) {
        
        selData <- subsetData(data, selectedRegion = selRegion, selectedWeek = NULL)
        
        if (input$datesType == "Date of reception")
          dates <- selData$Date.de.reception
        else
          dates <- selData$Date.isolement
        
        dateBreaks <- getDateBreaks(dates, input$year)
        
        # Visualise alerts
        pal <- c("red", "orange", "purple", "pink", "salmon")
        pal2 <- c("rgba(0, 0, 55, .3)", 
                  "rgba(0, 0, 55, .3)", 
                  "rgba(0, 0, 55, .3)", 
                  "rgba(0, 0, 55, .3)", 
                  "rgba(0, 0, 55, .3)")
        
        i <- 0
        for (method in c("bayes", "rki1", "farrington")) {
          i <- i + 1
          alerts <- findAlerts(dateBreaks, input$year, method)
          
          print(paste0("Ser alerts: ", alerts))
          if(!is.null(alerts)) {
            sel <- which(alerts == selectedWeek)
            
            if (length(sel) == 1) {
              p <- add_trace(p, x = df$Var1[which(df$Var1 == selRegion)],
                             y = df$Freq[which(df$Var1 == selRegion)],
                             marker = list(color = pal[i], size = 14 - i * 3),
                             inherit = FALSE,
                             type = 'scatter', mode = 'markers',
                             name = paste0("Alerts (", method, ")"), #yaxis = 'y2',
                             hoverinfo = "none"
              )
            } else {
              p <- add_trace(p, x = df$Var1[which(df$Var1 == selRegion)],
                             y = df$Freq[which(df$Var1 == selRegion)],
                             marker = list(color = pal2[i], size = 14 - i * 3),
                             inherit = FALSE,
                             type = 'scatter', mode = 'markers',
                             name = paste0("No alerts (", method, ")"), #yaxis = 'y2',
                             hoverinfo = "none"
              )
            }
          }
        }
      }
    }
    
    return(p)
  })
  
  
  
  
  # Serotype distribution
  output$serotypePlotly <- renderPlotly({
    selectedSerotype <- getSelectedSerotype()
    
    data <- loadedData()
    
    if (is.null(data))
      return(plotly_empty(type = "scatter", mode = "markers"))
    
    regData <- subsetData(data, selectedSerotype = "All serotypes", selectedSubtype = "All subtypes")
    
    # Use Serotype_FINAL
    df <- as.data.frame(table(regData$Serotype_FINAL))
    df$Var1 <- iconv(df$Var1, to = "UTF-8")
    
    colors <- rep("rgba(0,0,0,0)", length(df$Var1))
    sel <- which(df$Var1 == "Typhimurium" | df$Var1 == "4,12:i:- (monophasique)" |
                   df$Var1 == "4,5,12:i:- (monophasique)")
    colors[sel] <- rep("lightgreen", length(sel))
    colors[which(df$Var1 == "Enteritidis")] <- "orange"
    
    if (selectedSerotype != "All serotypes")
      colors[which(df$Var1 == selectedSerotype)] <- "red"
    
    p <- plot_ly(df, #labels = ~Var1, values = ~Freq, type = 'pie',
                 x = ~Var1, y = ~Freq, type = "bar", # barmode = "stack",
                 marker = list(color = "lightblue", line = list(color = colors, width = 1.5)),
                 source = "serPlot") %>%
      layout(title = paste0("Serotype distribution in ", getSelectedRegion(), 
                            # ifelse(getSelectedSubtype() == "All subtypes", "", paste0(" (subtype ", getSelectedSubtype(), ")")),
                            ifelse(getSelectedDepartment() == "All departments", "", paste0(" (dep ", getSelectedDepartment(), ")")),
                            ifelse(getSelectedAgeGroup() == "All age groups", "", paste0(" (age ", getSelectedAgeGroup(), ")")),
                            "<br> on week starting ", getSelectedWeek()),
             titlefont = smallFont,
             showlegend = FALSE,
             xaxis = list(title = "", showgrid = FALSE, zeroline = FALSE, showticklabels = TRUE, tickfont = smallFont, titlefont = smallFont),
             yaxis = list(title = "Number of samples", showgrid = TRUE, zeroline = FALSE, showticklabels = TRUE, tickfont = smallFont, titlefont = smallFont))
    
    
    if (input$detectSer) {
      print("Detecting alerts for serotypes:")
      print(df$Var1)
      for (selSerotype in df$Var1) {
        selData <- subsetData(data, selectedWeek = NULL, selectedSerotype = selSerotype)
        
        if (input$datesType == "Date of reception")
          dates <- selData$Date.de.reception
        else
          dates <- selData$Date.isolement
        
        dateBreaks <- getDateBreaks(dates, input$year)
        
        # Visualise alerts
        pal <- c("red", "orange", "purple", "pink", "salmon")
        pal2 <- c("rgba(0, 0, 55, .3)", 
                  "rgba(0, 0, 55, .3)", 
                  "rgba(0, 0, 55, .3)", 
                  "rgba(0, 0, 55, .3)", 
                  "rgba(0, 0, 55, .3)")
        
        i <- 0
        for (method in c("bayes", "rki1", "farrington")) {
          i <- i + 1
          
          print("Finding alerts:")
          print(input$year)
          print(dateBreaks)
          alerts <- findAlerts(dateBreaks, input$year, method)
          
          if (!is.null(alerts)) {
            sel <- which(alerts == selectedWeek)
            
            if (length(sel) == 1) {
              p <- add_trace(p, x = df$Var1[which(df$Var1 == selSerotype)],
                             y = df$Freq[which(df$Var1 == selSerotype)],
                             marker = list(color = pal[i], size = 14 - i * 3),
                             inherit = FALSE,
                             type = 'scatter', mode = 'markers',
                             name = paste0("Alerts (", method, ")"), #yaxis = 'y2',
                             hoverinfo = "none"
              )
            } else {
              p <- add_trace(p, x = df$Var1[which(df$Var1 == selSerotype)],
                             y = df$Freq[which(df$Var1 == selSerotype)],
                             marker = list(color = pal2[i], size = 14 - i * 3),
                             inherit = FALSE,
                             type = 'scatter', mode = 'markers',
                             name = paste0("No alerts (", method, ")"), #yaxis = 'y2',
                             hoverinfo = "none"
              )
            }
          }
        }
      }
    }
    
    return(p)
  })
  
  
  
  graphChoice <- ""
  strChoice <- ""
  
  
  
  getSubtypeColors <- function(subtypeNames, resolution) {
    data <- loadedData()
    
    MLST <- data$profil_MLST # as.character(regData$profil_MLST)
    MLST[which(MLST == "")] <- "NA" # Warning: invalid factor level, NA generated
    
    str <- as.character(MLST)
    MLST <- ifelse(startsWith(str, "ST"), str, paste0("ST", str))
    
    if (resolution == "ST") {
      df <- as.data.frame(table((MLST)), stringsAsFactors = FALSE)
      # df$Var1 <- iconv(df$Var1, to = "UTF-8")
      df$Var1[which(df$Var1 == "")] <- "Not seq."
      df$Var1[which(df$Var1 == "STNA")] <- "Not seq."
    }
    else if (resolution == "ST-CT") {
      CT <- data$CRISPOL.type
      CT[which(CT == "")] <- "NA"
      
      df <- as.data.frame(table(paste(MLST, CT, sep = "-CT")), stringsAsFactors = FALSE)
      # df$Var1 <- iconv(df$Var1, to = "UTF-8")
      # df$Var1[which(startsWith(as.character(df$Var1), "NA-"))] <- "Not seq." # ?
      df$Var1[which(as.character(df$Var1) == "STNA-CTNA")] <- "Not seq."
    }
    
    # require(RColorBrewer)
    # colors <- brewer.pal(length(df$Var1), "Spectral")
    
    colors <- rainbow(length(df$Var1))
    
    # require(colorRamps)
    # # colors <- primary.colors(length(df$Var1))
    # colors <- matlab.like(length(df$Var1))
    
    if (resolution == "ST-CT") {
      set.seed(100)
      colors <- colors[sample(1:length(colors), length(colors), replace = F)]
    }
    
    
    colors <- colors[match(subtypeNames, df$Var1)]
    
    
    colors[which(subtypeNames == "Not seq.")] <- "gray"
    
    return(colors)
  }
  
  # Regional/MLST/age distribution
  output$dynamicPlotly <- renderPlotly({
    data <- loadedData()
    
    #################################################################
    
    s <- event_data("plotly_hover", source = "serPlot")
    if (length(s)) {
      graphChoice <<- "ser"
      strChoice <- s$x
    }
    
    s <- event_data("plotly_hover", source = "mainPlot")
    if (length(s)) {
      graphChoice <<- "temp"
      strChoice <- s$x
    }
    
    s <- event_data("plotly_hover", source = "regPlot")
    if (length(s)) {
      graphChoice <<- "reg"
      strChoice <- s$x
    }
    
    #################################################################
    if (graphChoice == "ser") {
      if (strChoice == "")
        strChoice <- getSelectedSerotype()
      
      serotype <- strChoice
      
      regData <- subsetData(data, selectedSerotype = serotype, selectedSubtype = "All subtypes")
      
      MLST <- regData$profil_MLST # as.character(regData$profil_MLST)
      MLST[which(MLST == "")] <- "NA" # Warning: invalid factor level, NA generated
      
      str <- as.character(MLST)
      MLST <- ifelse(startsWith(str, "ST"), str, paste0("ST", str))
      # print(MLST)
      
      
      if (serotype == "Typhimurium" | serotype == "4,12:i:- (monophasique)" | serotype == "4,5,12:i:- (monophasique)"
          # | s$x == "Enteritidis"
      ) {
        CT <- regData$CRISPOL.type
        CT[which(CT == "")] <- "NA"
        
        df <- as.data.frame(table(paste(MLST, CT, sep = "-CT")), stringsAsFactors = FALSE)
        # df$Var1 <- iconv(df$Var1, to = "UTF-8")
        # df$Var1[which(startsWith(as.character(df$Var1), "NA-"))] <- "Not seq." # ?
        df$Var1[which(as.character(df$Var1) == "STNA-CTNA")] <- "Not seq."
        
        resolution <- "ST-CT"
      } else {
        df <- as.data.frame(table((MLST)), stringsAsFactors = FALSE)
        # df$Var1 <- iconv(df$Var1, to = "UTF-8")
        df$Var1[which(df$Var1 == "")] <- "Not seq."
        df$Var1[which(df$Var1 == "STNA")] <- "Not seq."
        
        resolution <- "ST"
      }
      
      colors <- getSubtypeColors(df$Var1, resolution)
      
      # # To check colors
      # colors2 <- brewer.pal(max(length(df$Var1), 3), "Spectral")
      # colors2[which(df$Var1 == "Not seq.")] <- "gray"
      # print("df$Var1: ")
      # print(df$Var1)
      # print(data.frame(x = colors, x2 = colors2))
      
      
      p <- plot_ly(df, labels = ~Var1, values = ~Freq, type = 'pie',
                   marker = list(colors = colors),
                   # textposition = 'outside',
                   textinfo = 'label+percent+value',
                   textfont = smallFont,
                   hole = 0.3,
                   showlegend = FALSE
                   #'font', 'title', 'titlefont', 'autosize', 'width', 'height', 'margin', 'paper_bgcolor', 'plot_bgcolor', 'separators', 'hidesources', 'smith', 'showlegend', 'dragmode', 'hovermode', 'xaxis', 'yaxis', 'scene', 'geo', 'legend', 'annotations', 'shapes', 'images', 'updatemenus', 'ternary', 'mapbox', 'radialaxis', 'angularaxis', 'direction', 'orientation', 'barmode', 'bargap', 'mapType'
      ) %>%
        layout(title = paste0("Breakdown of ", serotype, " by ", resolution, " <br> in ", getSelectedRegion(), 
                              ifelse(getSelectedDepartment() == "All departments", "", paste0(" (dep ", getSelectedDepartment(), ")")),
                              ifelse(getSelectedAgeGroup() == "All age groups", "", paste0(" (age ", getSelectedAgeGroup(), ")")),
                              " on week starting ", getSelectedWeek()),
               titlefont = smallFont,
               # legend = list(font = smallFont),
               xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
               yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE))
      
      return(p)
    }
    
    #################################################################
    if (graphChoice == "temp") {
      if (strChoice == "")
        strChoice <- getSelectedWeek()
      
      if (as.character(strChoice) == "All weeks") {
        selDate <- "All weeks" 
        regData <- subsetData(data, selectedWeek = "All weeks", selectedAgeGroup = "All age groups")
      } else {
        selDate <- as.Date(strChoice, origin = "1970-01-01")
        regData <- subsetData(data, selectedWeek = selDate, selectedAgeGroup = "All age groups")
      }
        
      ages <- regData$Age
      ages[which(ages == "< 1")] <- rep(0, length(which(ages == "< 1")))
      ages <- as.integer(ages)
      
      brks <- ageGroups
      ages <- cut(as.integer(ages), breaks = brks)
      df <- as.data.frame(table(ages))
      levels(df$ages) <- as.character(paste0("< ", brks[-1])) # TODO: actually, <=
      # df$ages <- factor(df$ages, levels = brks[-1])
      
      colors <- rep("rgba(0,0,0,0)", nrow(df))
      
      # Find the selected age group and colour it
      ages <- sapply(1:nrow(df), function(i) as.integer(strsplit(levels(df$ages)[i], " ")[[1]][2]))
      if (getSelectedAgeGroup() != "All age groups") {
        val <- as.integer(strsplit(getSelectedAgeGroup(), " ")[[1]][1]) - 1
        if (!is.na(val)) {
          if (val == -1)
            colors[which(ages == 1)] <- "red"
          else
            colors[which(ages == val) + 1] <- "red"
        }
      }
      
      p <- plot_ly(df, 
                   x = ~ages, y = ~Freq, type = "bar", # barmode = "stack", 
                   name = "Age", marker = list(color = "salmon", line = list(color = colors, width = 1.5)),
                   source = "agePlot"
      ) %>%
        layout(title = paste0("Age distr. of ", getSelectedSerotype(), 
                              ifelse(getSelectedSubtype() == "All subtypes", "", paste0(" (subtype ", getSelectedSubtype(), ")")),
                              " samples <br> in ", getSelectedRegion(), 
                              ifelse(getSelectedDepartment() == "All departments", "", paste0(" (dep ", getSelectedDepartment(), ")")),
                              " on",
                              ifelse(as.character(selDate) == "All weeks", " all weeks", paste0(" week starting ", selDate))), 
               titlefont = smallFont,
               showlegend = FALSE,
               xaxis = list(title = "Age", showgrid = FALSE, zeroline = FALSE, showticklabels = TRUE, titlefont = smallFont, tickfont = smallFont, tickangle = 0),
               yaxis = list(title = "Number of samples", showgrid = TRUE, zeroline = FALSE, showticklabels = TRUE, titlefont = smallFont, tickfont = smallFont))
      
      return(p)
    }
    
    #################################################################
    if (graphChoice == "reg") {
      if (strChoice == "")
        strChoice <- getSelectedRegion()
      
      regData <- subsetData(data, selectedRegion = strChoice, selectedDepartment = "All departments")
      
      df <- as.data.frame(table(regData$Departement))
      
      # # Exclude deps from other regions
      # deps <- unique(data[which(data$region == region), ]$Departement)
      # df <- df[-which(!(df$Var1 %in% deps)), ]
      
      # Exclude departments with zero cases
      df <- df[which(df$Freq != 0), ]
      
      if (length(df) == 0) {
        df <- data.frame(Var1 = 0, Freq = 0)
        colors <- c("rgba(0,0,0,0)")
      } else {
        df$Var1 <- factor(df$Var1, levels = unique(df$Var1))
        
        colors <- rep("rgba(0,0,0,0)", nrow(df))
        colors[which(df$Var1 == input$department)] <- "red"
        
        print(paste0("Dep graph update, colored: ", input$department))
      }
      
      p <- plot_ly(df, 
                   x = ~Var1, y = ~Freq, type = "bar",
                   marker = list(line = list(color = colors, width = 1.5)),
                   name = "Departements",
                   source = "depPlot"
      ) %>%
        layout(title = paste0("Dep. distr. of ", getSelectedSerotype(), 
                              ifelse(getSelectedSubtype() == "All subtypes", "", paste0(" (subtype ", getSelectedSubtype(), ")")),
                              " samples <br> in ", getSelectedRegion(), 
                              ifelse(getSelectedAgeGroup() == "All age groups", "", paste0(" (age ", getSelectedAgeGroup(), ")")),
                              " on week starting ", getSelectedWeek()),
               titlefont = smallFont,
               showlegend = FALSE,
               xaxis = list(title = "Departement", showgrid = FALSE, zeroline = FALSE, showticklabels = TRUE, titlefont = smallFont, tickfont = smallFont, tickangle = 0),
               yaxis = list(title = "Number of samples", showgrid = TRUE, zeroline = FALSE, showticklabels = TRUE, titlefont = smallFont, tickfont = smallFont))
      
      return(p)
    }
    
    plotly_empty(type = "scatter", mode = "markers") # TODO:
  })
  
  
  
  
  
  #######################################################################################################################
  
  
  # Projection plot
  output$projectionPlot <- renderPlotly({
    data <- loadedData()
    
    if (is.null(data)) {
      plotly_empty(type = "scatter", mode = "markers") %>%
        layout(title = "<- Please upload data", titlefont = list(
          # family = 'Courier New, monospace',
          size = 18, color = "red"#'#7f7f7f'
        ))
    } else {
      # data$region <- useNewRegions(as.vector(data$region))
      
      # print(sort(table(data$region), decreasing = TRUE)) # Print number of samples from each region
      
      h <- getHistForRegion(input$region, input$year, data, input$datesType)
      
      # print(h)
      
      p <- plot_ly(x = h$mids + 0.5, y = h$counts, type = 'bar', #line = list(color = '#45171D'),
                   hoverinfo = "text",
                   text = paste(h$counts, ' samples', " on week ", h$mids + 0.5, " (",
                                format(lastDayOfWeek(h$mids + 1.5, input$year) - 6, format = "%d-%b"), " to ",
                                format(lastDayOfWeek(h$mids + 1.5, input$year), format = "%d-%b"),
                                ")", sep = ""),
                   name = paste0(input$year, " (all samples: ", sum(h$counts), ")"),
                   source = "projPlot")
      
      # Histogram of samples with MLST results
      hm <- getHistForRegion(input$region, input$year, data[which(data$profil_MLST != ""), ], input$datesType)
      
      p <- add_trace(p, x = hm$mids + 0.5, y = hm$counts, type = 'bar', #line = list(color = '#45171D'),
                     hoverinfo = "text",
                     text = paste(hm$counts, ' samples', " on week ", hm$mids + 0.5, " (",
                                  format(lastDayOfWeek(hm$mids + 1.5, input$year) - 6, format = "%d-%b"), " to ",
                                  format(lastDayOfWeek(hm$mids + 1.5, input$year), format = "%d-%b"),
                                  ")", sep = ""),
                     name = paste0(input$year, " (WGS samples: ", sum(hm$counts), ")"))
      
      colorPal <- colorRamp(c("red", "orange", "yellow", "darkgreen", "blue"))(seq(0, 1,0.2))
      colorPal <- rgb(colorPal[, 1] / 255, colorPal[, 2] / 255, colorPal[, 3] / 255)
      
      for (yr in (as.integer(input$year) - 1):2010) {
        hPrev <- getHistForRegion(input$region, yr, data, input$datesType)
        
        do.call("<-", list(paste("h", yr, sep = ""), hPrev))
        
        p <- add_trace(p, x = hPrev$mids + 0.5, y = hPrev$counts, type = 'scatter', mode = 'lines', name = yr, #yaxis = 'y2',
                       line = list(color = colorPal[yr - 2010]), visible = "legendonly",
                       hoverinfo = "text", 
                       text = paste(hPrev$counts, 'samples'))
      }
      
      meanCounts <- apply(as.array(1:52), 1, function(x) mean(apply(as.array((as.integer(input$year) - 1):2010), 1, function(year) eval(as.name(paste("h", year, sep = "")))$counts[x]), na.rm = TRUE))
      print(meanCounts)
      p <- add_trace(p, x = 1:52 - 1, y = meanCounts, type = 'scatter', mode = 'lines', 
                     name = paste("Average ", 2010, "-", (as.integer(input$year) - 1), sep = ""), #yaxis = 'y2',
                     line = list(color = "black"),
                     hoverinfo = "text",
                     text = paste(signif(meanCounts, 2), 'samples on week', 1:52 - 1))
      
      
      
      
      
      curWeek <- data.table::isoweek(getSelectedWeek())
      
      
      s <- event_data("plotly_click", source = "projPlot")
      if (length(s)) {
        curWeek <- s$x[1]
        updateSelectInput(session, inputId = "week", selected = curWeek)
      }
      
      
      if (!is.na(curWeek)) {
        xLine <- as.integer(curWeek)
        p <- layout(p, shapes = list(type = 'line', 
                                     x0 = xLine, x1 = xLine, 
                                     y0 = 0, y1 = max(h$counts), 
                                     line = list(dash = 'dot', width = 2)))
        estVec <- c()
        for (i in curWeek:52) {
          if (i >= input$highSeasonRange[1] & i <= input$highSeasonRange[2])
            estVec <- c(estVec, floor(meanCounts[i] * input$coeffForHighSeason))
          else
            estVec <- c(estVec, floor(meanCounts[i] * input$coeffDefault))
        }
        
        # estimatedWGS <- floor(sum(meanCounts[curWeek:52]) * input$coeffDefault) + 
        #   floor(sum(meanCounts[input$highSeasonRange[1]:input$highSeasonRange[2]]) * (input$coeffForHighSeason - input$coeffDefault))
        estimatedWGS <- sum(estVec)
        
        p <- add_trace(p, x = as.integer(curWeek):52, y = estVec, type = 'scatter', mode = 'lines', 
                       name = paste("Projected", sep = ""), #yaxis = 'y2',
                       line = list(color = "orange"),
                       hoverinfo = "text",
                       text = estVec
        )
        
        p <- layout(p, annotations = list(
          x = curWeek,
          y = max(h$counts), 
          text = paste("Projected WGS: ", estimatedWGS, 
                       ", total: ", sum(hm$counts[1:curWeek]) + estimatedWGS, sep = ""),
          # font = smallFont,
          showarrow = FALSE,
          xanchor = "left"
        ))
        
        p <- layout(p, annotations = list(
          x = curWeek,
          y = max(h$counts) / 2, 
          text = paste("Actual: before: ", sum(hm$counts[1:curWeek]),
                       "(", signif(100 * sum(hm$counts[1:curWeek]) / sum(h$counts[1:curWeek]), 2), "%)", 
                       ", after: ", sum(hm$counts[(curWeek + 1):length(hm$counts)]), 
                       "(", signif(100 * sum(hm$counts[(curWeek + 1):length(hm$counts)]) / sum(h$counts[(curWeek + 1):length(hm$counts)]), 2), "%)",
                       sep = ""),
          # font = smallFont,
          showarrow = FALSE,
          xanchor = "left"
        ))
        
      }
      
      # print(as.integer(input$highSeasonRange[1]):as.integer(input$highSeasonRange[2]))
      # print(rep(max(h$counts) / 3, as.integer(input$highSeasonRange[2]) - as.integer(input$highSeasonRange[1]) + 1))
      
      p <- add_markers(p, x = as.integer(input$highSeasonRange[1]):as.integer(input$highSeasonRange[2]),
                       y = rep(floor(max(h$counts) / 3), as.integer(input$highSeasonRange[2]) - as.integer(input$highSeasonRange[1]) + 1),
                       type = 'scatter', mode = 'markers',
                       marker = list(color = "red"),
                       name = "High season",
                       inherit = FALSE,
                       hoverinfo = "text",
                       text = paste("Week ", as.integer(input$highSeasonRange[1]):as.integer(input$highSeasonRange[2]),
                                    # format(as.Date(h$breaks[-length(h$breaks)] + 1, origin = "1970-01-01"), format = "%d-%b"), " to ",
                                    # format(as.Date(h$breaks[-length(h$breaks)] + 7, origin = "1970-01-01"), format = "%d-%b"),
                                    sep = "")
      )
      
      p <- layout(p, title = paste0("Projection for ", input$region, " ", " samples in ", input$year, sep = ""), titlefont = smallFont,
                  legend = list(font = smallFont),
                  xaxis = list(title = "Week number", showgrid = FALSE, zeroline = FALSE, showticklabels = TRUE, tickfont = smallFont, titlefont = smallFont),
                  yaxis = list(side = 'left', title = 'Number of samples', showgrid = TRUE, zeroline = FALSE, showticklabels = TRUE, tickfont = smallFont, titlefont = smallFont))
      return(p)
    }
  })
  
  
  
  
  
  
  
  
  #######################################################################################################################
  
  findClusters <- function(dates, window, threshold) {
    if (length(dates) < threshold)
      return()
    
    dates <- sort(dates, decreasing = FALSE)
    
    starts <- c()
    for (i in 1:(length(dates) - threshold + 1)) {
      if (dates[i + threshold - 1] - dates[i] < window)
        starts <- c(starts, dates[i])
    }
    
    return(starts)
  }
  
  output$clusterPlot <- renderPlot({
    sData <- loadedData()
    
    if (is.null(sData)) {
      plot.new()
      title("Please upload data", col.main = "red")
    } else {
      year <- as.integer(input$year) ########################
      
      
      if (input$datesType == "Date of reception")
        timeData <- sData[which(sData$Date.de.reception >= as.Date(paste(year, "-01-01", sep = "")) & 
                                sData$Date.de.reception < as.Date(paste(year + 1, "-01-01", sep = ""))), ]
      else
        timeData <- sData[which(sData$Date.isolement >= as.Date(paste(year, "-01-01", sep = "")) & 
                                sData$Date.isolement < as.Date(paste(year + 1, "-01-01", sep = ""))), ]
      
      slices <- c()
      CTs <- c()
      for (CT in names(sort(table(timeData$CRISPOL.type, exclude = ""), decreasing = TRUE))) {
        xData <- timeData[which(timeData$CRISPOL.type == CT), ]
        xData$Departement <- as.integer(as.character(xData$Departement))
        sumStarts <- 0
        for (dep in sort(unique(xData$Departement))) {
          if (input$datesType == "Date of reception")
            dates <- xData$Date.de.reception[which(xData$Departement == dep)]
          else
            dates <- xData$Date.isolement[which(xData$Departement == dep)]
          
          window <- input$window # 30
          threshold <- input$threshold # 3
          starts <- findClusters(dates, window, threshold)
          if (length(starts) > 0) {
            # points(CT, length(starts))
            sumStarts <- sumStarts + length(starts)
          }
        }
        if (sumStarts > 0) {
          slices <- c(slices, sumStarts)
          CTs <- c(CTs, CT)
        }
      }
      
      # pie(slices, labels = CTs, #col = rainbow(length(lbls)),
      #     main = "Number of clusters per CT")
      
      # barplot(slices, main = "Number of clusters per CT", # horiz = TRUE,
      #         names.arg = CTs)
      
      # CT <- as.integer(input$CT) #############################
      
      s <- event_data("plotly_click", source = "clustPlot")
      
      CT <- CTs[1]
      
      if (length(s)) {
        # print(s)
        # print(s$x)
        CT <- s$x
      }
      
      xData <- timeData[which(timeData$CRISPOL.type == CT), ]
      
      
      xData$Departement <- as.integer(as.character(xData$Departement))
      
      if (input$datesType == "Date of reception")
        xDates <- xData$Date.de.reception
      else
        xDates <- xData$Date.isolement
      
      plot(xData$Departement, xDates, main = paste("CT_", CT, " in ", year, sep = ""), 
           xlab = "Department", ylab = "Date of isolation", las = 1,
           xlim = c(0, 100),
           pch = 3, cex = 1)
      abline(v = unique(xData$Departement), col = "gray", lty = 3)
      # points(xData$Departement, xData$Date.isolement, pch = 15, cex = 0.5)
      points(xData$Departement, xDates, pch = 3, cex = 1)
      axis(side = 1, at = 0:100, labels = FALSE, tck = -0.01)
      
      for (dep in sort(unique(xData$Departement))) {
        dates <- xDates[which(xData$Departement == dep)]
        window <- input$window # 30
        threshold <- input$threshold # 3
        starts <- findClusters(dates, window, threshold)
        if (length(starts) > 0) {
          for (i in 1:length(starts))
            lines(c(dep, dep), c(starts[i], starts[i] + window), col = rgb(1, 0.2, 0.1, 0.3), lwd = 5)
        }
      }
    }
    
  })
  
  
  
  output$clusterHist <- renderPlotly({
    sData <- loadedData()
    
    if (is.null(sData)) {
      plotly_empty(type = "scatter", mode = "markers") %>%
        layout(title = "<- Please upload data", titlefont = list(
          # family='Courier New, monospace',
          size = 18,
          color = "red"#'#7f7f7f'
        ))
    } else {
      year <- as.integer(input$year) ########################
      
      if (input$datesType == "Date of reception")
        timeData <- sData[which(sData$Date.de.reception >= as.Date(paste(year, "-01-01", sep = "")) & 
                                  sData$Date.de.reception < as.Date(paste(year + 1, "-01-01", sep = ""))), ]
      else
        timeData <- sData[which(sData$Date.isolement >= as.Date(paste(year, "-01-01", sep = "")) & 
                                sData$Date.isolement < as.Date(paste(year + 1, "-01-01", sep = ""))), ]
      slices <- c()
      CTs <- c()
      for (CT in names(sort(table(timeData$CRISPOL.type, exclude = ""), decreasing = TRUE))) {
        print(CT)
        
        xData <- timeData[which(timeData$CRISPOL.type == CT), ]
        xData$Departement <- as.integer(as.character(xData$Departement)) ##########
        sumStarts <- 0
        for (dep in sort(unique(xData$Departement))) {
          if (input$datesType == "Date of reception")
            dates <- xData$Date.de.reception[which(xData$Departement == dep)]
          else
            dates <- xData$Date.isolement[which(xData$Departement == dep)]
          
          window <- input$window # 30
          threshold <- input$threshold # 3
          starts <- findClusters(dates, window, threshold)
          if (length(starts) > 0) {
            # points(CT, length(starts))
            sumStarts <- sumStarts + length(starts)
          }
        }
        if (sumStarts > 0) {
          slices <- c(slices, sumStarts)
          CTs <- c(CTs, CT)
        }
      }
      
      # pie(slices, labels = CTs, #col = rainbow(length(lbls)),
      #     main = "Number of clusters per CT")
      
      # CTs <- CTs[order(as.integer(CTs))]
      # slices <- slices[order(as.integer(CTs))]
      
      
      selectedCT <- 1
      s <- event_data("plotly_click", source = "clustPlot")
      
      if (length(s)) {
        selectedCT <- s$x
      }
      
      
      colors <- rep("rgba(0,0,0,0)", length(CTs))
      colors[which(CTs == selectedCT)] <- "red"
      
      p <- plot_ly(x = CTs, y = slices, type = 'bar', #line = list(color = '#45171D'),
                   marker = list(#color = "lightblue", 
                     line = list(color = colors, width = 1.5)),
                   hoverinfo = "text",
                   text = paste0("CT-", CTs, ": ", slices, " clusters"),
                   name = "Clusters",
                   source = "clustPlot") %>%
        layout(title = "", 
               xaxis = list(title = "CT", showgrid = FALSE, zeroline = FALSE, showticklabels = TRUE, tickfont = smallFont, titlefont = smallFont), 
               yaxis = list(title = "Number of clusters", showgrid = TRUE, zeroline = FALSE, showticklabels = TRUE, tickfont = smallFont, titlefont = smallFont))
      
      return(p)
    }
    
  })
  
  
  
  output$clusterPlotly <- renderPlotly({
    sData <- loadedData()
    
    if (is.null(sData)) {
      # plot.new()
      # title("Please upload data", col.main = "red")
      plotly_empty(type = "scatter", mode = "markers") %>%
        layout(title = "", titlefont = list(
          # family='Courier New, monospace',
          size = 18,
          color = "red"#'#7f7f7f'
        ))
    } else {
      year <- as.integer(input$year) ########################
      
      
      
      if (input$datesType == "Date of reception")
        timeData <- sData[which(sData$Date.de.reception >= as.Date(paste(year, "-01-01", sep = "")) & 
                                  sData$Date.de.reception < as.Date(paste(year + 1, "-01-01", sep = ""))), ]
      else
        timeData <- sData[which(sData$Date.isolement >= as.Date(paste(year, "-01-01", sep = "")) & 
                                  sData$Date.isolement < as.Date(paste(year + 1, "-01-01", sep = ""))), ]
      
      slices <- c()
      CTs <- c()
      for (CT in names(sort(table(timeData$CRISPOL.type, exclude = ""), decreasing = TRUE))) {
        print(CT)
        
        xData <- timeData[which(timeData$CRISPOL.type == CT), ]
        xData$Departement <- as.integer(as.character(xData$Departement)) ##################
        sumStarts <- 0
        
        for (dep in sort(unique(xData$Departement))) {
          if (input$datesType == "Date of reception")
            dates <- xData$Date.de.reception[which(xData$Departement == dep)]
          else
            dates <- xData$Date.isolement[which(xData$Departement == dep)]
          
          window <- input$window # 30
          threshold <- input$threshold # 3
          starts <- findClusters(dates, window, threshold)
          if (length(starts) > 0) {
            # points(CT, length(starts))
            sumStarts <- sumStarts + length(starts)
          }
        }
        if (sumStarts > 0) {
          slices <- c(slices, sumStarts)
          CTs <- c(CTs, CT)
        }
      }
      # pie(slices, labels = CTs, #col = rainbow(length(lbls)),
      #     main = "Number of clusters per CT")
      
      # barplot(slices, main = "Number of clusters per CT", # horiz = TRUE,
      #         names.arg = CTs)
      
      # CT <- as.integer(input$CT) #############################
      
      s <- event_data("plotly_click", source = "clustPlot")
      
      CT <- CTs[1]
      
      if (length(s)) {
        CT <- s$x
      }
      
      
      xData <- timeData[which(timeData$CRISPOL.type == CT), ]
      xData$Departement <- as.integer(as.character(xData$Departement))
      
      # plot(xData$Departement, xData$Date.isolement, main = paste("CT_", CT, " in ", year, sep = ""), 
      #      xlab = "Department", ylab = "Date of isolation", las = 1,
      #      xlim = c(0, 100),
      #      pch = 3, cex = 1)
      # abline(v = unique(xData$Departement), col = "gray", lty = 3)
      # 
      # # points(xData$Departement, xData$Date.isolement, pch = 15, cex = 0.5)
      # points(xData$Departement, xData$Date.isolement, pch = 3, cex = 1)
      # 
      # axis(side = 1, at = 0:100, labels = FALSE, tck = -0.01)
      
      
      xData <- xData[which(xData$Departement < 100), ]
      
      if (input$datesType == "Date of reception")
        xDates <- xData$Date.de.reception
      else
        xDates <- xData$Date.isolement
      
      if (nrow(xData) != 0) {
        p <- plot_ly(x = xData$Departement, 
                     y = as.POSIXct(xDates) + 12 * 60 * 60 * (runif(length(xDates), 0, 1)), # Add jitter to dates
                     type = 'scatter', mode = "markers", 
                     marker = list(color = "black", size = 6, symbol = 22),
                     #line = list(color = '#45171D'),
                     hoverinfo = "text",
                     text = paste0("Department: ", xData$Departement, ", date: ", xDates, ""),
                     name = "Samples (isolated)",
                     source = "clustPlotly") %>%
          layout(title = paste("CT_", CT, " in ", year, sep = ""), titlefont = smallFont,
                 legend = list(font = smallFont),
                 xaxis = list(title = "Department", showgrid = TRUE, zeroline = FALSE, showticklabels = TRUE, tickfont = smallFont, titlefont = smallFont), 
                 yaxis = list(title = "", #"Date of isolation", 
                              showgrid = TRUE, zeroline = FALSE, showticklabels = TRUE, tickfont = smallFont, titlefont = smallFont))
        
        
        clustN <- 0
        for (dep in sort(unique(xData$Departement))) {
          dates <- xDates[which(xData$Departement == dep)]
          window <- input$window # 30
          threshold <- input$threshold # 3
          starts <- findClusters(dates, window, threshold)
          if (length(starts) > 0) {
            for (i in 1:length(starts)) {
              len <- window
              len <- max(dates[which(dates <= starts[i] + window)]) - starts[i]
              
              clustN <- clustN + 1
              p <- add_lines(p, c(dep, dep), c(as.Date(starts[i], origin = "1970-01-01"), 
                                               as.Date(starts[i] + len + 1, origin = "1970-01-01")), 
                             # col = rgb(1, 0.2, 0.1, 0.3), lwd = 5, 
                             name = paste0("Cluster ", clustN, " (dept. ", dep, ")"),
                             line = list(width = 3),
                             hoverinfo = "none",
                             inherit = FALSE)
            }
          }
        }
      } else {
        # TODO: remove duplicated code
        p <- plot_ly(x = c(), 
                     y = c(),
                     type = 'scatter', mode = "markers", 
                     marker = list(color = "black", size = 6, symbol = 22),
                     #line = list(color = '#45171D'),
                     # hoverinfo = "text",
                     # text = paste0("Department: ", xData$Departement, ", date: ", xData$Date.isolement, ""),
                     name = "Samples (isolated)",
                     source = "clustPlotly") %>%
          layout(#title = paste("CT_", CT, " in ", year, sep = ""), 
            titlefont = smallFont,
            legend = list(font = smallFont),
            xaxis = list(title = "Department", showgrid = TRUE, zeroline = FALSE, showticklabels = TRUE, tickfont = smallFont, titlefont = smallFont), 
            yaxis = list(title = "", showgrid = TRUE, zeroline = FALSE, showticklabels = TRUE, tickfont = smallFont, titlefont = smallFont))
        
      }
      
      
      return(p)
    }
    
  })
  
  # Info text
  output$info <- renderText({
    data <- loadedData()
    
    if (is.null(data)) {
      return("")
    } else {
      selDate <- getSelectedWeek()
      
      if (is.null(selDate))
        return("")
               
      if (as.character(selDate) == "All weeks")
        selDate <- "All weeks"
      
      s <- event_data("plotly_hover", source = "mainPlot")
      if (length(s)) {
        selDate <- as.Date(s$x, origin = "1970-01-01")
        extraInfo <- "See age distribution"
      }
      
      selRegion <- getSelectedRegion()
      s <- event_data("plotly_hover", source = "regPlot")
      if (length(s)) {
        selRegion <- s$x
        extraInfo <- "See departement distribution"
      }
      
      selSerotype <- getSelectedSerotype()
      s <- event_data("plotly_hover", source = "serPlot")
      if (length(s)) {
        selSerotype <- s$x
        extraInfo <- "See MLST/CT distribution"
      }
      
      print(paste0("sel dtae: ", selDate))
      
      # Subset data
      regData <- subsetData(data, selectedRegion = selRegion,
                            selectedWeek = selDate,
                            selectedSerotype = selSerotype)
      
      
      # selData <- regData[, c(1,
      #                        3, # Serotype_AGG
      #                        5, # Serotype_FINAL
      #                        6, 7, # MLST, CT
      #                        8, 9, 10, # dates iso-rec-val 
      #                        11, # Region
      #                        12, # Voayge
      #                        14, # Age
      #                        18, # Departement
      #                        20 # Epidemie
      # )]
      
      selData <- regData[, c("Numero.de.dossier", "Serotype_AGG", "Serotype_FINAL", "profil_MLST", 
                             "CRISPOL.type", "Date.isolement", "Date.de.reception", "Date.de.validation",
                             "region", "Voyage", "Age", "Departement", "Epidemie")]
      
      if (as.character(selDate) == "All weeks")
        output$table <- renderDataTable(selData[c(), ], options = list(#lengthMenu = c("All", 5, 10, 50, 100), 
          # scrollX = TRUE,
          pageLength = -1) # for renderDataTable only
        )
      else
      output$table <- renderDataTable(selData, options = list(#lengthMenu = c("All", 5, 10, 50, 100), 
        # scrollX = TRUE,
        pageLength = -1) # for renderDataTable only
      )
      
      extraInfo <- ""
      
      if (as.character(selDate) == "All weeks")
        numInfo <- paste0(nrow(regData), " samples (details are not shown)")
      else
        numInfo <- paste0(nrow(regData), " samples (see details below)")
      
      # dataInfo <- paste(head(regData[, c(5, 6, 7, 9, 11, 14)]), collapse = "\n")
      dataInfo <- ""
      
      res <- paste0("Selected week: ", ifelse(as.character(selDate) == "All weeks", 
                                              "all weeks", paste0(data.table::isoweek(selDate), " (", selDate, " - ", selDate + 6, ")")), "\n",
                    "Selected region: ", selRegion, "\n",
                    "Selected department: ", getSelectedDepartment(), "\n",
                    "Selected serotype: ", selSerotype, "\n",
                    "Selected subtype: ", getSelectedSubtype(), "\n",
                    "Selected age group: ", getSelectedAgeGroup(), "\n",
                    # extraInfo, "\n",
                    numInfo, "\n"
                    # dataInfo, "\n\n\n\n\n\n\n\n"
      )
      return(res)
    }
  })
  
  
  
  
  ###########################################################################################################################
  
  output$monthlyHist <- renderPlotly({
    data <- loadedData()
    
    if (is.null(data)) {
      plotly_empty(type = "scatter", mode = "markers") %>%
        layout(title = "<- Please upload data", titlefont = list(
          # family = 'Courier New, monospace',
          size = 18, color = "red"#'#7f7f7f'
        ))
    } else {
      xdata <- subsetData(data, selectedRegion = getSelectedRegion(), selectedWeek = "All weeks", 
                         selectedSerotype = getSelectedSerotype())
      
      h <- getMonthlyHist(xdata, input$datesType)
      
      # print(h)
      
      p <- plot_ly(x = h$mids + 0.5, y = h$counts, type = 'bar', #line = list(color = '#45171D'),
                   hoverinfo = "text",
                   text = paste(h$counts, ' samples', " in ", month.name[h$mids + 0.5], " ", input$year, sep = ""),
                   name = input$year,
                   source = "monthlyPlot")
      
      colorPal <- colorRamp(c("red", "orange", "yellow", "darkgreen", "blue"))(seq(0, 1, 0.2))
      colorPal <- rgb(colorPal[, 1] / 255, colorPal[, 2] / 255, colorPal[, 3] / 255)
      
      for (yr in (as.integer(input$year) - 1):2010) {
        xdata <- subsetData(data, selectedRegion = getSelectedRegion(), selectedWeek = "All weeks", 
                           selectedSerotype = getSelectedSerotype(), selectedYear = yr)
        hPrev <- getMonthlyHist(xdata, input$datesType)
        
        do.call("<-", list(paste("h", yr, sep = ""), hPrev))
        
        p <- add_trace(p, x = hPrev$mids + 0.5, y = hPrev$counts, 
                       # type = 'scatter', mode = 'lines',  
                       type = 'bar', 
                       name = yr, #yaxis = 'y2',
                       line = list(color = colorPal[yr - 2010]), visible = "legendonly",
                       hoverinfo = "text", 
                       text = paste0(hPrev$counts, ' samples in ', month.name[hPrev$mids + 0.5], " ", yr))
      }
      
      meanCounts <- apply(as.array(1:12), 1, function(x) mean(apply(as.array((as.integer(input$year) - 1):2010), 1, function(year) eval(as.name(paste("h", year, sep = "")))$counts[x]), na.rm = TRUE))
      print(meanCounts)
      p <- add_trace(p, x = 1:12, y = meanCounts, type = 'scatter', mode = 'lines', 
                     name = paste("Average ", 2010, "-", (as.integer(input$year) - 1), sep = ""), #yaxis = 'y2',
                     line = list(color = "black"),
                     hoverinfo = "text",
                     text = paste("Average", signif(meanCounts, 2), 'samples in', month.name[1:12]))
      
      p <- layout(p, title = paste0("Monthly case counts of ", getSelectedSerotype(), " for ", 
                                    input$region, " ", " samples in ", input$year, sep = ""), titlefont = smallFont,
                  legend = list(font = smallFont),
                  xaxis = list(title = "Month", type = "category", showgrid = FALSE, zeroline = FALSE, showticklabels = TRUE, tickfont = smallFont, titlefont = smallFont),
                  yaxis = list(side = 'left', title = 'Number of samples', showgrid = TRUE, zeroline = FALSE, showticklabels = TRUE, tickfont = smallFont, titlefont = smallFont))
      return(p)
    }
  })
})