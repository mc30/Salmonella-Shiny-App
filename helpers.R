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