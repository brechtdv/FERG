###=========================================================================#
### Functions to support Data Adjustment Activity
###
### Start  05.02.2014 Scott
### Update 14.02.2014 'impute' now also returns 'tau', 'mu.c', 'mu.tau'
### Update 14.02.2014 'extractIncidence' now recognizes 'sheet'
### Update 15.02.2014 'impute' now generates diagnostic plots
### Update 06.08.2014 'parseDistrDA' now recognizes U(2.5;97.5)
### Update 06.08.2014 'impute' gains options 'burnin', 'update'
###=========================================================================#

###=========================================================================#
###== FUNCTIONS ============================================================#
###-- extractIncidence ...... get incidence data from Excel database file
###---| readDatabase ........ read in Excel workbook using XLConnect
###-----| parseDistrDA ...... calculate mean incidence
###-- impute ................ impute missing values based on a RE model
###-- merge ................. merge observed with imputed values
###-- save_merged_db ........ save merged values to workbook


##--------------------------------------------------------------------------#
## GET INCIDENCE DATA FROM EXCEL FOR ALL COUNTRIES -------------------------#

extractIncidence <-
function(file, sheet) {
  ## read entire excel file
  wb <- readDatabase(file = file)

  ## check whether specified 'sheet' exists
  if (!(sheet %in% getSheets(wb))) 
    stop("sheet ", sQuote(sheet), " not recognised.")

  ## extract specified distribution
  dist <-
    readWorksheet(wb, sheet,
                  startRow = 7, endRow = 7, 
                  startCol = 14, endCol = 14,
                  header = FALSE)

  ## distribution must be 'Percentiles'
  if (dist != "Percentiles")
    stop("Distribution ", sQuote(dist), " not allowed.")

  ## extract data from sheet
  sheet_data <-
    readWorksheet(wb, sheet,
                  startRow = 8, endRow = 200, 
                  startCol = 2, endCol = 16)

  ## organize data into data.frame
  db192 <- sheet_data[, c(1, 7:15)]
  db192$WHOsub <- paste(sheet_data[, 2], sheet_data[, 3], sep = "")
  db192$Inc <- NA
  names(db192) <- c("country", paste("c", 1:9, sep = ""), "WHOsub", "Inc")

  ## identify missing values
  is_na <- apply(db192[, 2:10], 1, function(x) all(is.na(x)))

  ## calculate mean for countries with values
  for (i in seq(NumCountries_2005)[!is_na]) {
    db192$Inc[i] <- parseDistrDA(db192[i, 2:10], "Percentiles")[[2]]
  }

  ## reorganise 'db'
  db192 <- db192[, c("country", "WHOsub", "Inc", paste("c", 1:9, sep = ""))]

  ## seems to be necessary for some agent-databases...
  db192 <- transform(db192, Inc = as.numeric(db192$Inc))

  ## return data.frame
  return(db192)
}


##--------------------------------------------------------------------------#
## EXTRACT DESIRED DISTRIBUTION PARAMETERS ---------------------------------#

## Note: differs from parseDistr() in CalcToolActivity
parseDistrDA <-
function(cells, disttype) {
  # replace any ',' by '.' and convert to numeric
  cells <- as.numeric(gsub(",", ".", cells))

  # 0,5,25,50,75,95,100% all filled in  NORMALLY WOULD BE "CumulProb"
  # THIS IS ONLY TEMP SOL'N, FOR DATA ADJ ACTIVITY...
  if (disttype == "Percentiles" &
             !is.na(cells[1]) & !is.na(cells[3]) & !is.na(cells[4]) &
             !is.na(cells[5]) & !is.na(cells[6]) & !is.na(cells[7]) &
             !is.na(cells[9])){
    dist <- "Fixed"
    params <- c(cells[5])

  # 0,2.5,5,25,50,75,95,97.7,100% all filled in
  # THIS IS ONLY TEMP SOL'N, FOR DATA ADJ ACTIVITY...
  } else if (disttype == "Percentiles" &
             !is.na(cells[1]) & !is.na(cells[3]) & !is.na(cells[4]) &
             !is.na(cells[5]) & !is.na(cells[6]) & !is.na(cells[7]) &
             !is.na(cells[9]) & !is.na(cells[2]) & !is.na(cells[8])){
    dist <- "Fixed"
    params <- c(cells[5])

  # 2.5, 50 and 97.5% filled in
  # THIS IS ONLY TEMP SOL'N, FOR DATA ADJ ACTIVITY (b_cereus)
  } else if (disttype == "Percentiles" &
             !is.na(cells[2]) & !is.na(cells[8]) & !is.na(cells[5])){
    dist <- "Fixed"
    params <- c(cells[5])

  # 5, 50 and 95% filled in
  # THIS IS ONLY TEMP SOL'N, FOR DATA ADJ ACTIVITY (b_cereus)
  } else if (disttype == "Percentiles" &
             !is.na(cells[3]) & !is.na(cells[7]) & !is.na(cells[5])){
    dist <- "Fixed"
    params <- c(cells[5])

  # 50% only filled in
  # THIS IS ONLY TEMP SOL'N, FOR DATA ADJ ACTIVITY (Peanut)
  } else if (disttype == "Percentiles" &
             !is.na(cells[5])){
    dist <- "Fixed"
    params <- c(cells[5])

  # 5 and 95% filled in
  # TEMP SOL'N, FOR DATA ADJ ACTIVITY (Brucella)
  } else if (disttype == "Percentiles" &
             !is.na(cells[3]) & !is.na(cells[7])){
    dist <- "Fixed"
    params <- c(mean(c(as.numeric(cells[3]), as.numeric(cells[7]))))

  # 2.5 and 97.5% filled in
  } else if (disttype == "Percentiles" &
             !is.na(cells[2]) & !is.na(cells[8])){
    dist <- "Fixed"
    params <- c(mean(c(as.numeric(cells[2]), as.numeric(cells[8]))))

  # 0 and 100% filled in
  # TEMP SOL'N, FOR DATA ADJ ACTIVITY (Peanut)
  } else if (disttype == "Percentiles" &
             !is.na(cells[1]) & !is.na(cells[9])){
    dist <- "Fixed"
    params <- c(mean(c(as.numeric(cells[1]), as.numeric(cells[9]))))

  # if something else, return error
  } else {
    stop("Distribution not recognised.")
  }

  return(list(dist, as.numeric(params)))
}


##--------------------------------------------------------------------------#
## IMPUTE ------------------------------------------------------------------#

impute <-
function(inc, WHOsub, model = "lognormal.txt",
         inits = NULL, burnin = 5000, update = 5000) {
  ## define 'data' list
  data <-
    list(N = length(inc),
         NR = length(levels(WHOsub)),
         reg = WHOsub,
         log_y = log(inc))

  ## define & initialize model
  mod <-
    jags.model(file = system.file(model, package = "FERG"),
               data = data,
               inits = inits,
               n.chains = 2)

  ## burn-in
  update(mod, n.iter = burnin)

  ## get samples
  samples <-
    coda.samples(mod, c("log_y", "tau", "mu.c", "mu.tau"), n.iter = update)

  ## diagnostic plots
  dev.new(8, 5.5)
  par(mfrow = c(2, 3))
  traceplot(samples[, c("tau", "mu.tau", "mu.c")])
  densplot(samples[, c("tau", "mu.tau", "mu.c")])

  ## remove existing data
  is_na <- is.na(inc)  
  predictions <- samples[, c(is_na, rep(FALSE, 3))]

  ## BGR diagnostic
  cat("\nMultivariate psrf\n")
  cat(gelman.diag(samples[, c(is_na, rep(TRUE, 3))])$mpsrf, "\n")

  ## merge chains
  log_y <- rbind(predictions[[1]], predictions[[2]])
  y <- exp(log_y)

  ## return predictions and parameters
  return(list(y      = y,
              tau    = unlist(samples[, "tau"]),
              mu.c   = unlist(samples[, "mu.c"]),
              mu.tau = unlist(samples[, "mu.tau"])))
}


##--------------------------------------------------------------------------#
## MERGE -------------------------------------------------------------------#

merge <-
function(db, db_imputed) {
  ## identify missing values
  is_na <- is.na(db$Inc)

  ## prepare merged matrix
  db_merged <- sapply(db[, 4:12], as.numeric)

  ## calculate quantiles for log_y samples
  quantiles_y <-
    t(apply(db_imputed$y, 2,
            function(x) {
              quantile(x, c(0.025, 0.05, .25, .50, .75, .95, .975))
            }))

  ## merge data and fitted values
  db_merged[is_na, ] <- cbind(NA, quantiles_y, NA)

  ## return merged matrix
  return(db_merged)
}


##--------------------------------------------------------------------------#
## SAVE MERGED VALUES ------------------------------------------------------#

save_merged_db <-
function(db_merged, file, sheet) {
  wb <- loadWorkbook(file)
  if (!(sheet %in% getSheets(wb)))
    stop("sheet ", sQuote(sheet), " not recognised.")
  writeWorksheet(wb, db_merged, sheet,
                 startRow = 9, startCol = 8, header = FALSE, rownames = NULL)
  saveWorkbook(wb)
}
