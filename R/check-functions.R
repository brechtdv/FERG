###=========================================================================#
### Functions to support Calculation Tool Activity
###
### Start  08.12.2013 Brecht
### Update 14.02.2014 add 'check_db'
### Update 17.02.2014 add 'check_convergence'
### Update 13.07.2014 'check_convergence' understands regions
### Update 13.07.2014 'check_convergence' checks median & upper/lower UL
### Update 13.07.2014 'cumulative_moving_fun' replaces '.._average'
###=========================================================================#

###=========================================================================#
###== FUNCTIONS ============================================================#
###-- check_XL ...................... check if sheets are correctly ordered
###-- check_db ...................... check raw imputation database
###-- check_convergence ............. check MC convergence
###---| cumulative_moving_fun ....... calculate CM function


##--------------------------------------------------------------------------#
## Check if sheets are correctly ordered -----------------------------------#

check_XL <- 
function(XL){
  ## extract sheet names & disease model
  sheets <- getSheets(XL)
  DM <- readWorksheet(XL, sheet = 1, header = TRUE)

  ## first sheets should be DM and KEY
  if (sheets[1] != "DM")
    stop("First sheet is not a 'DM' sheet.")
  if (sheets[2] != "KEY")
    stop("Second sheet is not a 'KEY' sheet.")

  ## sheets should follow DM order
  if (!all(sheets[-c(1,2)] == DM[, 1]))
    stop("Sheet names do not correspond to DM.")

  ## if all checks were successful
  message("All checks were successful.")
}


##--------------------------------------------------------------------------#
## Check raw imputation database -------------------------------------------#

check_db <-
function(db) {
  ## check length of object
  if (length(db) != 12) 
    stop("INC data.frame does not contain 12 elements")

  ## 'Inc' shoudl be numeric
  if (!is.numeric(db$Inc))
    stop("'db$Inc' is not numeric.")

  ## check if countries and subregions are specified in correct order
  if (!all(db$country == crpop_2005$Country))
    stop("'db$country' does not correspond to 'crpop_2005$Country'.")
  if (!all(db$WHOsub ==
           with(crpop_2005, paste(WHORegion, WHOsubRegion, sep = ""))))
    stop("'db$WHOsub' does not correspond to 'crpop_2005'.")

  ## if all checks were successful
  message("All checks were successful.")
}


##--------------------------------------------------------------------------#
## CHECK MC CONVERGENCE ----------------------------------------------------#

check_convergence <-
function(x, n = 100L, denom = 1e5) {
  daly <- 0
  n_nodes <- length(x[[1]][[1]]@samples)
  n_iterations <- nrow(x[[2]][[1]])
  contrib <- get_contrib(x[[1]])

  n_regions <- length(x[[2]])

  for (i in seq(n_regions)) {
    for (j in seq(n_nodes)) {
      daly <-
        daly +
        rowSums(x[[1]][[i]]@samples[[j]][[contrib[j]]])
    }
  }

  daly <- daly / denom

  cm_mean <- cumulative_moving_fun(daly, mean)
  cm_med <- cumulative_moving_fun(daly, median)
  cm_lwr <- cumulative_moving_fun(daly, function(x) quantile(x, 0.025))
  cm_upr <- cumulative_moving_fun(daly, function(x) quantile(x, 0.975))

  cm_mean_ad <- abs_dev(cm_mean)
  cm_med_ad <- abs_dev(cm_med)
  cm_lwr_ad <- abs_dev(cm_lwr)
  cm_upr_ad <- abs_dev(cm_upr)

  par(mfrow = c(2, 1))
  par(mar = c(4, 4, 0, 0) + .5)
  plot_mx <-
    cbind(tail(cm_upr_ad, n),
          tail(cm_lwr_ad, n),
          tail(cm_mean_ad, n),
          tail(cm_med_ad, n))
  matplot(plot_mx,
          ylim = c(0, max(0.01, max(plot_mx))),
          ylab = "absolute deviance", xlab = "iteration",
          type = "l", lty = 1, axes = FALSE)
  add_to_ticks <- n - tail(axTicks(1), 1)
  axis(1,
       at = c(0, axTicks(1) + add_to_ticks),
       labels = c(0, axTicks(1) + add_to_ticks) +
                n_iterations - n)
  axis(2)
  box()
  abline(h = 0.01, col = "grey", lwd = 2)
  legend("topright",
         legend = c("mean", "median", "2.5%", "97.5%"),
         col = c(3, 4, 2, 1), lty = 1)

  plot_mx <-
    cbind(tail(cm_upr, n),
          tail(cm_lwr, n),
          tail(cm_mean, n),
          tail(cm_med, n))
  matplot(plot_mx,
          ylab = "DALY", xlab = "iteration",
          type = "l", lty = 1, log = "y", axes = FALSE)
  add_to_ticks <- n - tail(axTicks(1), 1)
  axis(1,
       at = c(0, axTicks(1) + add_to_ticks),
       labels = c(0, axTicks(1) + add_to_ticks) +
                n_iterations - n)
  axis(2)
  box()
  legend("topright",
         legend = c("mean", "median", "2.5%", "97.5%"),
         col = c(3, 4, 2, 1), lty = 1)

  out <-
    data.frame(mean = c(mean(tail(cm_mean_ad, n)),
                        mean(tail(cm_med_ad, n)),
                        mean(tail(cm_lwr_ad, n)),
                        mean(tail(cm_upr_ad, n))),
               max  = c(max(tail(cm_mean_ad, n)),
                        max(tail(cm_med_ad, n)),
                        max(tail(cm_lwr_ad, n)),
                        max(tail(cm_upr_ad, n))))
  rownames(out) <- c("mean", "median", "2.5%", "97.5%")
  return(round(t(out), 4))
}

cumulative_moving_fun <-
function(x, fun) {
  cmf <-
    mapply(function(a, b) fun(a[1:b]),
           b = seq_along(x),
           MoreArgs = list(a = x))
  return(cmf)
}

abs_dev <-
function(x) {
  abs(tail(x, -1) - head(x, -1)) / tail(x, -1)
}
