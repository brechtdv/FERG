###=========================================================================#
### Source Attribution functions
### First  28.07.2014 Brecht
### Update 03.08.2014 Add 'incidence' to 'apply_prop_fb()' 
### Update 03.08.2014 'get_global()' gains argument 'by_region'
### Update 11.08.2014 fix bug in 'apply_prop_fb()'
### Update 11.08.2014 add 'dots' to 'summary_stats()'
### Update 16.07.2015 change sa_plot colnames
###=========================================================================#

###=========================================================================#
###== FUNCTIONS ============================================================#
###-- generate_samples ........ generate normalized samples
###-- normalize_by_food ....... normalize samples by proportion foodborne
###-- summary_stats ........... summarize samples
###-- sa_plot ................. generate source attribution line plot
###-- sa_barplot .............. generate source attribution barplot
###-- sa_report ............... generate source attribution report
###-- apply_prop_fb ........... apply proportion foodborne
###-- get_global .............. get global results

##-------------------------------------------------------------------------##
## GENERATE NORMALIZED SAMPLES --------------------------------------------##

generate_samples <-
function(n, file, cols, names) {
  ## read .dis file
  x <-
    read.table(
      paste0("dis/", file),
      header = T, sep = "\t", skip = 1)

  ## save percentiles column
  p <- x[, 1]

  ## extract relevant columns
  x <- x[, -1]
  x <- x[, cols]
  names(x) <- names

  ## constrain quantiles
  x[x < 0] <- 0
  x[x > 1] <- 1

  ## sample n values from each column
  samples <-
    apply(x, 2,
          function(x) {
            approx(x = c(0, p, 100),
                   y = c(0, x, 1),
                   xout = runif(n, 0, 100),
                   rule = 1)$y
          })

  ## normalize sampled values row-wise
  samples_norm <-
    t(apply(samples, 1, function(x) x / sum(x)))

  ## set NaN values to zero (obtained when rowsum == 0)
  samples_norm[is.nan(samples_norm)] <- 0

  ## return normalized samples
  return(samples_norm)
}


##-------------------------------------------------------------------------##
## NORMALIZE SAMPLES BY PROPORTION FOODBORNE ------------------------------##

normalize_by_food <-
function(food, exp) {
  ## normalize by proportion foodborne
  food_norm <- food
  exp_food <- sapply(exp, function(x) x[, "food"])
  for (r in names(food)) {
    food_norm[[r]] <- food[[r]] * exp_food[, r]
  }

  ## check if normalization worked
  check <- logical(14)
  for (i in seq(14)) {
    check[i] <- 
      isTRUE(
        all.equal(
          rowSums(food_norm[[i]]),
          exp[[i]][, "food"]))
  }
  if (!all(check)) {
    stop("Normalization failed for subregion ",
         paste0(which(!check), collapse = ","))
  }

  ## return normalized samples
  return(food_norm)
}


##-------------------------------------------------------------------------##
## SUMMARIZE SAMPLES ------------------------------------------------------##

summary_stats <-
function(x, dig = 2, ...) {
  return(
    round(c(mean = mean(x, ...),
            median = median(x, ...),
            quantile(x, c(0.025, 0.050, 0.950, 0.975), ...)),
          dig))
}


##-------------------------------------------------------------------------##
## GENERATE LINE PLOT -----------------------------------------------------##

sa_plot <-
function(samples, what = "food", main = NULL) {
  if (class(samples) == "list") {
    out <- sapply(samples, function(x) summary_stats(x[, what]))
    names <- names(samples)

  } else if (class(samples) == "matrix") {
    out <- apply(samples, 2, summary_stats)
    names <- colnames(samples)

  } else {
    stop("'samples' should be either a list or a matrix.")
  }

  names <- factor(names, levels = rev(names))
  df <- cbind(names, as.data.frame(t(out)))
  colnames(df) <- c("names", "mean", "median", "X2.5", "X5", "X95", "X97.5")

  ggplot(df, aes_string(x = "names", y = "median")) +
    geom_linerange(aes_string(ymin = "X2.5", ymax = "X97.5"),
                   size = 1, col = "grey") +
    geom_linerange(aes_string(ymin = "X5", ymax = "X95"), size = 1) +
    geom_point(size = 2.5, colour = "black") +
    coord_flip() +
    scale_y_continuous("Proportion", limits = c(0, 1)) +
    scale_x_discrete("") +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank()) +
    ggtitle(main)
}

##-------------------------------------------------------------------------##
## GENERATE BARPLOT -------------------------------------------------------##

sa_barplot <-
function(samples) {
  ## calculate average proportion per route, for each subregion
  p <- sapply(samples, colMeans)

  ## create ggplot2 data.frame
  df <-
    data.frame(subregion = factor(rep(colnames(p), each = nrow(p)),
                                  rev(colnames(p))),
               exposure = factor(rep(rownames(p), ncol(p)),
                                 rownames(p)),
               p = c(p))

  ## defined colors
  col <- brewer.pal(nrow(p), "Spectral")

  ## create ggplot
  ggplot(df, aes_string(x = "subregion", y = "p", fill = "exposure")) +
    geom_histogram(stat = "identity", position = "stack") +
    scale_fill_manual(values = col) +
    scale_y_continuous("Proportion") +
    scale_x_discrete("") +
    coord_flip() +
    theme_bw()
}


##-------------------------------------------------------------------------##
## GENERATE PDF REPORT ----------------------------------------------------##

sa_report <-
function(hazard, exp_samples, food_samples, exp_source, food_source) {
  knit2pdf(system.file("sa_report.Rnw", package = "FERG"))
  file.rename("sa_report.pdf", paste0(hazard, "-sa_report.pdf"))
  unlink("figure", recursive = T)
  file.remove(paste0("sa_report.",
                     c("tex", "toc", "aux", "log", "out")))
}


##-------------------------------------------------------------------------##
## APPLY PROPORTION FOODBORNE ---------------------------------------------##

apply_prop_fb <-
function(agent, date_DALY, exp) {
  ## define today
  today <- format(Sys.time(), "%Y%m%d")

  ## define contribution
  file <-
    paste0("DALY_", agent, "_", crpop$ISO3[1], "-", date_DALY, ".RData")
  load(file)
  contrib <- get_contrib(list(DALYrun))

  ## define progress bar
  pb <- txtProgressBar(max = NumCountries, style = 3)

  for (i in seq(NumCountries)) {
    ## load 'DALYrun' & 'incidence' per country
    file <-
      paste0("DALY_", agent, "_", crpop$ISO3[i], "-", date_DALY, ".RData")
    load(file)

    ## which subregion?
    r <- as.character(countryRegion$WHOSubregion[i])

    ## define vector of fb proportion
    fb <- exp[[r]][, "food"]

    ## MULTIPLY 'DALYrun' WITH FOODBORNE PROPORTION
    for (j in seq_along(contrib)) {
      if (contrib[j] == "yld") {
        DALYrun@samples[[j]]$inc <-
          DALYrun@samples[[j]]$inc * fb

        DALYrun@samples[[j]]$cases <-
          DALYrun@samples[[j]]$cases * fb

        DALYrun@samples[[j]]$yld <-
          DALYrun@samples[[j]]$yld * fb

      } else if (contrib[j] == "yll") {
        DALYrun@samples[[j]]$inc <-
          DALYrun@samples[[j]]$inc * fb

        DALYrun@samples[[j]]$deaths <-
          DALYrun@samples[[j]]$deaths * fb

        DALYrun@samples[[j]]$yll <-
          DALYrun@samples[[j]]$yll * fb
      }
    }

    ## MULTIPLY 'incidence' WITH FOODBORNE PROPORTION
    incidence$inc <- incidence$inc * fb
    incidence$cases <- incidence$cases * fb

    ## save DALYrun object
    file_new <-
      paste0("FB_DALY_", agent, "_", crpop$ISO3[i], "-", today, ".RData")
    save(DALYrun, incidence, file = paste0("fb/", file_new))

    ## set progress bar
    setTxtProgressBar(pb, i)
  }

  cat("\n")
}


##-------------------------------------------------------------------------##
## GET GLOBAL RESULTS -----------------------------------------------------##

get_global <-
function(x, what, by_region = FALSE,
         denom = 1e5, rate = FALSE, names = "all", age = "all") {
  ## derive regions
  if (length(x[[1]]) == 6) {
    regions <- countryRegion$WHORegion
  } else if (length(x[[1]]) == 14) {
    regions <- countryRegion$WHOSubregion
  } else {
    stop("regions not recognized.")
  }

  ## derive contributions per node
  contrib <- get_contrib(x[[1]])

  ## check 'names'
  if (length(names) == 1 && names == "all") {
    names <- names(contrib)

  } else {
    for (i in seq_along(names))
      names[i] <- match.arg(names[i], names(contrib))
  }

  ## prepare objects
  global <- global_inc <- 0
  what_yld <- ifelse(what %in% c("daly", "daly_case"), "yld", what)
  what_yll <- ifelse(what %in% c("daly", "daly_case"), "yll", what)

  ## define age group columns
  cols <-
    switch(age,
           "all"   = 1:22,
           "<5"    = c( 1: 2, 12:13),
           "5+"    = c( 3:11, 14:22),
           "5-14"  = c( 3,    14),
           "15-54" = c( 4: 7, 15:18),
           "55+"   = c( 8:11, 19:22))

  ## aggregate results over regions
  region_out <- NULL
  for (r in seq_along(unique(regions))) {

    if (what %in% c("incidence", "daly_case")) {
      region_inc <- rowSums(x[[2]][[r]][, cols])
      global_inc <- region_inc + global_inc
    }

    if (what != "incidence") {
      region <- 0
      for (j in seq_along(contrib)) {
        if (contrib[j] == "yld" &
            what %in% c("cases", "yld", "daly", "daly_case") &
            names(contrib[j]) %in% names) {
          region <- region +
                    rowSums(x[[1]][[r]]@samples[[j]][[what_yld]][, cols])
        } else if (contrib[j] == "yll" &
                   what %in% c("deaths", "yll", "daly", "daly_case") &
                   names(contrib[j]) %in% names) {
          region <- region +
                    rowSums(x[[1]][[r]]@samples[[j]][[what_yll]][, cols])
        }
      }
      global <- global + region

      if (what == "daly_case") {
        region_daly_case <- region / (region_inc / denom)
      }
    }

    if (what == "incidence") {
      region <- region_inc
    } else if (what == "daly_case") {
      region <- region_daly_case
    }

    region_out <- cbind(region_out, region / denom)
  }

  if (what == "incidence") {
    global <- global_inc
  } else if (what == "daly_case") {
    global <- global / (global_inc / denom)
  }

  if (by_region) {
    return(region_out)
  } else {
    return(global / denom)
  }
}