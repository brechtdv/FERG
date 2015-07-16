###=========================================================================#
### Functions to support Calculation Tool Activity
###
### Start  22.11.2013 Brecht & Scott
### Update 27.11.2013 merge_nodes, rpert
### Update 30.11.2013 pop_by_region
### Update 08.12.2013 'sim' gains argument 'type'
### Update 08.12.2013 <DM> gains column 'LOCAL'
###                   |-> 'read_XL' understands 'n_row'
###                   |-> 'get_DALY_model' understands 'TYPE'
### Update 08.12.2013 'stuff_DALY_model' understands dur == 'Lifelong'
### Update 08.12.2013 'pre_sample' understands 'LOCAL'
### Update 08.12.2013 'stratify' understands 'LOCAL'
### Update 23.12.2013 introduce 'sim_Beta' and 'sim_Gamma'
### Update 23.12.2013 fix bug in 'sim_percentiles'
### Update 26.12.2013 fix bug in 'stratify' (target, target_pop)
### Update 05.01.2014 'getDALY_per_country' gains argument 'n_samples'
### Update 05.01.2014 'merge_nodes' derives 'n_samples'
### Update 05.01.2014 <DM> gains column 'INCIDENCE'
###                   |-> 'get_DALY_model' sees 6 columns
###                   |-> 'stratify' and 'merge_nodes' use all nodes
###                   |-> new function 'get_incidence'
###                   |-> new function 'sum_list'
### Update 05.01.2014 fix bug in 'stuff_DALY_model' (dur: DB[[node]])
### Update 14.02.2014 add lognormal distribution to 'sim_percentiles'
### Update 14.02.2014 new function 'sim_lognorm'
### Update 14.02.2014 new function 'optim_lognorm'
### Update 17.02.2014 'split_age' uses country specific age dist
### Update 04.06.2014 'stratify' now understands 10 Groups
### Update 06.06.2014 'pre_sample' converts 'pars' to numeric
### Update 06.06.2014 fix bug in 'stratify'
### Update 06.06.2014 fix bug in 'get_DALY_model' when strat needed
### Update 13.07.2014 allow user to select LE table
### Update 13.07.2014 'getDALY_per_country' measures duration
### Update 28.07.2014 'stuff_DALY_model' uses local_LE when dur == 'Lifelong'
###=========================================================================#

###=========================================================================#
###== FUNCTIONS ============================================================#
###-- readDatabase .......... read entire excel workbook using XLConnect
###-- locY, locX ............ extract row/column indices
###-- read_XL ............... read excel and transform to 'DB' list
###-- pre_sample ............ pre-sampler for INC and PROB nodes
###---| sim ................. generic sampling function
###-----| sim_percentiles ... distribution specific sampling function
###-----| sim_mean .......... distribution specific sampling function
###-----| sim_mode .......... distribution specific sampling function
###-----| sim_gamma ......... distribution specific sampling function
###-----| sim_beta .......... distribution specific sampling function
###-------| sim_fixed ....... actual sampler, fixed
###-------| sim_unif ........ actual sampler, uniform
###-------| sim_pert ........ actual sampler, beta-pert
###-------| sim_Beta ........ actual sampler, beta
###-------| sim_Gamma ....... actual sampler, gamma
###-------| sim_lognorm ..... actual sampler, log-normal
###-------| optim_unif ...... find uniform pars through optimization
###-------| optim_gamma ..... find gamma pars through optimization
###-------| optim_beta ...... find beta pars through optimization
###-- multiply .............. combine INC and PROB to get INC per outcome
###-- stratify .............. stratify INC according to age and sex
###---| split_age ........... split age groups according to FERG groups
###-- merge_nodes ........... merge nodes belonging to same outcome
###-- save_as_DALY .......... save samples as DALY S4 object
###---| get_DALY_model ...... translate 'DM' into DALY model
###---| get_pop_size ........ extract population sizes for current country
###---| stuff_DALY_model .... add data to 'DALYmodel'
###-----| as_num ............ convert character to numeric
###-- getDALY_per_country ... calculate DALYs for each country
###-- get_incidence ......... calculate incidence and incident cases
###---| sum_list ............ sum all elements of a list


##--------------------------------------------------------------------------#
## Read in entire Excel workbook using XLConnect library -------------------#

readDatabase <- 
function(file = NULL) {
  ## evaluate 'file': should end with '.xls' or 'xlsx'
  if (!grepl(".xls$|.xlsx$", file, ignore.case = TRUE))
     stop("The file you selected is not an '.xls[x]' file", call. = FALSE)

  wb <- loadWorkbook(file)
  print(summary(wb))
  return(wb)
}


##--------------------------------------------------------------------------#
## Extract row/column indices ----------------------------------------------#

## extract row index from single cell reference (eg. AF6 => 6)
locY <-
function(cell) {
  return(cref2idx(cell)[1])
}

## extract column index from single cell reference (eg. AF6 => 32)
locX <-
function(cell) {
  return(cref2idx(cell)[2])
}


##--------------------------------------------------------------------------#
## Read Excel file for a given agent and transform to 'DB' list ------------#

read_XL <-
function(XL_wb, DM, KEY) {
  ## settings
  n_nodes <- nrow(DM)
  DB <- vector("list", n_nodes)

  ## possible distributions for excel data
  dists <- c("Percentiles", "Mean & 95%CI", "Min/Mode/Max", "Gamma", "Beta")

  for (node in seq(n_nodes)) {
    ## node info
    node_info <- DM[node, ]

    ## data (inc, prob)
    data <- list()
    data$dist <- KEY[5, node + 1]

    Y <- locY(KEY[1, node + 1])
    X <- locX(KEY[1, node + 1])
    n_row <- ifelse(DM[node, "LOCAL"] == "TRUE", NumCountries, 1)
    n_col <- c(9, 3, 3, 2, 2)[which(data$dist == dists)]

    data$data <-
      readWorksheet(XL_wb, sheet = node + 2, header = FALSE,
                    autofitRow = FALSE, autofitCol = FALSE,
                    startRow = Y, startCol = X,
                    endRow = Y + n_row - 1, endCol = X + n_col - 1)

    if (nrow(data$data) != n_row)
      stop("Something went wrong when reading 'data' @ node ", node, ".")
    if (ncol(data$data) != n_col)
      stop("Something went wrong when reading 'data' @ node ", node, ".")

    data$groups <-
      readWorksheet(XL_wb, sheet = node + 2, header = FALSE,
                    autofitRow = FALSE, autofitCol = FALSE,
                    startRow = Y, startCol = X + 9,
                    endRow = Y + n_row - 1, endCol = X + 11)
    colnames(data$groups) <- c("age", "sex", "duration")

    if (nrow(data$groups) != n_row)
      stop("Something went wrong when reading 'groups' @ node ", node, ".")
    if (ncol(data$groups) != 3)
      stop("Something went wrong when reading 'groups' @ node ", node, ".")

    ## age distribution
    n_age <- length(unique(data$groups$age))
    age <- vector("list", n_age)

    for (j in seq(n_age)) {
      age[[j]]$dist <- KEY[6, node + 1]

      Y <- locY(KEY[2, node + 1]) + (j-1) * 4
      X <- locX(KEY[2, node + 1])

      age[[j]]$grps <-
        readWorksheet(XL_wb, sheet = node + 2, header = FALSE,
                      startRow = Y, startCol = X,
                      endRow = Y, endCol = X + 20)
      age[[j]]$data <-
        readWorksheet(XL_wb, sheet = node + 2, header = FALSE,
                      startRow = Y + 1, startCol = X,
                      endRow = Y + 1, endCol = X + 20)
    }

    ## sex distribution
    n_sex <- length(unique(data$groups$sex))
    sex <- vector("list", n_sex)

    for (j in seq(n_sex)) {
      sex[[j]]$dist <- KEY[7, node + 1]

      Y <- locY(KEY[3, node + 1]) + (j-1) * 4
      X <- locX(KEY[3, node + 1])
      n_row <- ifelse(sex[[j]]$dist == "Mean", 0, 2)

      sex[[j]]$data <-
        readWorksheet(XL_wb, sheet = node + 2, header = FALSE,
                      startRow = Y, startCol = X,
                      endRow = Y + n_row, endCol = X)
    }

    ## duration distribution
    n_dur <- length(unique(data$groups$dur))
    dur <- vector("list", n_dur)

    for (j in seq(n_dur)) {
      dur[[j]]$dist <- KEY[8, node + 1]
      dur[[j]]$strt <- KEY[9, node + 1]

      Y <- locY(KEY[4, node + 1]) + (j-1) * 12
      X <- locX(KEY[4, node + 1])

      dur[[j]]$data <-
        readWorksheet(XL_wb, sheet = node + 2, header = FALSE,
                      startRow = Y, startCol = X,
                      endRow = Y + 10, endCol = X + 5)
    }

    ## disability weight
    dsw <- ifelse(KEY[10, node + 1] == "NA", NA, KEY[10, node + 1])

    ## compile all data
    DB[[node]] <-
      list(node = node_info, data = data,
           age = age, sex = sex,
           dur = dur, dsw = dsw)
  }

  return(DB)
}


##--------------------------------------------------------------------------#
## Pre-sampler for INC and PROB nodes --------------------------------------#

pre_sample <-
function(country, DB, DM, n_samples) {
  n_nodes <- nrow(DM)

  samples <- vector("list", n_nodes)

  for (node in seq(n_nodes)) {
    dist <- DB[[node]]$data$dist
    pars <-
      if (unname(unlist(DB[[node]]$node["LOCAL"])) == "TRUE") {
        unname(unlist(DB[[node]]$data$data[country, ]))
      } else {
        unname(unlist(DB[[node]]$data$data))
      }
    pars <- as.numeric(pars)
    type <- unname(unlist(DB[[node]]$node["TYPE"]))
    samples[[node]] <- sim(n_samples, dist, pars, type)
  }

  return(samples)
}


##--------------------------------------------------------------------------#
## Generic sampling function -----------------------------------------------#

sim <-
function(n, dist, pars, type) {
  ## define seed
  ## -> deals with stratification uncertainty
  ## -> makes results reproducible
  set.seed(264)

  ## generate samples
  samples <-
    switch(dist,
           "Percentiles" = sim_percentiles(n, pars, type),
           "Mean & 95%CI" = sim_mean(n, pars, type),
           "Min/Mode/Max" = sim_mode(n, pars, type),
           "Gamma" = sim_gamma(n, pars, type),
           "Beta" = sim_beta(n, pars, type))
  return(samples)
}


##--------------------------------------------------------------------------#
## Distribution specific sampling functions --------------------------------#
## -> check whether input is OK, and pass on to actual sampler -------------#

sim_percentiles <-
function(n, pars, type) {
  # all nine filled in -> CUMULATIVE
  if (all(!is.na(pars[1:9]))) {
    samples <- sim_cum(n, pars)

  # only 2.5%, 5%, 25%, 50%, 75%, 95%, 97.5% filled in -> LOG-NORMAL
  } else if (all(!is.na(pars[c(2:8)])) &
             all(is.na(pars[c(1, 9)]))) {
    pars_lognorm <- optim_lognorm(pars[c(2:8)])
    samples <- sim_lognorm(n, pars_lognorm)

  # only 0%, 50%, 100% filled in -> PERT
  } else if (all(!is.na(pars[c(1, 5, 9)])) &
             all(is.na(pars[c(2:4, 6:8)]))) {
    samples <- sim_pert(n, pars[c(1, 5, 9)])

  # only 2.5%, 50%, 97.5% filled in -> GAMMA/BETA
  } else if (all(!is.na(pars[c(2, 5, 8)])) &
             all(is.na(pars[c(1, 3:4, 6:7, 9)]))) {
    if (type == "INC") {
      pars_gamma <- optim_gamma(pars[c(2, 5, 8)], 0.95)
      samples <- sim_Gamma(n, pars_gamma)
    } else if (type == "PROB") {
      pars_beta <- optim_beta(pars[c(2, 5, 8)], 0.95)
      samples <- sim_Beta(n, pars_beta)
    } else {
      stop("sim_percentiles() could not find a proper type.")
    }

  # only 5%, 50%, 95% filled in -> GAMMA/BETA
  } else if (all(!is.na(pars[c(3, 5, 7)])) &
             all(is.na(pars[c(1:2, 4, 6, 8:9)]))) {
    if (type == "INC") {
      pars_gamma <- optim_gamma(pars[c(3, 5, 7)], 0.90)
      samples <- sim_Gamma(n, pars_gamma)
    } else if (type == "PROB") {
      pars_beta <- optim_beta(pars[c(3, 5, 7)], 0.90)
      samples <- sim_Beta(n, pars_beta)
    } else {
      stop("sim_percentiles() could not find a proper type.")
    }

  # only 25%, 50%, 75% filled in -> GAMMA/BETA
  } else if (all(!is.na(pars[c(4, 5, 6)])) &
             all(is.na(pars[c(1:3, 7:9)]))) {
    if (type == "INC") {
      pars_gamma <- optim_gamma(pars[c(4, 5, 6)], 0.50)
      samples <- sim_Gamma(n, pars_gamma)
    } else if (type == "PROB") {
      pars_beta <- optim_beta(pars[c(4, 5, 6)], 0.50)
      samples <- sim_Beta(n, pars_beta)
    } else {
      stop("sim_percentiles() could not find a proper type.")
    }

  # only 0% and 100% filled in -> UNIFORM
  } else if (all(!is.na(pars[c(1, 9)]))&
             all(is.na(pars[c(2:8)]))) {
    samples <- sim_unif(n, pars[c(1, 9)])

  # only 2.5% and 97.5% filled in -> UNIFORM
  } else if (all(!is.na(pars[c(2, 8)]))&
             all(is.na(pars[c(1, 3:7, 9)]))) {
    pars_unif <- optim_unif(pars[c(2, 8)], 0.95)
    samples <- sim_unif(n, pars_unif)

  # only 5% and 95% filled in -> UNIFORM
  } else if (all(!is.na(pars[c(3, 7)]))&
             all(is.na(pars[c(1:2, 4:6, 8:9)]))) {
    pars_unif <- optim_unif(pars[c(3, 7)], 0.90)
    samples <- sim_unif(n, pars_unif)

  # only 25% and 75% filled in -> UNIFORM
  } else if (all(!is.na(pars[c(4, 6)]))&
             all(is.na(pars[c(1:3, 5, 7:9)]))) {
    pars_unif <- optim_unif(pars[c(4, 6)], 0.50)
    samples <- sim_unif(n, pars_unif)

  # only 50% filled in -> FIXED
  } else if (!is.na(pars[5])&
             all(is.na(pars[c(1:4, 6:9)]))) {
    samples <- sim_fixed(n, pars[5])

  # something else -> error!
  } else {
    stop("sim_percentiles() could not find a proper distribution.")
  }
}

sim_mean <-
function(n, pars, type) {
  # all three filled in -> GAMMA/BETA
  if (all(!is.na(pars))) {
    if (type == "INC") {
      pars_gamma <- optim_gamma(pars, 0.95)
      samples <- sim_Gamma(n, pars_gamma)
    } else if (type == "PROB") {
      pars_beta <- optim_beta(pars, 0.95)
      samples <- sim_Beta(n, pars_beta)
    } else {
      stop("sim_mean() could not find a proper type.")
    }

  # only min and max filled in -> UNIFORM
  } else if (all(!is.na(pars[c(1, 3)])) & is.na(pars[2])) {
    samples <- sim_unif(n, pars[c(1, 3)])

  # only mode filled in -> FIXED
  } else if (all(is.na(pars[c(1, 3)])) & !is.na(pars[2])) {
    samples <- sim_fixed(n, pars[2])

  # something else -> error!
  } else {
    stop("sim_mean() could not find a proper distribution.")
  }
}

sim_mode <-
function(n, pars, type) {
  # all three filled in -> PERT
  if (all(!is.na(pars))) {
    samples <- sim_pert(n, pars)

  # only min and max filled in -> UNIFORM
  } else if (all(!is.na(pars[c(1, 3)])) & is.na(pars[2])) {
    samples <- sim_unif(n, pars[c(1, 3)])

  # only mode filled in -> FIXED
  } else if (all(is.na(pars[c(1, 3)])) & !is.na(pars[2])) {
    samples <- sim_fixed(n, pars[2])

  # something else -> error!
  } else {
    stop("sim_mode() could not find a proper distribution.")
  }
}

sim_gamma <-
function(n, pars, type) {
  # error message if type = 'PROB'
  if (type == "PROB")
    stop("sim_gamma() is not appropriate for PROB nodes.")

  # all two filled in -> GAMMA
  if (all(!is.na(pars))) {
    samples <- sim_Gamma(n, pars)

  # something else -> error!
  } else {
    stop("sim_gamma() could not find a proper distribution.")
  }
}

sim_beta <-
function(n, pars, type) {
  # error message if type = 'INC'
  if (type == "INC")
    stop("sim_beta() is not appropriate for INC nodes.")

  # all two filled in -> BETA
  if (all(!is.na(pars))) {
    samples <- sim_Beta(n, pars)

  # something else -> error!
  } else {
    stop("sim_beta() could not find a proper distribution.")
  }
}


##--------------------------------------------------------------------------#
## Actual samplers ---------------------------------------------------------#

sim_fixed <-
function(n, par) {
  return(rep(par, n))
}

sim_pert <-
function(n, par) {
  bp <- betaPERT(a = par[1], m = par[2], b = par[3])
  return(rbeta(n, bp$alpha, bp$beta) * (bp$b - bp$a) + bp$a)
}

sim_unif <-
function(n, par) {
  return(runif(n, par[1], par[2]))
}

sim_Beta <-
function(n, par) {
  return(rbeta(n, par[1], par[2]))
}

sim_Gamma <-
function(n, par) {
  return(rgamma(n, par[1], par[2]))
}

sim_lognorm <-
function(n, par) {
  return(rlnorm(n, par[1], par[2]))
}

sim_cum <-
function(n, par) {
  stop("'sim_cum()' not yet implemented.")
}

##--------------------------------------------------------------------------#
## Find Gamma & Beta pars through optimization -----------------------------#

optim_unif <-
function(par, p) {
  margin <- ((diff(par) / p) - diff(par)) / 2
  return(c(par[1] - margin, par[2] + margin))
}

optim_gamma <-
function(par, p) {
  target <- par[c(1, 3)]
  p <- c(0, p) + (1 - p)/2

  f <-
    function(x, mean, p, target) {
      dev <- qgamma(p = p, shape = x, rate = x / mean)
      return(sum((dev - target) ^ 2))
    }

  shape <-
    optimize(f, c(0, 1000), mean = par[2], p = p, target = target)$minimum
  rate <-
    shape / par[2]

  return(c(shape, rate))
}

optim_beta <-
function(par, p) {
  target <- par[c(1, 3)]
  p <- c(0, p) + (1 - p)/2

  f <-
    function(x, mean, p, target) {
      dev <- qbeta(p = p, shape1 = x, shape2 = (x * (1 - mean))/mean)
      return(sum((dev - target) ^ 2))
    }

  shape1 <-
    optimize(f, c(0, 1000), mean = par[2], p = p, target = target)$minimum
  shape2 <-
    (shape1 * (1 - par[2])) / par[2]

  return(c(shape1, shape2))
}

optim_lognorm <-
function(target) {
  f <-
  function(par, target) {
    sum((qlnorm(c(.025, .05, .25, .5, .75, .95, .975),
                par[1], par[2]) - target) ^ 2)
  }

  par <- optim(c(1, 2), f, target = target)$par

  return(par)
}

##--------------------------------------------------------------------------#
## Combine INC and PROB samples to get INC per outcome ---------------------#

multiply <-
function(samples, DM) {
  n_nodes <- nrow(DM)
  INC_samples <- samples

  while(!all(DM[, "PARENT"] == "NA")) {
    for (node in seq(n_nodes)) {
      if (DM[node, "PARENT"] != "NA") {
        parent <- which(DM[, "NODE"] == DM[node, "PARENT"])
        INC_samples[[node]] <- INC_samples[[node]] * INC_samples[[parent]]
        DM[node, "PARENT"] <- DM[parent, "PARENT"]
      }
    }
  }

  return(INC_samples)
}


##--------------------------------------------------------------------------#
## Stratify INC per outcome according to age and sex -----------------------#
## -> need to split age groups according to FERG 11 age groups! ------------#

stratify <-
function(country, INC, DB, DM) {
  n_nodes <- length(INC)
  INC_samples_strat <- vector("list", sum(n_nodes))

  for (node in seq(n_nodes)) {
    if (unname(unlist(DB[[node]]$node["LOCAL"])) == "TRUE") {
      groups <- unlist(DB[[node]]$data$groups[country, ])
    } else {
      groups <- unlist(DB[[node]]$data$groups)
    }
    groups <- match(groups, paste0("G", 1:10))

    ## obtain age distribution
    #age_dist <- DB[[node]]$age[[groups[1]]]$dist
    age_grps <- unlist(DB[[node]]$age[[groups[1]]]$grps)
    age_data <- unlist(DB[[node]]$age[[groups[1]]]$data)
    age_prob <- split_age(age_grps, age_data, country)
    if (!isTRUE(all.equal(sum(age_prob), 1)))
      stop("stratify() found an incorrect age distribution @ node ",
           node, ".")

    ## obtain sex distribution
    #sex_dist <- DB[[node]]$sex[[groups[2]]]$dist
    sex_data <- unlist(DB[[node]]$sex[[groups[2]]]$data)
    sex_prob <- c(sex_data, 1 - sex_data)
    if (!isTRUE(all.equal(sum(sex_prob), 1)))
      stop("stratify() found an incorrect sex distribution @ node ",
           node, ".")

    ## combine age and sex distribution
    age_sex <- c(matrix(age_prob, ncol = 1) %*% sex_prob)

    ## apply age and sex distribution to INC samples
    ## note: target population different for age/sex sheets!
    if (grepl("@", DM[node, "NODE"])) {
      target <- age_sex != 0
    } else {
      target <- TRUE
    }

    target_pop <- numeric(22)
    target_pop[target] <- c(get_pop_size(country)[, -1])[target]

    CASES <- INC[[node]] * sum(target_pop)
    CASES_strat <- matrix(CASES, ncol = 1) %*% age_sex
    INC_strat <- t(t(CASES_strat) / target_pop)

    ## replace NaN by 0
    INC_strat <- replace(INC_strat, is.nan(INC_strat), 0)

    INC_samples_strat[[node]] <- INC_strat
  }

  ## return stratified nodes
  return(INC_samples_strat)
}


##--------------------------------------------------------------------------#
## Split age groups according to FERG 11 age groups ------------------------#

split_age <-
function(grps, data, country) {
  ## see file 'MS_pop_2005-2010.xlsx'
  #world_2005 <-
  #  c(0.01941, 0.07584, 0.18710, 0.18076, 0.15474, 0.13536, 0.10551, 
  #    0.06848, 0.04554, 0.02214, 0.00511)
  #world_2010 <- age_dist <-
  #  c(0.01897, 0.07339, 0.17575, 0.17601, 0.15377, 0.13766, 0.10955, 
  #    0.07896, 0.04613, 0.02368, 0.00612)

  #> or use country-specific dist?? difference is not significant!
  #> length(country) == 194
  age_dist <-
    tapply(apply(Popul[, , country], 1, sum),
           c(1, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10, 10, 11),
           sum) / sum(Popul[, , country])

  ## define vector of 'split' data
  probs <- numeric()

  ## split '<' group
  if (grepl("<", grps[1], fixed = TRUE)) {
    upr <- as.numeric(strsplit(grps[1], "<", fixed = TRUE)[[1]][2])
    which_upr <- which(upr == FERG_agrps) - 1
    split <- age_dist[seq(which_upr)] / sum(age_dist[seq(which_upr)])
    probs[seq(which_upr)] <- data[1] * split
  }

  ## split '-' groups
  for (i in seq_along(grps)) {
    if (grepl("-", grps[i], fixed = TRUE)) {
      lwr <- as.numeric(strsplit(grps[i], "-", fixed = TRUE)[[1]][1])
      upr <- as.numeric(strsplit(grps[i], "-", fixed = TRUE)[[1]][2]) + 1
      which_lwr <- which(lwr == FERG_agrps)
      which_upr <- which(upr == FERG_agrps) - 1
      split <- age_dist[seq(which_lwr, which_upr)] /
                 sum(age_dist[seq(which_lwr, which_upr)])
      probs <- c(probs, data[i] * split)
    }
  }

  ## split '+' group
  if (grepl("+", tail(grps, 1L), fixed = TRUE)) {
    lwr <- as.numeric(strsplit(tail(grps, 1L), "+", fixed = TRUE)[[1]])
    which_lwr <- which(lwr == FERG_agrps)
    split <- age_dist[seq(which_lwr, 11L)] /
               sum(age_dist[seq(which_lwr, 11L)])
    probs <- c(probs, tail(data, 1L) * split)
  }

  ## return FERG-compatable probabilities
  return(probs)
}


##--------------------------------------------------------------------------#
## Merge nodes belonging to same outcome -----------------------------------#

merge_nodes <-
function(INC, DM) {
  n_nodes <- length(INC)
  n_samples <- nrow(INC[[1]])
  strat <- grepl("@", DM[, "NODE"])

  node_names <-
    unique(sapply(DM[strat, "NODE"],
                  function(x) strsplit(x, "@")[[1]][1]))
  merge_INC <- list()
  for (node in seq_along(node_names)) {
    merge_INC[[node_names[node]]] <- matrix(0, ncol = 22, nrow = n_samples)
  }

  for (node in seq(n_nodes)[strat]) {
    node_name <- strsplit(DM[node, "NODE"], "@")[[1]][1]
    merge_INC[[node_name]] <- merge_INC[[node_name]] + INC[[node]]
  }

  ## return combination of merged and unstratified nodes
  INC_merge_all <- c(merge_INC, INC[!strat])
  names(INC_merge_all) <- c(names(merge_INC), DM[!strat, "NODE"])
  return(INC_merge_all)
}


##--------------------------------------------------------------------------#
## Save samples as DALY S4 object ------------------------------------------#

save_as_DALY <-
function(i, INC, DB, DM, agent, LE) {
  ## name of the health cause
  disease <- agent

  ## get population size
  population <- get_pop_size(i)

  ## transform DM to DALY model (be careful with age/sex sheets!)
  model <- get_DALY_model(DM)

  ## compile 'DALYmodel' object
  DALY <- setDALYmodel(disease, population, LE, model)

  ## define list of 'DALYdata' objects
  DALY <- stuff_DALY_model(i, DALY, INC, DB, DM, LE)

  ## return complete 'DALYmodel' object
  return(DALY)
}


##--------------------------------------------------------------------------#
## Translate 'DM' into DALY model  -----------------------------------------#

get_DALY_model <-
function(DM) {
  DM_merged <-
    data.frame("NODE" = character(),
               "PARENT" = character(),
               "TYPE" = character(),
               "LOCAL" = character(),
               "CONTRIBUTION" = character(),
               "INCIDENCE" = character())
  n_contrib <- sum(DM[, "CONTRIBUTION"] != "NA")
  model <- vector("list", n_contrib)

  ## merge nodes (same order as 'INC_strat_samples_merged'!)
  DM_contrib <- DM[DM[, "CONTRIBUTION"] != "NA", ]
  strat <- grepl("@", DM_contrib[, "NODE"])

  if (any(strat)) {
    for (node in seq(nrow(DM_contrib))[strat]) {
      node_name <- strsplit(DM_contrib[node, "NODE"], "@")[[1]][1]
      if (!(node_name %in% DM_merged[, "NODE"]))
        DM_merged <-
          rbind(DM_merged,
                cbind(NODE = node_name, DM_contrib[node, 2:6]))
    }
    DM_merged <- rbind(DM_merged, DM_contrib[!strat, ])

  } else {
    DM_merged <- DM_contrib
  }

  ## identify nodes with contribution
  contrib <- DM_merged[, "CONTRIBUTION"] != "NA"

  ## define node names
  names(model) <- DM_merged[, "NODE"]

  ## add parent, type, contribution
  for (node in seq(n_contrib)) {
    model[[node]] <-
      c(ifelse(DM_merged[node, 2] == "NA", NA, DM_merged[node, 2]),
        "inc",
        ifelse(DM_merged[node, 5] == "NA", NA, tolower(DM_merged[node, 5])))
  }

  ## return disease model, with only the contributing & merged nodes
  return(model[seq_along(contrib)])
}


##--------------------------------------------------------------------------#
## Extract population sizes for current country  ---------------------------#

get_pop_size <-
function(country) {
  ## find country index if country name specified
  if (class(country) == "character")
    country <- which(crpop$Country == country)

  ## get full population table
  pop_full <- Popul[, , country]

  ## merge FERG bounds; collapsed male pop; collapsed female pop
  pop <-
    cbind(FERG_agrps,
          tapply(pop_full[, 1], c(1:2, rep(3:10, each = 2), 11), sum),
          tapply(pop_full[, 2], c(1:2, rep(3:10, each = 2), 11), sum))

  ## return merged population size
  return(unname(pop))
}


##--------------------------------------------------------------------------#
## Add data to 'DALYmodel' -------------------------------------------------#

stuff_DALY_model <-
function(country, DALY, INC, DB, DM, LE) {
  ## calculate number of nodes
  n_nodes <- length(DALY@data)

  ## set data per node
  for (node in seq(n_nodes)) {
    ## which 'DB' node corresponds to current 'DALY' node?
    DB_node <- which(names(INC) == names(DALY@model)[node])

    ## stuff DALY if node contributes YLDs
    if (DALY@model[[node]][3] == "yld") {

      # incidence
      DALY@data[[node]]$inc <- INC[[DB_node]]

      # disability weight
      DALY@data[[node]]$dsw <-
        setDALYdata("fixed", "none", as_num(DB[[DB_node]]$dsw))

      # onset
      DALY@data[[node]]$ons <-
        setDALYdata("fixed", "age", FERG_means)

      # duration
      if (unname(unlist(DB[[DB_node]]$node["LOCAL"])) == "TRUE") {
        group <- unlist(DB[[DB_node]]$data$groups[country, 3])
      } else {
        group <- unlist(DB[[DB_node]]$data$groups[1, 3])
      }
      group <- match(group, c("G1", "G2", "G3", "G4", "G5"))
      dur_dist_DB <- DB[[DB_node]]$dur[[group]]$dist
      if (dur_dist_DB == "Mean") {
        dur_dist <- "fixed"
      } else if (dur_dist_DB == "Min/Mode/Max") {
        dur_dist <- "pert"
      } else {
        stop("duration distribution '", dur_dist_DB, "' not yet implemented.")
      }

      dur_pars <- as.matrix(DB[[DB_node]]$dur[[group]]$data) ##> check!
      dur_strt <- DB[[DB_node]]$dur[[group]]$strt

      if (length(DB[[DB_node]]$dur[[group]]$data) != 0 &&
          any(DB[[DB_node]]$dur[[group]]$data == "Lifelong")) {
        if (dur_dist != "fixed")
          stop("'Lifelong' can only occur in 'fixed' distribution")
        dur_pars <-
          c(getLE(x = matrix(rep(FERG_means, 2), nrow = 1),
                  LE = local_LE[[country]],
                  n_agegroups = length(FERG_means)))
        dur_strt <- "full"
      }

      DALY@data[[node]]$dur <-
        setDALYdata(dur_dist, dur_strt, dur_pars)  

    ## stuff DALY if node contributes YLLs
    } else if (DALY@model[[node]][3] == "yll") {

      # incidence
      DALY@data[[node]]$inc <- INC[[DB_node]]

      # age at death
      DALY@data[[node]]$aad  <- setDALYdata("fixed", "age", FERG_means)

    ## if something went wrong..
    } else {
      stop("stuff_DALY_model() found an error.")
    }
  }

  ## return stuffed DALY
  return(DALY)
}


##--------------------------------------------------------------------------#
## Convert character to numeric --------------------------------------------#

as_num <-
function(x) {
  as.numeric(gsub(",", ".", x))
}


##--------------------------------------------------------------------------#
## Perform DALY calculations per country -----------------------------------#

getDALY_per_country <-
function(agent, DB, DM, today,
         LE = c("WHO", "GBD", "CD"), n_samples = 1000) {
  ## check 'LE'
  LE <- match.arg(LE)
  std_LE <-
    switch(LE,
           WHO = std_LE_WHO,
           GBD = std_LE_GBD,
           CD = std_LE_CD)

  ## define start time
  start_time <- Sys.time()

  ## define progress bar
  pb <- txtProgressBar(max = NumCountries, style = 3)

  for (i in seq(NumCountries)) {
    ## sample from INC and PROB distributions
    samples <- pre_sample(i, DB, DM, n_samples)

    ## combine INC and PROB samples to get INC per node
    INC_samples <- multiply(samples, DM)

    ## stratify INC per node according to age and sex
    INC_strat_samples <- stratify(i, INC_samples, DB, DM)

    ## merge nodes belonging to same outcome
    INC_strat_samples_merged <- merge_nodes(INC_strat_samples, DM)

    ## calculate incidence
    incidence <- get_incidence(i, INC_strat_samples_merged, DM)

    ## create 'DALYmodel' object
    DALYmodel <-
      save_as_DALY(i, INC_strat_samples_merged, DB, DM, agent, std_LE)

    ## calculate DALYs
    DALYrun <- getDALY(DALYmodel, n_samples, 0, 0)

    ## save DALYs
    save(DALYrun, incidence,
         file = paste0("DALY_", agent, "_", crpop$ISO3[i], "-", today,
                       ".RData"))

    ## set progress bar
    setTxtProgressBar(pb, i)
  }

  cat("\n")

  ## show elapsed time
  diff <- Sys.time() - start_time
  cat("DALY calculation took", round(diff), attr(diff, "units"), "\n\n")
}


##--------------------------------------------------------------------------#
## Calculate incidence and incident cases ----------------------------------#

get_incidence <-
function(country, INC, DM) {
  ## derive number of samples
  n_samples <- nrow(INC[[1]])

  ## define incidence contributors
  is_inc_node <- DM[, "INCIDENCE"] == TRUE
  INC_inc_nodes <- INC[is_inc_node]

  ## calculate incident cases
  pop <- c(get_pop_size(country)[, -1])
  INC_inc <- sum_list(INC_inc_nodes)
  CASES_inc <- t(t(INC_inc) * pop)

  ## return incidence rate and incident cases
  return(list(inc = INC_inc, cases = CASES_inc))
}


##--------------------------------------------------------------------------#
## Sum all elements of a list ----------------------------------------------#

sum_list <-
function (list) {
  n_nodes <- length(list)
  sum <-
    eval(parse(text = paste("list[[", seq(n_nodes), "]]",
                            collapse = " + ")))
  return(sum)
}
