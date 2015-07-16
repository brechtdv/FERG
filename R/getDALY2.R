###=========================================================================#
### DALY CALCULATOR functions
### First  22.11.2013 Brecht
### Update 11.12.2013 analytical solution of DALY integral
### Update 23.12.2013 define C and beta in 'burden'
### Update 28.12.2013 dur -> min(dur, LE)
###=========================================================================#

###=========================================================================#
###== CLASSES ==============================================================#
###-- DALYmodel ................ DALY calculation model
###-- DALYdata ................. DALY input data (dist, strat, pars)
###-- DALYrun .................. DALY simulation output

###=========================================================================#
###== FUNCTIONS ============================================================#
###-- setDALYdata .............. create new 'DALYdata' object
###-- setDALYmodel ............. create new 'DALYmodel' object
###-- burden ................... calculate YLD/YLL
###-- getLE .................... linear imputation of life expectancy
###-- getDALY .................. calculate DALYs from 'DALYmodel' object
###-- rpert .................... generate random Beta-PERT deviates
###-- getSamples ............... generate random deviates from 'DALYdata' object
###-- getParam ................. format parameters array


##-------------------------------------------------------------------------##
## S4 CLASSES -------------------------------------------------------------##

## Define S4 classes
setClass("DALYmodel",
         representation(disease = "character",
                        population = "matrix",
                        LE_table = "matrix",
                        model = "list",
                        data = "list"))

setClassUnion("par", c("numeric", "matrix", "array"))

setClass("DALYdata",
         representation(dist = "character",
                        strat = "character",
                        param = "par"))

setClass("DALYrun",
         contains = "DALYmodel",
         representation(samples = "list",
                        K = "numeric",
                        r = "numeric"))

## Create a new 'DALYdata' object
setDALYdata <-
function(dist, strat, param){
  new("DALYdata", dist = dist, strat = strat, param = param)
}

## Create a new 'DALYmodel' object
setDALYmodel <-
function(disease, population, LE_table, model){
  ## number of nodes
  n_nodes <- length(model)

  ## make skeleton DALYmodel
  x <- new("DALYmodel",
           disease = disease,
           population = population,
           LE_table = LE_table,
           model = model,
           data = vector("list", n_nodes))
  names(x@data) <- names(model)

  ## add skeleton DALYdata objects
  for (i in seq(n_nodes)){
    x@data[[i]] <- list()
    table <- model[[i]][2]
    x@data[[i]][[table]] <- new("DALYdata")
    if(!is.na(model[[i]][3])){
      if(model[[i]][3] == "yld"){
        x@data[[i]]$dsw <- new("DALYdata")
        x@data[[i]]$ons <- new("DALYdata")
        x@data[[i]]$dur <- new("DALYdata")
      } else {
        x@data[[i]]$aad <- new("DALYdata")
        x@data[[i]]$lfe <- new("matrix")
      }
    }
  }

  ## return DALYmodel
  return(x) 
}

## functions ===============================================================#

burden <-
function(N, DW, A, L, K, r){
  ## N = number of cases
  ## DW = disability weight
  ## A = onset
  ## L = duration
  ## K = age weighting constant (0,1)
  ## r = time discount rate
     C    <- 0.1658
     beta <- 0.04

  ## burden calculation formula
  ## formula from WHO 'DALY calculation template'
  burden <-
    if(r == 0){
      N * DW *
      (K * C *
       (exp(-beta * A) / beta^2) *
       (exp(-beta * L) * (-beta * (L + A) - 1) - (-beta * A - 1)) +
       ((1 - K) * L)
      )
    } else {
      N * DW *
        (K *
         ((C * exp(r * A)) / ((-(beta + r))^2)) *
         ((exp((-(beta + r)) * (L + A)) * ((-(beta + r)) * (L + A) - 1)) -
          (exp((-(beta + r)) * A) * ((-(beta + r)) * A - 1))) +
         ((1 - K) / r) * ((1 - exp(-r * L)))
        )
    }

  ## return burden estimates
  return(burden)
}


getLE <-
function(x, LE, n_agegroups){
  LE_m <- approx(LE[, 1], LE[, 2], x[, seq(1, n_agegroups)])$y
  LE_f <- approx(LE[, 1], LE[, 3], x[, seq(n_agegroups + 1, ncol(x))])$y
  LE_m <- matrix(LE_m, ncol = n_agegroups)
  LE_f <- matrix(LE_f, ncol = n_agegroups)    
  return(cbind(LE_m, LE_f))
}

getDALY <-
function(x, n, K, r){
  ##(1) initialize nodes ===================================================
  ## settings
  n_agegroups <- nrow(x@population)
  nodes <- names(x@model)

  ## create empty 'DALYrun' object
  y <- new("DALYrun", K = K, r = r)

  ## stuff 'DALYrun' object per node
  for (i in seq_along(nodes)){
    y@samples[[i]] <- list()
    for (j in seq_along(x@data[[i]]))
      ## if 'DALYdata' -> generate samples
      if (class(x@data[[i]][[j]]) == "DALYdata"){
        y@samples[[i]][[j]] <-
          getSamples(x@data[[i]][[j]],
                     n = n,
                     n_agegroups = n_agegroups)

      ## if 'matrix' -> copy samples
      } else if (class(x@data[[i]][[j]]) == "matrix"){
        y@samples[[i]][[j]] <- x@data[[i]][[j]]
      }

    ## assign data element names
    names(y@samples[[i]]) <- names(x@data[[i]])
  }

  ## assign node names
  names(y@samples) <- names(x@model)

  ##(2) combine nodes according to model ===================================
  for (i in seq_along(nodes)){

    ## calculate 'inc' for all nodes ---------------------------------------
    ## multiply current node 'prob' with parent node 'inc'
    if (x@model[[nodes[[i]]]][2] != "inc"){
      parent <- x@model[[nodes[[i]]]][1]
      inc <- y@samples[[parent]]$inc
      prb <- y@samples[[nodes[[i]]]]$prob
      y@samples[[nodes[[i]]]]$inc <- inc * prb
    }

    ## calculate 'cases' and 'deaths' for all nodes -----------------------
    inc <- y@samples[[nodes[[i]]]]$inc
    pop <- x@population[, -1]
    if (!is.na(x@model[[nodes[[i]]]][3])){
      if (x@model[[nodes[[i]]]][3] == "yld"){
        y@samples[[nodes[[i]]]]$cases <- t(t(inc) * c(pop))
      } else {
        y@samples[[nodes[[i]]]]$deaths <- t(t(inc) * c(pop))
      }
    }

    ## obtain life expectancy for all yll contributors -------------------
    ## = linear interpolation from 'LE_table'
    if (!is.na(x@model[[nodes[[i]]]][3]) &&
          x@model[[nodes[[i]]]][3] == "yll"){
       y@samples[[nodes[[i]]]]$lfe <-
         getLE(x = y@samples[[nodes[[i]]]]$aad,
               LE = x@LE_table,
               n_agegroups = n_agegroups)
    }

    ## set duration to min(dur, LE) for all yld contributors -------------
    if (!is.na(x@model[[nodes[[i]]]][3]) &&
          x@model[[nodes[[i]]]][3] == "yld"){
       y@samples[[nodes[[i]]]]$dur <-
         pmin(y@samples[[nodes[[i]]]]$dur,
              getLE(x = y@samples[[nodes[[i]]]]$ons,
                    LE = x@LE_table,
                    n_agegroups = n_agegroups))
    }

    ## calculate 'yld' and 'yll' for all nodes ---------------------------
    if (!is.na(x@model[[nodes[[i]]]][3])){
      ## calculate YLDs
      if (x@model[[nodes[[i]]]][3] == "yld"){
        y@samples[[nodes[[i]]]]$yld <-
          burden(N = y@samples[[nodes[[i]]]]$cases,
                 DW = y@samples[[nodes[[i]]]]$dsw,
                 A = y@samples[[nodes[[i]]]]$ons,
                 L = y@samples[[nodes[[i]]]]$dur,
                 K = K,
                 r = r)

      ## calculate YLLs
      } else {
        y@samples[[nodes[[i]]]]$yll <-
          burden(N = y@samples[[nodes[[i]]]]$deaths,
                 DW = 1,
                 A = y@samples[[nodes[[i]]]]$aad,
                 L = y@samples[[nodes[[i]]]]$lfe,
                 K = K,
                 r = r)
      }
    }
  }

  ## done!
  return(y)
}

rpert <-
function(n, min, mode, max){
  par <- betaPERT(a = min, m = mode, b = max)
  rbeta(n, par$alpha, par$beta) * (max - min) + min
}

getSamples <-
function(x, n, n_agegroups){
  ## generate random deviates from DALYdata
  ## switch -> select distribution
  ## mapply -> generate for all strata
  samples <- 
    switch(tolower(x@dist),
           fixed = mapply(rep,
                          x = getParam(x, 1),
                          MoreArgs = list(length = n)),
           uniform = mapply(runif,
                            min = getParam(x, 1),
                            max = getParam(x, 2),
                            MoreArgs = list(n = n)),
           pert = mapply(rpert,
                         min = getParam(x, 1),
                         mode = getParam(x, 2),
                         max = getParam(x, 3),
                         MoreArgs = list(n = n)),
           gamma = mapply(rgamma, 
                          shape = getParam(x, 1),
                          rate = getParam(x, 2),
                          MoreArgs = list(n = n)),
           beta = mapply(rbeta,
                         shape1 = getParam(x, 1),
                         shape2 = getParam(x, 2),
                         MoreArgs = list(n = n)))

  ## make a full matrix to allow combinations
  ## note: this might take a lot of memory
  ##       + other options are available if needed
  if (tolower(x@strat) == "sex"){
    samples <- t(apply(samples, 1, rep, each = n_agegroups))
  }
  if (tolower(x@strat) == "age"){
    samples <- t(apply(samples, 1, rep, times = 2))
  }
  if (tolower(x@strat) == "none"){
    samples <- t(apply(samples, 1, rep, times = 2 * n_agegroups))
  }
  return(samples)
}

getParam <-
function(x, par){
  if (tolower(x@dist) == "fixed"){
    param <- c(x@param)
  } else {
    if (tolower(x@strat) == "age" | tolower(x@strat) == "none"){
      param <- c(x@param[, par])
    } else {
      param <- c(x@param[, par, ])
    }
  }
  return(param)
}