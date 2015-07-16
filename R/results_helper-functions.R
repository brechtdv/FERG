## FERG RESULTS :: HELPER FUNCTIONS

## load_results ................ load (fb_)results_per_region/subregion
## get_total ................... get summary per (sub)region for 1 agent
## get_samples ................. get DALY samples by (sub)region for 1 agent
## get_all_samples ............. get DALY samples by (sub)region for >1 agent
## get_samples_per_outcome ..... idem + per outcome
## sum_agents .................. sum samples of different agents
## summarize_all ............... get summary for all agents 
## extract ..................... extract row/col from summaries
## get_local_LE ................ get age-specific LE for a given country

## load (fb_)results_per_region/subregion
load_results <-
function(agent, by, fb) {
  ## return NA if fb not available
  if (fb & length(agent) == 3) return(NA)

  ## should fb results be returned?
  this_fb <- fb & !is.na(agent[4])

  ## define data
  date <- ifelse(this_fb, agent[4], agent[3])
  dir <- agent[2]
  agent <- agent[1]

  ## define prefix
  pre <- ifelse(this_fb, "/FB_DALY_", "/DALY_")

  ## check if '(FB_)DALY_subregion' exists
  subregion_exists <-
    any(grepl(paste0(ifelse(this_fb, "FB_", ""), "DALY_subregion"), dir(dir)))

  ## load required results
  if (by == "subregion") {
    if (subregion_exists) {
      load(paste0(dir, pre, "subregion_", agent, "-", date, ".RData"))

    } else {
      load(paste0(dir, pre, "region_", agent, "-", date, ".RData"))
    }

    results <-
     ifelse(this_fb, "fb_results_per_subregion", "results_per_subregion")

  } else {
    load(paste0(dir, pre, "region_", agent, "-", date, ".RData"))
    results <-
     ifelse(this_fb, "fb_results_per_region", "results_per_region")
  }

  ## return required results
  return(get(results))
}


## get summary per (sub)region for one agent
get_total <-
function(agent, what = "daly", rate = FALSE,
         by = c("subregion", "region"), age = c("all", "<5", "5+"),
         fb = FALSE) {
  ## check 'by' and 'age'
  by <- match.arg(by)
  age <- match.arg(age)

  ## load required results
  results <- load_results(agent, by, fb)

  ## return summary
  return(summarize(results, what, rate, age))
}

## get samples by (sub)region for one agent
get_samples <-
function(agent,
         what = c("daly", "yld", "yll", "incidence", "cases", "deaths"),
         by = c("subregion", "region"), fb = FALSE, ...) {
  print(agent[1])

  ## check 'what', 'by'
  what <- match.arg(what)
  by <- match.arg(by)

  ## load required results
  results <- load_results(agent, by, fb)

  ## get samples
  samples <-
    if (!any(is.na(results))) {
      get_global(results, what, ...)
    } else {
      NA
    }

  ## correction for Trichinella
  if (agent[1] == "Trichinella") {
    samples <- samples / 1e4
  }

  ## return samples
  return(list(samples = samples, agent = agent))
}


## get samples by (sub)region for more than one agent
get_all_samples <-
function(agents,
         what = c("daly", "yld", "yll", "incidence", "cases", "deaths"),
         by = c("subregion", "region"), fb = FALSE, ...) {
  ## check 'what', 'by'
  what <- match.arg(what)
  by <- match.arg(by)

  ## get samples
  samples <-
    lapply(agents,
           function(x)
             do.call(get_samples,
                     list(agent = x, what = what, by = by, fb = fb, ...))$samples)

  ## return samples
  return(samples)
}

## get DALY samples by (sub)region for one agent per outcome
get_samples_per_outcome <-
function(agent, dir, date, denom = 1e5,
         by = c("subregion", "region"), outcomes) {
  by <- match.arg(by)

  ## placeholders for 'results_per_region' and 'results_per_subregion'
  results_per_region <- NULL
  results_per_subregion <- NULL

  load(paste(dir, "/DALY_", by, "_", agent, "-", date, ".RData", sep = ""))
  if (by == "subregion") {
    results <- results_per_subregion
    n_regions <- 14
  } else if (by == "region"){
    results <- results_per_region
    n_regions <- 6
  }

  contrib <- get_contrib(results[[1]])
  if (missing(outcomes)) outcomes <- names(contrib)

  samples <- vector("list", length(contrib) + 1)
  samples <- lapply(samples, function(x) numeric(1000))
  names(samples) <- c(names(contrib), agent)

  for (r in seq(n_regions)) {
    for (n in seq(contrib)) {
      samples[[n]] <-
        samples[[n]] +
        rowSums(results[[1]][[r]]@samples[[n]][[contrib[n]]])
      samples[[length(contrib) + 1]] <-
        samples[[length(contrib) + 1]] +
        rowSums(results[[1]][[r]]@samples[[n]][[contrib[n]]])
    }
  }

  samples <- lapply(samples, function(x) x / denom)

  return(list(samples = samples, agent = agent, outcomes = outcomes))
}

## sum samples of different agents
sum_agents <-
function(...) {
  agents <- list(...)

  mx <-
    matrix(unlist(lapply(agents, function(x) x[[1]])), ncol = length(agents))
  return(rowSums(mx))
}

## ------------------------------------------------------------------------#
## SUMMARIZE AGENTS - HELPER FUNCTIONS ------------------------------------#

summarize_all <-
function(agents,
         what = "daly", rate = F, by = "region", age = "all", fb = F) {
  out <-
    lapply(agents,
           get_total,
           what = what,
           rate = rate,
           by = by,
           age = age,
           fb = fb)

  ## correction for Trichinella
  if (what != "daly_case" & any(names(out) == "p_tri")) {
    out$p_tri <- out$p_tri / 1e4
  }

  return(out)
}

extract <-
function(x, rows = NULL, cols = NULL) {
  if (is.null(rows)) rows <- seq(nrow(x[[1]]))
  if (is.null(cols)) cols <- seq(ncol(x[[1]]))
  return(sapply(x, function(x) x[rows, cols]))
}


## ------------------------------------------------------------------------#
## GET AGE-SPECIFIC LE FOR A GIVEN COUNTRY --------------------------------#

get_local_LE <-
function(country) {
  if (class(country) == "character") 
    country <- which(crpop$Country == country)
  return(local_LE[[country]])
}