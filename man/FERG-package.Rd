\name{FERG-package}
\alias{FERG-package}
\alias{FERG}

\docType{package}

\title{DALY Calculation Framework for WHO/FERG}

\description{
The \pkg{FERG} package provides functions for calculating and reporting DALYs in the framework of WHO's Foodborne Disease Burden Epidemiology Reference Group work towards the global burden of foodborne disease. The DALY calculation framework of the \pkg{FERG} package is based on the \pkg{DALY} package.
}

\details{
  \tabular{ll}{
    Package:    \tab  FERG\cr
    Type:       \tab  Package\cr
    Version:    \tab  0.1.0\cr
    Date:       \tab  2015-07-16\cr
    Authors:    \tab  Brecht Devleesschauwer, Scott McDonald\cr
    Maintainer: \tab  Brecht Devleesschauwer <\email{brechtdv@gmail.com}>\cr
    Depends:    \tab  R (>= 3.0.0), methods\cr
    Imports:    \tab  XLConnect, knitr, xtable, lme4, rjags, ggplot2, prevalence, rworldmap, RColorBrewer\cr
    License:    \tab  GPL (>= 2)\cr
    LazyLoad:   \tab  yes\cr
}

Available functions in the \pkg{FERG} package:
  \tabular{ll}{
    \strong{Imputation Activity}\cr
    \code{\link{extractIncidence}}     \tab Get incidence data from Excel database file using \code{XLConnect}.\cr
    \code{\link{check_db}}             \tab Check raw imputation database.\cr
    \code{\link{impute}}               \tab Impute missing values based on a Random Effects model.\cr
    \code{\link{merge}}                \tab Merge observed with imputed values.\cr
    \code{\link{save_merged_db}}       \tab Save merged values to workbook.\cr
    \code{\link{plot_missing}}         \tab Plot missing values using \code{rworldmap}.\cr
    \code{\link{summarize_imputation}} \tab Summarize imputation results.\cr
    \code{\link{get_total_cases}}      \tab Obtain cases by WHO subregion & total.\cr
    \code{\link{imputation_report}}    \tab Generate report of imputation results.\cr\cr

    \strong{DALY Calculation Activity}\cr
    \code{\link{readDatabase}}         \tab Read entire Excel workbook using \code{XLConnect}.\cr
    \code{\link{check_XL}}             \tab Check if sheets are correctly ordered.\cr
    \code{\link{read_XL}}              \tab Read Excel and transform to 'DB' list.\cr
    \code{\link{getDALY_per_country}}  \tab Calculate DALYs for each country.\cr
    \code{\link{getDALY_per_region}}   \tab Aggregate DALYs per region.\cr
    \code{\link{getDALY_per_subregion}}\tab Aggregate DALYs per subregion.\cr
    \code{\link{check_convergence}}    \tab Check MC convergence.\cr
    \code{\link{DALY_map}}             \tab Produce global burden of disease map.\cr
    \code{\link{FERG_report}}          \tab Produce FERG PDF report.\cr
    \code{\link{get_contrib}}          \tab Derive node contributions (ie, YLD or YLL)\cr
    \code{\link{get_local_LE}}         \tab Obtain the life expectancy table of a given country\cr
    \code{\link{mean_age}}             \tab Calculate mean age from age distribution.\cr\cr

    \strong{Source Attribution Activity}\cr
    \code{\link{generate_samples}}     \tab Generate normalized samples.\cr
    \code{\link{normalize_by_food}}    \tab Normalize samples by proportion foodborne.\cr
    \code{\link{summary_stats}}        \tab Summarize samples.\cr
    \code{\link{sa_barplot}}           \tab Generate source attribution barplot.\cr
    \code{\link{sa_plot}}              \tab Generate source attribution line plot.\cr
    \code{\link{sa_report}}            \tab Generate source attribution report.\cr
    \code{\link{apply_prop_fb}}        \tab Apply proportion foodborne.\cr
    \code{\link{get_global}}           \tab Get global results.\cr\cr

    \strong{FERG Results :: helper functions}\cr
    \code{\link{load_results}}         \tab Load (fb_)results_per_region/subregion.\cr
    \code{\link{get_total}}            \tab Get summary per (sub)region for 1 agent.\cr
    \code{\link{get_samples}}          \tab Get DALY samples by (sub)region for 1 agent.\cr
    \code{\link{get_all_samples}}      \tab Get DALY samples by (sub)region for >1 agent.\cr
    \code{\link{get_samples_per_outcome}}   \tab Get all samples per outcome.\cr
    \code{\link{summarize}}            \tab Summarize 'DALYrun_region' object.\cr
    \code{\link{summarize_all}}        \tab Get summary for all agents.\cr\cr  
    
    \strong{FERG Results :: plot functions}\cr
    \code{\link{line_plot}}            \tab Median and UI, per agent.\cr
    \code{\link{full_barplot}}         \tab Median rate, per agent & region.\cr
    \code{\link{region_barplot}}       \tab YLD/YLL, ages, scaled to 100\%, per agent.\cr
    \code{\link{scatter_plot}}         \tab DALY/case ~ DALY/100,000, per agent.\cr\cr    
  }
}

\keyword{package}
