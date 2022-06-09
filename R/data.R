#' @importFrom tibble tibble
NULL

#' Simulated masked data 1
#'
#' Simulated masked data for exponential series system satisfying conditions
#' C1, C2, and C3. See source url for information about the parameters of the
#' simulation.
#'
#' @format A data frame with 200 rows and 10 variables, only 5 of which are observable.
#' \describe{
#'   \item{t}{Real unobservable variable, system lifetime}
#'   \item{s}{Real observable variable, right-censored system lifetime}
#'   \item{delta}{Boolean observable variable, indicates whether system lifetime has been censored}
#'   \item{k}{Integer unobservable variable, index of the failed component}
#'   \item{t1}{Real unobservable variable, lifetime of component 1}
#'   \item{t2}{Real unobservable variable, lifetime of component 2}
#'   \item{t3}{Real unobservable variable, lifetime of component 3}
#'   \item{x1}{Boolean observable variable, TRUE indicates component 1 is in candidate set}
#'   \item{x2}{Boolean observable variable, TRUE indicates component 2 is in candidate set}
#'   \item{x3}{Boolean observable variable, TRUE indicates component 3 is in candidate set}
#' }
#' @source \url{https://github.com/queelius/series_system_estimation_masked_data/blob/master/data-raw/exp_series_md_1.R}
"exp_series_md_1"

#' Real-world masked data. The source is from the paper,
#' "Estimating Component Reliabilities from Incomplete System Failure Data",
#' Table 2 (Example Data for a Series System).
#'
#' They assume the data is compatible with a series system of 3 components
#' w/Weibull distributed lifetimes and candidate sets that satisfy conditions
#' C1, C2, and C3.
#'
#' @format A data frame.
#' \describe{
#'   \item{t}{Real unobservable variable, system lifetime}
#'   \item{x1}{Boolean observable variable, TRUE indicates component 1 is in candidate set}
#'   \item{x2}{Boolean observable variable, TRUE indicates component 2 is in candidate set}
#'   \item{x3}{Boolean observable variable, TRUE indicates component 3 is in candidate set}
#' }
#' @source \url{https://github.com/queelius/series_system_estimation_masked_data/blob/master/data-raw/guo_weibull_series_md.R}
"guo_weibull_series_md"
