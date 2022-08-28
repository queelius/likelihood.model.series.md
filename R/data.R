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

#' Bootstrapped MLE sampling distribution statistics from masked data, compared
#' with the asymptotic theory, for exponential series system with 3 components.
#'
#' @format A data frame.
#' \describe{
#'   \item{n}{sample size}
#'   \item{asymptotic.mse}{asymptotic mean squared error of the MLE},
#'   \item{boot.mse}{estimate of mean squared error using Bootstrap method}
#'   \item{asymptotic.rate1.se}{asymptotic standard error of MLE for parameter rate1}
#'   \item{boot.rate1.se}{estimate of standard error of MLE for parameter rate1 using Bootstrap method}
#'   \item{asymptotic.rate2.se}{asymptotic standard error of MLE for parameter rate2}
#'   \item{boot.rate2.se}{estimate of standard error of MLE for parameter rate2 using Bootstrap method}
#'   \item{asymptotic.rate3.se}{asymptotic standard error of MLE for parameter rate3}
#'   \item{boot.rate3.se}{estimate of standard error of MLE for parameter rate3 using Bootstrap method}
#'   \item{rate1}{MLE for rate1}
#'   \item{rate1.bias}{estimate of bias of rate1 using Boostrap method}
#'   \item{rate2}{MLE for rate2}
#'   \item{rate2.bias}{estimate of bias of rate2 using Boostrap method}
#'   \item{rate3}{MLE for rate3}
#'   \item{rate3.bias}{estimate of bias of rate3 using Boostrap method}
#' }
#' @source \url{https://github.com/queelius/series_system_estimation_masked_data/blob/master/data-raw/exp_series_stats_1.R}
"exp_series_stats_1"

#' MLE sampling distribution statistics, compared with the
#' asymptotic theory, for exponential series system with 3 components.
#'
#' @format A data frame.
#' \describe{
#'   \item{n}{sample size}
#'   \item{mse}{asymptotic mean squared error of the MLE},
#'   \item{frob}{frobenius norm of the difference vcov(theta.hat) - vcov(theta) given `md`}
#'   \item{rate1.lb}{lower-bound 95% confidence interval for rate1}
#'   \item{rate1.ub}{upper-bound 95% confidence interval for rate1}
#'   \item{rate1.in}{indicator of whether confidence interval contains rate1}
#'   \item{rate2.lb}{lower-bound 95% confidence interval for rate2}
#'   \item{rate2.ub}{upper-bound 95% confidence interval for rate2}
#'   \item{rate2.in}{indicator of whether confidence interval contains rate2}
#'   \item{rate3.lb}{lower-bound 95% confidence interval for rate3}
#'   \item{rate3.ub}{upper-bound 95% confidence interval for rate3}
#'   \item{rate3.in}{indicator of whether confidence interval contains rate3}
#'   \item{rate4.lb}{lower-bound 95% confidence interval for rate4}
#'   \item{rate4.ub}{upper-bound 95% confidence interval for rate4}
#'   \item{rate4.in}{indicator of whether confidence interval contains rate4}
#' }
#' @source \url{https://github.com/queelius/series_system_estimation_masked_data/blob/master/data-raw/exp_series_stats_2.R}
"exp_series_stats_2"


#' Simulated masked data 1 for weibull series
#'
#' Simulated masked data for weibull series system satisfying conditions
#' C1, C2, and C3. See source url for information about the parameters of the
#' simulation.
#'
#' @format A data frame with 300 rows and 10 variables, only 5 of which are observable.
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
#' @source \url{https://github.com/queelius/series_system_estimation_masked_data/blob/master/data-raw/weibull_series_md_1.R}
"weibull_series_md_1"
