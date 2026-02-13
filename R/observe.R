#' Right-censoring observation scheme
#'
#' Creates an observation functor that applies right-censoring at time \code{tau}.
#' Systems that fail before \code{tau} are observed exactly; systems surviving
#' past \code{tau} are right-censored.
#'
#' @param tau censoring time (positive numeric)
#' @return A function with signature \code{function(t_true)} returning a list
#'   with components:
#'   \describe{
#'     \item{t}{observed time}
#'     \item{omega}{"exact" or "right"}
#'     \item{t_upper}{NA (not used for this scheme)}
#'   }
#' @export
#' @examples
#' obs <- observe_right_censor(tau = 100)
#' obs(50)   # exact: list(t = 50, omega = "exact", t_upper = NA)
#' obs(150)  # right-censored: list(t = 100, omega = "right", t_upper = NA)
observe_right_censor <- function(tau) {
  function(t_true) {
    if (t_true <= tau) {
      list(t = t_true, omega = "exact", t_upper = NA_real_)
    } else {
      list(t = tau, omega = "right", t_upper = NA_real_)
    }
  }
}


#' Left-censoring observation scheme (single inspection)
#'
#' Creates an observation functor for a single-inspection design. If the system
#' has already failed by inspection time \code{tau}, we know it failed before
#' \code{tau} but not exactly when (left-censored). If it is still running, we
#' know it survived past \code{tau} (right-censored).
#'
#' @param tau inspection time (positive numeric)
#' @return A function with signature \code{function(t_true)} returning a list
#'   with components:
#'   \describe{
#'     \item{t}{inspection time \code{tau}}
#'     \item{omega}{"left" if failed before \code{tau}, "right" otherwise}
#'     \item{t_upper}{NA (not used for this scheme)}
#'   }
#' @export
#' @examples
#' obs <- observe_left_censor(tau = 100)
#' obs(50)   # left-censored: list(t = 100, omega = "left", t_upper = NA)
#' obs(150)  # right-censored: list(t = 100, omega = "right", t_upper = NA)
observe_left_censor <- function(tau) {
  function(t_true) {
    if (t_true <= tau) {
      list(t = tau, omega = "left", t_upper = NA_real_)
    } else {
      list(t = tau, omega = "right", t_upper = NA_real_)
    }
  }
}


#' Periodic inspection observation scheme
#'
#' Creates an observation functor for periodic inspections at intervals of
#' \code{delta}. Failures are bracketed between the last inspection before
#' failure and the first inspection after failure (interval-censored). Systems
#' surviving past \code{tau} are right-censored.
#'
#' @param delta inspection interval (positive numeric)
#' @param tau study end time (positive numeric or Inf for no right-censoring)
#' @return A function with signature \code{function(t_true)} returning a list
#'   with components:
#'   \describe{
#'     \item{t}{lower bound of interval (or \code{tau} if right-censored)}
#'     \item{omega}{"interval" or "right"}
#'     \item{t_upper}{upper bound of interval (NA if right-censored)}
#'   }
#' @export
#' @examples
#' obs <- observe_periodic(delta = 10, tau = 100)
#' obs(25)   # interval: list(t = 20, omega = "interval", t_upper = 30)
#' obs(150)  # right-censored: list(t = 100, omega = "right", t_upper = NA)
observe_periodic <- function(delta, tau = Inf) {
  function(t_true) {
    if (t_true > tau) {
      return(list(t = tau, omega = "right", t_upper = NA_real_))
    }
    lower <- floor(t_true / delta) * delta
    upper <- lower + delta
    list(t = lower, omega = "interval", t_upper = upper)
  }
}


#' Mixture of observation schemes
#'
#' Creates an observation functor that randomly selects from a set of
#' observation schemes for each observation. This models heterogeneous
#' monitoring environments where different units are observed differently.
#'
#' @param ... observation functors (created by \code{observe_*} functions)
#' @param weights mixing probabilities (numeric vector). If NULL, uniform
#'   weights are used. Weights are normalized to sum to 1.
#' @return A function with signature \code{function(t_true)} returning a list
#'   from one of the constituent schemes, selected randomly according to
#'   \code{weights}.
#' @export
#' @examples
#' obs <- observe_mixture(
#'   observe_right_censor(tau = 100),
#'   observe_left_censor(tau = 50),
#'   weights = c(0.7, 0.3)
#' )
#' set.seed(42)
#' obs(30)  # randomly selects one of the two schemes
observe_mixture <- function(..., weights = NULL) {
  schemes <- list(...)
  if (length(schemes) == 0) stop("at least one observation scheme is required")
  if (is.null(weights)) weights <- rep(1, length(schemes))
  if (length(weights) != length(schemes)) {
    stop("weights must have the same length as the number of schemes")
  }
  weights <- weights / sum(weights)
  function(t_true) {
    idx <- sample.int(length(schemes), 1L, prob = weights)
    schemes[[idx]](t_true)
  }
}
