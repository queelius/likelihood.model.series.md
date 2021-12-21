
#' Generate masked data from series system data
#'
#' @param data a data frame with column 's' system failure time,
#'                               column 'k' failed component, and
#'                               column 't' vector of component lifetimes.
#' @param w a vector of sizes for candidate sets
#' @param alpha a vector of probabilities
#'
#' @return masked data
#' @export
#'
#' @examples
#' n = 100
#' t = rseries.exp(n,c(1,2,3))
#' w = rep(2,n)
#' alpha = rep(.9,n)
#' data = rmasked.data(t,w,alpha)
rmasked.data = function(data,w=NULL,alpha=NULL)
{
  m = ncol(data)-2
  t = data[,2:m]
  n = nrow(t)
  K = data$K
  S = data$S
  C = matrix(nrow=n,ncol=m,F)

  if (is.null(alpha))
    alpha = rep(1,n)

  if (is.null(w))
    w = rep(m,n)

  test = stats::runif(n) < alpha

  for (i in 1:n)
  {
    if (test[i])
    {
      C[i,K[i]] = T
      C[i,sample((1:m)[-K[i]],size=w[i]-1,replace=F)] = T
    }
    else
      C[i,sample((1:m)[-K[i]],size=w[i],replace=F)] = T
  }

  df = data.frame(S,alpha,C)
  names(df) = c("S","alpha",paste("C",sep="_",1:m))
  df
}
