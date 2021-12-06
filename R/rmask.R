rmask = function(t,w=NULL,alpha=NULL)
{
  n = nrow(t)
  m = ncol(t)
  K = apply(t,1,which.min)
  S = apply(t,1,min)
  C = matrix(nrow=n,ncol=m,F)

  if (is.null(alpha))
    alpha = rep(1,n)

  if (is.null(w))
    w = rep(m,n)

  test = runif(n) < alpha

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
