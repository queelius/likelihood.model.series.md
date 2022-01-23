# Retrieve the function arguments.
md_func_args <- function(...)
{
    call <- evalq(match.call(expand.dots = F), parent.frame(1))
    formals <- evalq(formals(), parent.frame(1))

    for(i in setdiff(names(formals), names(call)))
        call[i] <- list( formals[[i]] )

    match.call(sys.function(sys.parent()), call)
}

# pseudo-inverse
pinv <- function(m)
{
    s <- svd(m);
    d <- rep(0, length(s$d));
    for (i in 1:length(d)) { d[i] <- ifelse(s$d[i] == 0, 0, 1/s$d[i]);  }
    s$v %*% diag(d) %*% t(s$u); }
