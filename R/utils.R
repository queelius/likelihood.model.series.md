# Retrieve the function arguments.
md_func_args <- function(...)
{
    call <- evalq(match.call(expand.dots = F), parent.frame(1))
    formals <- evalq(formals(), parent.frame(1))

    for(i in setdiff(names(formals), names(call)))
        call[i] <- list( formals[[i]] )

    match.call(sys.function(sys.parent()), call)
}

