#' Write masked data data frame (tibble) object to a CSV (comma separated file),
#' optionally writing associated meta-data to a JSON file. In particular, meta-data
#' in this case is defined as the attributes of the data frame object.
#'
#' @param md a masked data frame
#' @param filename filename for csv
#' @param write.metadata write a separate
#'
#' @export
md_write_csv <- function(md, filename, write.metadata=T)
{
    readr::write_csv(md, filename)
    if (write.metadata)
    {
        metadata <- attributes(md)
        metadata[["dataset"]] <- c(basename(filename))
        metadata[["row.names"]] <- NULL
        metadata[["names"]] <- NULL
        metadata[["class"]] <- NULL

        metadata.out <- paste(tools::file_path_sans_ext(filename),"json",sep=".")
        jsonlite::write_json(metadata, metadata.out, pretty=T)
    }
}

#' Read masked data from a JSON file.
#' If the JSON file has a 'dataset' field,
#' then each member of this field is assumed
#' to refer to a CSV file to read a masked
#' data sample from.
#'
#' Any metadata in the JSON
#' file is inserted into the attributes
#' of the masked data samples.
#'
#' @param filename filename for csv
#' @return list of masked data objects
#'
#' @export
md_read_json <- function(filename)
{
    metadata <- jsonlite::read_json(filename)
    dataset <- metadata[["dataset"]]

    mds <- list()
    for (data in dataset)
    {
        data.path <- file.path(dirname(filename),data)
        md <- readr::read_csv(data.path,col_types=list(k="i",w="i"))
        tmp <- metadata
        tmp[["dataset"]] <- c(data)
        attributes(md) <- c(attributes(md),tmp)
        class(md) <- c("tbl_md",class(md))

        mds[[data]] <- md
    }
    mds
}

#' Convert the columns corresponding to the
#' candidate matrix to a matrix object.
#'
#' @param md masked data
#' @return Candidate sets represented as a Boolean matrix
#'
#' @export
md_candidates_as_matrix <- function(md)
{
    c <- dplyr::select(md, dplyr::starts_with("c."))
    if (ncol(c) == 0)
        NA
    else
    {
        names(c) <- NULL
        as.matrix(c)
    }
}

#' Convert the columns corresponding to the
#' node times matrix to a matrix object.
#'
#' @param md masked data
#' @return Node times represented as a real matrix
#'
#' @export
md_node_times_as_matrix <- function(md)
{
    t <- dplyr::select(md, dplyr::starts_with("t."))
    if (ncol(t) == 0)
        NA
    else
    {
        names(t) <- NULL
        as.matrix(t)
    }
}

#' Retrieve the number of nodes implicitly
#' defined by the masked data input 'md'.
#'
#' @param md masked data
#' @return number of nodes in the series system
#'
#' @export
md_num_nodes <- function(md)
{
    if (!is.null(attr(md,"m")))
        as.integer(attr(md,"m"))
    else
        ncol(md_candidates_as_matrix(md))
}

#' Test whether \code{x} is masked data
#'
#' An object is considered to be masked data if
#' it is a type of data frame (e.g., tibble)
#' and it has at least two columns for candidate
#' sets named \code{c.1} and \code{c.2} and a column for
#' system failure time named \code{s}.
#'
#' @param x object to test
#' @export
md_is_masked_data <- function(x)
{
    is.data.frame(x) && all(c("c.1","c.2","s") %in% colnames(x))
}

#' Generates masked data for a series system with the given node failure times
#' \code{t}, candidate set model \code{candidate_model}, and candidate
#' set sizes \code{w}.
#'
#' @param t matrix of node failure times
#' @param w Integer vector. For the ith observation, generate w_j candidates.
#' @param candidate_model Function that accepts masked data as an argument.
#'                        The candidate model, defaults to md_candidate_m0.
#'                        If set to NULL, then do not generate a candidate set.
#'                        md_mle_exp_series will treat such masked data as
#'                        a sample that includes every node as candidates.
#' @return masked data, a data frame of n observations, (s,k,t1,...,tm,c1,...,cm)
#'         where k, t, and c are covariates (or predictors) of s,k,t1,...,tm.
#' @importFrom dplyr %>%
#' @export
md_series_data <- function(t,w,candidate_model=md_candidate_m0)
{
    m <- ncol(t)
    t <- tibble::as_tibble(t)
    names(t) <- paste("t",1:m,sep=".")

    md <- tibble::tibble(
        s = apply(t,1,min),
        k = apply(t,1,which.min),
        w = w)
    md <- dplyr::bind_cols(md,t)
    if (!is.null(candidate_model))
        md <- candidate_model(md,m)

    class(md) <- c("tbl_md",class(md))
    md
}

#' Candidate matrix to stringified vector of integers
#'
#' @param md masked data
#' @export
md_candidates_to_strings <- function(md)
{
    stopifnot(md_is_masked_data(md))

    C <- md_candidates_as_matrix(md)
    m <- md_num_nodes(md)

    cand_str <- character(nrow(md))
    for (i in 1:nrow(md))
        cand_str[i] <- toString((1:m)[C[i,]])

    cand_str
}

#' Print method for masked data (tbl_md).
#'
#' @param x masked data to print
#' @param pprint Boolean, show candidates as a string column
#' @param drop_latent Boolean, drop the latent random variables
#' @importFrom dplyr %>%
#' @export
print.tbl_md <- function(x,pprint=F,drop_latent=F,...)
{
    if (drop_latent)
        x <- x %>% dplyr::select(-c("k")) %>%
                   dplyr::select(-dplyr::starts_with("t.")) %>%
                   dplyr::select(-attr(x,"latent"))

    if (pprint)
        x <- x %>% dplyr::mutate("C" = md_candidates_to_strings(x)) %>%
                   dplyr::select(-dplyr::starts_with("c."))

    class(x) <- class(x)[class(x) != "tbl_md"]
    print(x,...)
}
