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
    write_csv(md, filename)
    if (write.metadata)
    {
        metadata <- attributes(md)
        metadata[["dataset"]] <- c(basename(filename))
        metadata[["row.names"]] <- NULL
        metadata[["names"]] <- NULL
        metadata[["class"]] <- NULL

        metadata.out <- paste(tools::file_path_sans_ext(filename),"json",sep=".")
        write_json(metadata, metadata.out, pretty=T)
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
    metadata <- read_json(filename)
    dataset <- metadata[["dataset"]]

    mds <- list()
    for (data in dataset)
    {
        data.path <- file.path(dirname(filename),data)
        md <- read_csv(data.path,col_types=list(k="i",w="i"))
        tmp <- metadata
        tmp[["dataset"]] <- c(data)
        attributes(md) <- c(attributes(md),tmp)
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
    c <- select(md, starts_with("c."))
    if (ncol(c) == 0)
        NA
    else
        as.matrix(c)
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
    t <- select(md, starts_with("t."))
    if (ncol(t) == 0)
        NA
    else
        as.matrix(t)
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
    ncol(md_candidates_matrix(md))
}
