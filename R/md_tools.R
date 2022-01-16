#' Write masked data data frame (tibble) object to a CSV (comma separated file),
#' optionally writing associated meta-data to a JSON file. In particular, meta-data
#' in this case is defined as the attributes of the data frame object.
#'
#' @param md a masked data frame
#' @param filename filename for csv
#' @param write.metadata write a separate
#'
#' @export
#'
#' @examples
#' n <- 100
#' w <- rep(2,n)
#' md <- rmd.exp.series.m0(n,theta,w)
#' md.write(md,"data.csv")
#' md.test <- md.read.from.json("data.json")
#' assert(md == md.test)
md_write_csv <- function(md, filename, write.metadata=T)
{
    library(jsonlite)
    library(tidyverse)

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

md_read_json <- function(filename)
{
    library(tidyverse)
    library(jsonlite)
    meta <- read_json(filename)
    dataset <- meta[["dataset"]]

    mds <- list()
    for (data in dataset)
    {
        data.path <- file.path(dirname(filename),data)
        md <- read_csv(data.path,col_types=list(k="i",w="i"))
        meta.tmp <- meta
        meta.tmp[["dataset"]] <- c(data)
        attributes(md) <- c(attributes(md),meta.tmp)
        mds[[data]] <- md
    }
    mds
}

md_candidates_as_matrix <- function(md)
{
    c <- select(md, starts_with("c."))
    if (ncol(c) == 0) NA
    else              as.matrix(c)
}

md_node_times_as_matrix <- function(md)
{
    t <- select(md, starts_with("t."))
    if (ncol(t) == 0) NA
    else              as.matrix(t)
}

md_num_nodes <- function(md)
{
    ncol(md_candidates_matrix(md))
}
