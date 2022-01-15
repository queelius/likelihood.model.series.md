#' Write masked data data frame object to a CSV file, optionally writing
#' an associated metadata json file where the metadata is the contents of
#' the attributes of the data frame object.
#'
#' @param md a masked data frame
#' @param filename filename for csv
#' @param write.meta write a separate
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
md.write <- function(md, filename, write.meta=T)
{
    library(jsonlite)
    library(tidyverse)

    write_csv(md, filename)
    if (write.meta)
    {
        metadata <- attributes(md)
        metadata[["dataset"]] <- c(basename(filename))
        metadata[["row.names"]] <- NULL
        metadata[["names"]] <- NULL
        metadata[["class"]] <- NULL

        meta.out <- paste(tools::file_path_sans_ext(filename),"json",sep=".")
        write_json(metadata, meta.out, pretty=T)
    }
}

md.read.from.json <- function(filename)
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

md.candidates.matrix <- function(md)
{
    c <- select(md, starts_with("c."))
    if (ncol(c) == 0) NA
    else              as.matrix(c)
}

md.node.times.matrix <- function(md)
{
    t <- select(md, starts_with("t."))
    if (ncol(t) == 0) NA
    else              as.matrix(t)
}

md.nnodes <- function(md)
{
    ncol(md.candidates.matrix(md))
}
