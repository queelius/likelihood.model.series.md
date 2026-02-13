#' Decode a matrix from prefixed columns in a data frame
#'
#' Extracts columns matching the pattern `var1, var2, ...` or `var.1, var.2, ...`
#' from a data frame and returns them as a numeric matrix, ordered by index.
#'
#' @param df data frame containing the matrix columns
#' @param var character prefix for the column names
#' @return a matrix, or NULL if no matching columns are found
#' @keywords internal
md_decode_matrix <- function(df, var) {
  stopifnot(is.data.frame(df), is.character(var))
  pat <- paste0("^", var, "\\.?(\\d+)$")
  cols <- grep(pat, colnames(df), value = TRUE)
  if (length(cols) == 0L) return(NULL)
  rank <- as.integer(sub(pat, "\\1", cols))
  as.matrix(df[cols[order(rank)]])
}

#' Encode a matrix as a data frame with prefixed column names
#'
#' Converts a matrix to a data frame with columns named `var1, var2, ...`.
#'
#' @param mat matrix to encode
#' @param var character prefix for the column names
#' @return a data frame with named columns
#' @export
md_encode_matrix <- function(mat, var) {
  stopifnot(is.matrix(mat))
  df <- as.data.frame(mat)
  names(df) <- paste0(var, seq_len(ncol(mat)))
  df
}

#' Mark columns as latent in a masked data frame
#'
#' Sets the `"latent"` attribute on a data frame, recording which columns
#' represent unobserved (latent) variables.
#'
#' @param md data frame to modify
#' @param vars character vector of column names to mark as latent
#' @return the data frame with updated latent attribute
#' @keywords internal
md_mark_latent <- function(md, vars) {
  stopifnot(is.character(vars))
  attr(md, "latent") <- union(vars, attr(md, "latent"))
  md
}

#' Convert Boolean candidate set columns to character set notation
#'
#' Replaces Boolean matrix columns (e.g., `x1, x2, x3`) with a single
#' character column showing set notation like `{1, 3}`.
#'
#' @param df data frame containing Boolean matrix columns
#' @param setvar column prefix for the Boolean matrix (default `"x"`)
#' @param cname name for the new character column (default: same as `setvar`)
#' @param drop_set if TRUE, remove the original Boolean columns (default FALSE)
#' @return data frame with character set column added
#' @export
md_boolean_matrix_to_charsets <- function(df, setvar = "x", cname = NULL,
                                          drop_set = FALSE) {
  if (is.null(cname)) cname <- setvar
  mat <- md_decode_matrix(df, setvar)
  stopifnot(!is.null(mat))

  if (drop_set) {
    drop_cols <- grep(paste0("^", setvar, "\\.?\\d+$"), colnames(df), value = TRUE)
    df <- df[, !(colnames(df) %in% drop_cols), drop = FALSE]
  }

  df[[cname]] <- apply(mat, 1, function(row) {
    paste0("{", paste(which(row), collapse = ", "), "}")
  })
  df
}
