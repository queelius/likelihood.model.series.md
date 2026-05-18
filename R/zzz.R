.onLoad <- function(libname, pkgname) {
    # Force eager load of dist.structure now, then re-register our
    # polymorphic mean.exp_series so it wins the S3 dispatch table.
    # Without this: dist.structure stays lazy (it is in Imports but not
    # @importFrom-ed), then the first dist.structure:: call from any
    # maskedcauses function loads dist.structure for real, whose
    # NAMESPACE re-registers a monomorphic mean.exp_series that breaks
    # the legacy structure(rates, class = "exp_series") path.
    requireNamespace("dist.structure", quietly = TRUE)
    registerS3method("mean", "exp_series", mean.exp_series,
                     envir = asNamespace(pkgname))
}
