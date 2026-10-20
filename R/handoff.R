#' Plan rows as extents
#'
#' Returns one named numeric vector `c(xmin, xmax, ymin, ymax)` per plan row,
#' the form `terra::ext()` accepts and the form the 'geotargets' tile helpers
#' return, so a plan can feed `tar_terra_tiles()` as it stands.
#'
#' @param plan a plan with an extent group.
#' @param role `"owned"` (the default) or `"context"`.
#'
#' @return a list of named numeric vectors.
#' @export
#' @examples
#' g <- kerf_grid_spec(c(100, 80), c(0, 1000, 0, 800))
#' kerf_extents(kerf_plan(g, kerf_grid(2, 1)))
kerf_extents <- function(plan, role = c("owned", "context")) {
  role <- match.arg(role)
  prefix <- if (role == "context") "ctx_" else ""
  cols <- paste0(prefix, c("xmin", "xmax", "ymin", "ymax"))
  if (!all(cols %in% names(plan))) {
    stop(sprintf("'plan' has no %s extent group", role), call. = FALSE)
  }
  m <- as.matrix(plan[cols])
  colnames(m) <- c("xmin", "xmax", "ymin", "ymax")
  lapply(seq_len(nrow(m)), function(i) m[i, ])
}

#' Plan row as GDAL command line arguments
#'
#' Turns one plan row into arguments for the GDAL command line interface: a
#' window becomes `--window` for `gdal raster clip`; an extent with no
#' window becomes `--bbox`, accepted by `gdal raster clip` and
#' `gdal vector filter`; an FID range, an FID set and a predicate become one
#' `--where` for `gdal vector filter` (an FID set is written out as
#' `FID IN (...)`, so very large sets make long arguments).
#'
#' @param plan a plan with exactly one row.
#' @param role `"owned"` (the default) or `"context"`.
#'
#' @return a character vector of arguments.
#' @export
#' @examples
#' g <- kerf_grid_spec(c(100, 80), c(0, 1000, 0, 800))
#' p <- kerf_plan(g, kerf_grid(2, 2), context = kerf_halo(4))
#' kerf_gdal_args(p[1, ])
#' kerf_gdal_args(p[1, ], role = "context")
kerf_gdal_args <- function(plan, role = c("owned", "context")) {
  role <- match.arg(role)
  if (nrow(plan) != 1L) stop("'plan' must have exactly one row", call. = FALSE)
  prefix <- if (role == "context") "ctx_" else ""
  col <- function(nm) plan[[paste0(prefix, nm)]]
  has <- function(nms) all(paste0(prefix, nms) %in% names(plan)) &&
    !is.na(col(nms[1L])[1L])
  args <- character()
  if (has(c("xoff", "yoff", "xsize", "ysize"))) {
    args <- c(args, "--window",
              paste(col("xoff"), col("yoff"), col("xsize"), col("ysize"),
                    sep = ","))
  } else if (has(c("xmin", "ymin", "xmax", "ymax"))) {
    args <- c(args, "--bbox",
              paste(sprintf("%.17g", c(col("xmin"), col("ymin"),
                                       col("xmax"), col("ymax"))),
                    collapse = ","))
  }
  where <- owned_sql(plan, prefix)
  if (!is.null(where)) args <- c(args, "--where", where)
  if (!length(args)) {
    stop(sprintf("the %s selection has no window, extent, fid or where group",
                 role), call. = FALSE)
  }
  args
}
