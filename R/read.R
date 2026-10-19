#' Read one batch
#'
#' Materialises one plan row. The plan decides *what* is read; the backend
#' decides *how*.
#'
#' * `backend = "terra"` returns a `SpatRaster` with its window set to the
#'   batch, so nothing is read until values are needed.
#' * `backend = "gdalraster"` reads values now and returns a numeric array
#'   `[nrow, ncol, band]` carrying a `"kerf_window"` attribute
#'   `c(xoff, yoff, xsize, ysize)`. Arrays serialise without special
#'   handling, which suits pipeline targets.
#'
#' @param x the source: a `SpatRaster` (terra backend only), a character data
#'   source name, or a `"kerf_grid"` with a `dsn`.
#' @param plan a plan with exactly one row and a window group.
#' @param role `"context"` (the default) reads the context window when the
#'   plan has one, otherwise the owned window; `"owned"` reads the owned
#'   window.
#' @param backend `"terra"` or `"gdalraster"`.
#' @param band for the gdalraster backend, the bands to read; default all.
#'
#' @return a `SpatRaster` or an array, see Details.
#' @export
#' @examples
#' if (requireNamespace("terra", quietly = TRUE)) {
#'   f <- system.file("ex/elev.tif", package = "terra")
#'   p <- kerf_plan(f, kerf_grid(2, 2), context = kerf_halo(2))
#'   kerf_read(f, p[1, ])
#' }
kerf_read <- function(x, plan, role = c("context", "owned"),
                      backend = c("terra", "gdalraster"), band = NULL) {
  role <- match.arg(role)
  backend <- match.arg(backend)
  w <- plan_window(plan, role)
  if (backend == "terra") {
    r <- terra_source(x)
    src <- kerf_raster(r)
    e <- window_extent(src, w[1L], w[2L], w[3L], w[4L])
    out <- c(r)
    terra::window(out) <- terra::ext(e$xmin, e$xmax, e$ymin, e$ymax)
    return(out)
  }
  dsn <- source_dsn(x)
  ds <- methods::new(gdalraster::GDALRaster, dsn, read_only = TRUE)
  on.exit(ds$close(), add = TRUE)
  if (is.null(band)) band <- seq_len(ds$getRasterCount())
  out <- array(NA_real_, c(w[4L], w[3L], length(band)))
  for (i in seq_along(band)) {
    v <- ds$read(band[i], w[1L], w[2L], w[3L], w[4L], w[3L], w[4L])
    out[, , i] <- matrix(as.numeric(v), nrow = w[4L], ncol = w[3L],
                         byrow = TRUE)
  }
  attr(out, "kerf_window") <- w
  out
}

#' Trim a batch result to its owned window
#'
#' Cuts a result read with context back to the cells the batch owns. A result
#' that was read without context is returned unchanged.
#'
#' @param x a `SpatRaster` or an array from [kerf_read()], or any result
#'   with the same geometry as the context read.
#' @param plan a plan with exactly one row and a window group.
#'
#' @return `x` cut to the owned window.
#' @export
#' @examples
#' if (requireNamespace("terra", quietly = TRUE)) {
#'   f <- system.file("ex/elev.tif", package = "terra")
#'   p <- kerf_plan(f, kerf_grid(2, 2), context = kerf_halo(2))
#'   r <- kerf_read(f, p[1, ])
#'   kerf_trim(r, p[1, ])
#' }
kerf_trim <- function(x, plan) {
  if (!"ctx_xoff" %in% names(plan)) return(x)
  own <- plan_window(plan, "owned")
  ctx <- plan_window(plan, "context")
  dx <- own[1L] - ctx[1L]
  dy <- own[2L] - ctx[2L]
  if (inherits(x, "SpatRaster")) {
    ex <- as.vector(terra::ext(x))
    rs <- terra::res(x)
    e <- terra::ext(ex[1L] + dx * rs[1L], ex[1L] + (dx + own[3L]) * rs[1L],
                    ex[4L] - (dy + own[4L]) * rs[2L], ex[4L] - dy * rs[2L])
    return(terra::crop(x, e, snap = "near"))
  }
  if (is.array(x)) {
    out <- x[dy + seq_len(own[4L]), dx + seq_len(own[3L]), , drop = FALSE]
    attr(out, "kerf_window") <- own
    return(out)
  }
  stop(sprintf("don't know how to trim an object of class '%s'",
               class(x)[1L]), call. = FALSE)
}

#' Recombine batch results
#'
#' Trims each part to its owned window and assembles the parts. Because owned
#' windows partition the grid, halos never double count.
#'
#' @param parts a list of results, one per plan row and in plan order:
#'   `SpatRaster`s or arrays as returned by [kerf_read()].
#' @param plan the plan the parts were made from.
#'
#' @return for `SpatRaster` parts, a single merged `SpatRaster`; for arrays,
#'   one array covering the union of owned windows, with a `"kerf_window"`
#'   attribute.
#' @export
kerf_reduce <- function(parts, plan) {
  if (length(parts) != nrow(plan)) {
    stop("'parts' must have one element per plan row", call. = FALSE)
  }
  trimmed <- lapply(seq_along(parts), function(i) {
    kerf_trim(parts[[i]], plan[i, , drop = FALSE])
  })
  if (all(vapply(trimmed, inherits, logical(1), "SpatRaster"))) {
    if (length(trimmed) == 1L) return(trimmed[[1L]])
    return(do.call(terra::merge, unname(trimmed)))
  }
  if (all(vapply(trimmed, is.array, logical(1)))) {
    x0 <- min(plan$xoff)
    y0 <- min(plan$yoff)
    nx <- max(plan$xoff + plan$xsize) - x0
    ny <- max(plan$yoff + plan$ysize) - y0
    nb <- dim(trimmed[[1L]])[3L]
    out <- array(NA_real_, c(ny, nx, nb))
    for (i in seq_along(trimmed)) {
      ii <- plan$yoff[i] - y0 + seq_len(plan$ysize[i])
      jj <- plan$xoff[i] - x0 + seq_len(plan$xsize[i])
      out[ii, jj, ] <- trimmed[[i]]
    }
    attr(out, "kerf_window") <- c(x0, y0, nx, ny)
    return(out)
  }
  stop("'parts' must all be SpatRasters or all arrays", call. = FALSE)
}

## ---- internals -------------------------------------------------------------

## c(xoff, yoff, xsize, ysize) for one plan row
plan_window <- function(plan, role) {
  if (nrow(plan) != 1L) stop("'plan' must have exactly one row", call. = FALSE)
  prefix <- if (role == "context" && "ctx_xoff" %in% names(plan)) "ctx_" else ""
  cols <- paste0(prefix, c("xoff", "yoff", "xsize", "ysize"))
  if (!all(cols %in% names(plan))) {
    stop("'plan' needs a window group (xoff, yoff, xsize, ysize)",
         call. = FALSE)
  }
  as.integer(unlist(plan[1L, cols], use.names = FALSE))
}

terra_source <- function(x) {
  if (inherits(x, "SpatRaster")) return(x)
  terra::rast(source_dsn(x))
}

source_dsn <- function(x) {
  if (is.character(x) && length(x) == 1L) return(x)
  if (inherits(x, "kerf_grid") && !is.na(x$dsn)) return(x$dsn)
  if (inherits(x, "SpatRaster")) {
    s <- terra::sources(x)
    if (length(s) && nzchar(s[1L])) return(s[1L])
    stop("this SpatRaster is in memory; use backend = \"terra\"",
         call. = FALSE)
  }
  stop("'x' must be a data source name, a kerf_grid with a dsn, ",
       "or a SpatRaster", call. = FALSE)
}
