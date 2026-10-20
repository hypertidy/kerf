#' Read one batch
#'
#' Materialises one plan row. The plan decides *what* is read; the backend
#' decides *how*.
#'
#' For a raster plan (one with a window group):
#'
#' * `backend = "terra"` returns a `SpatRaster` with its window set to the
#'   batch, so nothing is read until values are needed.
#' * `backend = "gdalraster"` reads values now and returns a numeric array
#'   `[nrow, ncol, band]` carrying a `"kerf_window"` attribute
#'   `c(xoff, yoff, xsize, ysize)`. Arrays serialise without special
#'   handling, which suits pipeline targets.
#'
#' For a vector plan, the owned features are read through an attribute filter
#' built from the plan's FID range, FID set or `where` predicate. With
#' `role = "context"` and a context extent, every feature whose bounding box
#' meets that extent is read instead. Either way the result carries a logical
#' column `kerf_owned` marking the batch's own features, which [kerf_trim()]
#' keeps.
#'
#' * `backend = "terra"` returns a `SpatVector`, with the FID in `kerf_fid`.
#' * `backend = "gdalraster"` returns the data frame from `GDALVector$fetch()`,
#'   with the FID in `FID` and geometry as WKB. It serialises without special
#'   handling.
#'
#' @param x the source: for rasters a `SpatRaster` (terra backend only), a
#'   character data source name, or a `"kerf_grid"` with a `dsn`; for vectors
#'   a data source name or a `"kerf_vector"` made from one.
#' @param plan a plan with exactly one row.
#' @param role `"context"` (the default) reads the context when the row has
#'   one, otherwise the owned selection; `"owned"` reads the owned selection.
#' @param backend `"terra"` or `"gdalraster"`.
#' @param band for the gdalraster raster backend, the bands to read; default
#'   all.
#'
#' @return a `SpatRaster`, array, `SpatVector` or data frame, see Details.
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
  if (!is.data.frame(plan) || nrow(plan) != 1L) {
    stop("'plan' must have exactly one row", call. = FALSE)
  }
  if (!"xoff" %in% names(plan)) {
    return(read_vector(x, plan, role, backend))
  }
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

#' Trim a batch result to what it owns
#'
#' Cuts a result read with context back to the batch's own members: the owned
#' window of a raster, or the features marked `kerf_owned` in a vector result
#' (the marker column is dropped). A raster result read without context is
#' returned unchanged.
#'
#' @param x a result from [kerf_read()], or any result with the same geometry
#'   or the same `kerf_owned` column.
#' @param plan a plan with exactly one row.
#'
#' @return `x` cut to what the batch owns.
#' @export
#' @examples
#' if (requireNamespace("terra", quietly = TRUE)) {
#'   f <- system.file("ex/elev.tif", package = "terra")
#'   p <- kerf_plan(f, kerf_grid(2, 2), context = kerf_halo(2))
#'   r <- kerf_read(f, p[1, ])
#'   kerf_trim(r, p[1, ])
#' }
kerf_trim <- function(x, plan) {
  if (is_vector_result(x)) {
    keep <- as.logical(x$kerf_owned)
    keep[is.na(keep)] <- FALSE
    out <- x[keep, ]
    out$kerf_owned <- NULL
    return(out)
  }
  if (!"ctx_xoff" %in% names(plan) || is.na(plan$ctx_xoff[1L])) return(x)
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
#' Trims each part to what its batch owns and assembles the parts. Because
#' owned selections partition the input, context never double counts.
#'
#' @param parts a list of results, one per plan row and in plan order:
#'   `SpatRaster`s or arrays for raster plans, `SpatVector`s or data frames
#'   for vector plans, as returned by [kerf_read()].
#' @param plan the plan the parts were made from.
#'
#' @return a single merged `SpatRaster`; one array covering the union of owned
#'   windows, with a `"kerf_window"` attribute; or the owned rows of every
#'   part bound together, as a `SpatVector` or data frame.
#' @export
kerf_reduce <- function(parts, plan) {
  if (length(parts) != nrow(plan)) {
    stop("'parts' must have one element per plan row", call. = FALSE)
  }
  trimmed <- lapply(seq_along(parts), function(i) {
    kerf_trim(parts[[i]], plan[i, , drop = FALSE])
  })
  if (all(vapply(trimmed, inherits, logical(1), "SpatVector"))) {
    return(do.call(rbind, unname(trimmed)))
  }
  if (all(vapply(trimmed, is.data.frame, logical(1)))) {
    out <- do.call(rbind, unname(trimmed))
    rownames(out) <- NULL
    return(out)
  }
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
  stop("'parts' must all be of one kind: SpatRasters, arrays, ",
       "SpatVectors or data frames", call. = FALSE)
}

## ---- internals -------------------------------------------------------------

## c(xoff, yoff, xsize, ysize) for one plan row
plan_window <- function(plan, role) {
  if (nrow(plan) != 1L) stop("'plan' must have exactly one row", call. = FALSE)
  use_ctx <- role == "context" && "ctx_xoff" %in% names(plan) &&
    !is.na(plan$ctx_xoff[1L])
  prefix <- if (use_ctx) "ctx_" else ""
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

is_vector_result <- function(x) {
  (inherits(x, "SpatVector") || is.data.frame(x)) &&
    "kerf_owned" %in% names(x)
}

has_ctx_extent <- function(plan) {
  all(c("ctx_xmin", "ctx_ymin", "ctx_xmax", "ctx_ymax") %in% names(plan)) &&
    !is.na(plan$ctx_xmin[1L])
}

## OGR SQL for the owned selection of one plan row, or NULL
owned_sql <- function(plan, prefix = "") {
  col <- function(nm) plan[[paste0(prefix, nm)]]
  has <- function(nms) all(paste0(prefix, nms) %in% names(plan)) &&
    !all(is_na_cell(col(nms[1L])))
  parts <- character()
  if (has(c("fid_lo", "fid_hi"))) {
    parts <- c(parts, sprintf("FID >= %s AND FID <= %s",
                              format_num_sql(col("fid_lo")[1L]),
                              format_num_sql(col("fid_hi")[1L])))
  }
  if (has("fids")) {
    f <- col("fids")[[1L]]
    parts <- c(parts, if (length(f)) {
      sprintf("FID IN (%s)", paste(format_num_sql(f), collapse = ", "))
    } else {
      "1 = 0"
    })
  }
  if (has("where")) parts <- c(parts, col("where")[1L])
  if (!length(parts)) return(NULL)
  if (length(parts) > 1L) parts <- sprintf("(%s)", parts)
  paste(parts, collapse = " AND ")
}

read_vector <- function(x, plan, role, backend) {
  src <- kerf_vector(x)
  if (is.na(src$dsn)) {
    stop("reading features needs a data source name; a feature index ",
         "holds no geometry", call. = FALSE)
  }
  use_ctx <- role == "context" && has_ctx_extent(plan)
  where <- owned_sql(plan)
  if (is.null(where)) {
    stop("'plan' has no owned selection kerf can read (fid, fids or where)",
         call. = FALSE)
  }
  box <- if (use_ctx) {
    c(plan$ctx_xmin, plan$ctx_ymin, plan$ctx_xmax, plan$ctx_ymax)
  }
  if (backend == "gdalraster") {
    lyr <- open_layer(src$dsn, src$layer)
    on.exit(lyr$close(), add = TRUE)
    if (use_ctx) {
      lyr$setSpatialFilterRect(box)
    } else {
      lyr$setAttributeFilter(where)
    }
    out <- lyr$fetch(-1)
    fid <- as.numeric(out$FID)
  } else {
    q <- sprintf("SELECT FID AS kerf_fid, * FROM \"%s\"", src$layer)
    if (use_ctx) {
      out <- terra::vect(src$dsn, query = q,
                         extent = terra::ext(box[1L], box[3L], box[2L],
                                             box[4L]))
    } else {
      out <- terra::vect(src$dsn, query = paste(q, "WHERE", where))
    }
    fid <- as.numeric(out$kerf_fid)
  }
  out$kerf_owned <- if (use_ctx) owned_fid(fid, plan, src, where) else
    rep(TRUE, length(fid))
  out
}

## which of fid the plan row owns
owned_fid <- function(fid, plan, src, where) {
  if ("fids" %in% names(plan)) return(fid %in% plan$fids[[1L]])
  if (all(c("fid_lo", "fid_hi") %in% names(plan))) {
    own <- fid >= plan$fid_lo & fid <= plan$fid_hi
    if (!"where" %in% names(plan)) return(own)
  }
  ## a predicate: ask GDAL which FIDs it selects
  lyr <- open_layer(src$dsn, src$layer)
  on.exit(lyr$close(), add = TRUE)
  ## skip geometry; attributes stay readable because the predicate needs them
  g <- lyr$getGeometryColumn()
  lyr$setIgnoredFields(if (nzchar(g)) g else "OGR_GEOMETRY")
  lyr$setAttributeFilter(where)
  fid %in% as.numeric(lyr$fetch(-1)$FID)
}
