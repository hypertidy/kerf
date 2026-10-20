#' Vector source handle from metadata
#'
#' A vector source is a data source name and layer, read through
#' 'gdalraster', or a data frame that already is a feature index (one row per
#' feature with `fid, xmin, ymin, xmax, ymax`, plus any attribute columns a
#' plan needs). A data frame index lets plans be made from any tool that
#' emits feature bounding boxes, for example a slim `gdal vector sql`
#' projection.
#'
#' @param x a character data source name, a feature index data frame, or a
#'   `"kerf_vector"`.
#' @param layer layer name or index; default the first layer.
#' @param ... ignored.
#'
#' @return a `"kerf_vector"` list.
#' @export
#' @examples
#' idx <- data.frame(fid = 1:3, xmin = c(0, 5, 9), ymin = 0,
#'                   xmax = c(1, 6, 10), ymax = 1)
#' kerf_vector(idx)
kerf_vector <- function(x, layer = NULL, ...) {
  UseMethod("kerf_vector")
}

#' @export
kerf_vector.kerf_vector <- function(x, layer = NULL, ...) {
  x
}

#' @export
kerf_vector.data.frame <- function(x, layer = NULL, ...) {
  need <- c("fid", "xmin", "ymin", "xmax", "ymax")
  miss <- setdiff(need, names(x))
  if (length(miss)) {
    stop("a feature index needs columns: ", paste(miss, collapse = ", "),
         call. = FALSE)
  }
  x$fid <- as.numeric(x$fid)
  if (anyNA(x$fid) || anyDuplicated(x$fid)) {
    stop("'fid' must be unique and not missing", call. = FALSE)
  }
  ok <- is.finite(x$xmin) & is.finite(x$ymin) &
    is.finite(x$xmax) & is.finite(x$ymax)
  ext <- if (any(ok)) {
    c(min(x$xmin[ok]), max(x$xmax[ok]), min(x$ymin[ok]), max(x$ymax[ok]))
  } else {
    rep(NA_real_, 4L)
  }
  structure(
    list(dsn = NA_character_, layer = NA_character_, crs = NA_character_,
         n = nrow(x), extent = ext, index = x),
    class = "kerf_vector"
  )
}

#' @export
kerf_vector.character <- function(x, layer = NULL, ...) {
  if (length(x) != 1L || is.na(x)) {
    stop("'x' must be a single data source name", call. = FALSE)
  }
  lyr <- open_layer(x, layer)
  on.exit(lyr$close(), add = TRUE)
  bb <- lyr$bbox()
  structure(
    list(dsn = x, layer = lyr$getName(),
         crs = lyr$getSpatialRef(),
         n = as.numeric(lyr$getFeatureCount()),
         extent = c(bb[1L], bb[3L], bb[2L], bb[4L]),
         index = NULL),
    class = "kerf_vector"
  )
}

#' @export
kerf_vector.default <- function(x, layer = NULL, ...) {
  stop(sprintf("no kerf_vector() method for class '%s'", class(x)[1L]),
       call. = FALSE)
}

#' @export
print.kerf_vector <- function(x, ...) {
  cat("<kerf_vector>\n")
  if (!is.na(x$dsn)) {
    cat(sprintf("        dsn: %s\n", x$dsn))
    cat(sprintf("      layer: %s\n", x$layer))
  } else {
    cat("     source: feature index data frame\n")
  }
  cat(sprintf("   features: %s\n", format(x$n, big.mark = ",")))
  cat(sprintf("     extent: %s  (xmin, xmax, ymin, ymax)\n",
              paste(format(x$extent), collapse = ", ")))
  invisible(x)
}

#' Feature index: FID and bounding box per feature
#'
#' The planning read for vector sources. Only FIDs and geometry are read
#' (attribute fields are ignored), and each geometry is reduced to its
#' bounding box. Features with empty or missing geometry get `NA` boxes.
#'
#' @param x a vector source, anything [kerf_vector()] accepts.
#' @param fields attribute fields to read as well, for example the field an
#'   attribute plan groups by.
#'
#' @return a data frame with columns `fid, xmin, ymin, xmax, ymax` and any
#'   requested fields, ordered by `fid`.
#' @export
kerf_vector_index <- function(x, fields = NULL) {
  src <- kerf_vector(x)
  if (!is.null(src$index)) {
    miss <- setdiff(fields, names(src$index))
    if (length(miss)) {
      stop("the feature index has no field(s): ", paste(miss, collapse = ", "),
           call. = FALSE)
    }
    idx <- src$index
    return(idx[order(idx$fid), , drop = FALSE])
  }
  lyr <- open_layer(src$dsn, src$layer)
  on.exit(lyr$close(), add = TRUE)
  all_fields <- layer_fields(lyr)
  miss <- setdiff(fields, all_fields)
  if (length(miss)) {
    stop("the layer has no field(s): ", paste(miss, collapse = ", "),
         call. = FALSE)
  }
  ignore <- setdiff(all_fields, fields)
  if (length(ignore)) lyr$setIgnoredFields(ignore)
  lyr$returnGeomAs <- "WKB"
  d <- lyr$fetch(-1)
  gcol <- attr(d, "gis")$geom_column
  if (!length(gcol) || !nzchar(gcol[1L])) gcol <- lyr$defaultGeomColName
  env <- unclass(wk::wk_envelope(wk::wkb(d[[gcol]])))
  idx <- data.frame(fid = as.numeric(d$FID),
                    xmin = env$xmin, ymin = env$ymin,
                    xmax = env$xmax, ymax = env$ymax)
  bad <- !(is.finite(idx$xmin) & is.finite(idx$ymin) &
             is.finite(idx$xmax) & is.finite(idx$ymax))
  idx[bad, c("xmin", "ymin", "xmax", "ymax")] <- NA_real_
  for (f in fields) idx[[f]] <- d[[f]]
  idx[order(idx$fid), , drop = FALSE]
}

## ---- internals -------------------------------------------------------------

open_layer <- function(dsn, layer = NULL) {
  if (!requireNamespace("gdalraster", quietly = TRUE)) {
    stop("reading vector sources needs 'gdalraster'", call. = FALSE)
  }
  if (is.null(layer) || is.na(layer)) {
    methods::new(gdalraster::GDALVector, dsn)
  } else {
    methods::new(gdalraster::GDALVector, dsn, layer)
  }
}

## attribute field names of a layer (geometry fields excluded)
layer_fields <- function(lyr) {
  defn <- lyr$getLayerDefn()
  is_geom <- vapply(defn, function(d) isTRUE(d$is_geom), logical(1))
  names(defn)[!is_geom]
}

## Hilbert index of integer cell coordinates on a 2^level grid
hilbert_index <- function(x, y, level) {
  n <- 2^level
  x <- as.integer(x)
  y <- as.integer(y)
  d <- numeric(length(x))
  s <- as.integer(n / 2)
  while (s >= 1L) {
    rx <- as.integer(bitwAnd(x, s) > 0L)
    ry <- as.integer(bitwAnd(y, s) > 0L)
    d <- d + as.numeric(s) * s * bitwXor(3L * rx, ry)
    ## rotate the quadrant
    flip <- ry == 0L & rx == 1L
    x[flip] <- as.integer(n - 1 - x[flip])
    y[flip] <- as.integer(n - 1 - y[flip])
    swap <- ry == 0L
    tmp <- x[swap]
    x[swap] <- y[swap]
    y[swap] <- tmp
    s <- s %/% 2L
  }
  d
}

## Hilbert keys of bounding-box centres within a key extent
## c(xmin, xmax, ymin, ymax); centres outside the extent are clamped
centre_keys <- function(idx, extent, level) {
  cx <- (idx$xmin + idx$xmax) / 2
  cy <- (idx$ymin + idx$ymax) / 2
  n <- 2^level
  w <- extent[2L] - extent[1L]
  h <- extent[4L] - extent[3L]
  ix <- if (w > 0) floor((cx - extent[1L]) / w * n) else rep(0, length(cx))
  iy <- if (h > 0) floor((cy - extent[3L]) / h * n) else rep(0, length(cy))
  ix <- pmin(pmax(ix, 0), n - 1)
  iy <- pmin(pmax(iy, 0), n - 1)
  key <- rep(NA_real_, length(cx))
  ok <- !is.na(ix) & !is.na(iy)
  key[ok] <- hilbert_index(ix[ok], iy[ok], level)
  key
}
