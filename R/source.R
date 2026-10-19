#' Describe a raster grid
#'
#' A grid spec is the metadata a plan is made from: dimension, extent, and
#' optionally the CRS, the native block size and the data source name. No
#' pixel values are involved.
#'
#' @param dimension integer `c(ncol, nrow)`.
#' @param extent numeric `c(xmin, xmax, ymin, ymax)`. Defaults to pixel
#'   coordinates `c(0, ncol, 0, nrow)`.
#' @param crs character CRS description, stored but not used.
#' @param blocksize integer `c(block_ncol, block_nrow)`, the native block size
#'   of the source, or `NULL` if unknown.
#' @param dsn character data source name, or `NA`.
#'
#' @return a `"kerf_grid"` list.
#' @export
#' @examples
#' kerf_grid_spec(c(100, 80), c(0, 100, 0, 80), blocksize = c(32, 16))
kerf_grid_spec <- function(dimension, extent = NULL, crs = NA_character_,
                           blocksize = NULL, dsn = NA_character_) {
  dimension <- as.integer(dimension)
  if (length(dimension) != 2L || anyNA(dimension) || any(dimension < 1L)) {
    stop("'dimension' must be two positive integers, c(ncol, nrow)",
         call. = FALSE)
  }
  if (is.null(extent)) extent <- c(0, dimension[1L], 0, dimension[2L])
  extent <- as.numeric(extent)
  if (length(extent) != 4L || anyNA(extent) ||
      extent[2L] <= extent[1L] || extent[4L] <= extent[3L]) {
    stop("'extent' must be c(xmin, xmax, ymin, ymax) with xmax > xmin ",
         "and ymax > ymin", call. = FALSE)
  }
  if (!is.null(blocksize)) {
    blocksize <- as.integer(blocksize)
    if (length(blocksize) != 2L || anyNA(blocksize) || any(blocksize < 1L)) {
      stop("'blocksize' must be two positive integers, ",
           "c(block_ncol, block_nrow)", call. = FALSE)
    }
  }
  structure(
    list(dimension = dimension, extent = extent,
         crs = as.character(crs), blocksize = blocksize,
         dsn = as.character(dsn)),
    class = "kerf_grid"
  )
}

#' @export
print.kerf_grid <- function(x, ...) {
  cat("<kerf_grid>\n")
  cat(sprintf("  dimension: %i x %i  (ncol, nrow)\n",
              x$dimension[1L], x$dimension[2L]))
  cat(sprintf("     extent: %s  (xmin, xmax, ymin, ymax)\n",
              paste(format(x$extent), collapse = ", ")))
  if (!is.null(x$blocksize)) {
    cat(sprintf("  blocksize: %i x %i\n", x$blocksize[1L], x$blocksize[2L]))
  }
  if (!is.na(x$dsn)) cat(sprintf("        dsn: %s\n", x$dsn))
  invisible(x)
}

#' Raster source handle from metadata
#'
#' Reads only metadata. A `SpatRaster` uses 'terra'; a character data source
#' name uses 'gdalraster' when it is installed, otherwise 'terra'.
#'
#' @param x a `SpatRaster`, a character data source name, or a `"kerf_grid"`.
#' @param ... ignored.
#'
#' @return a `"kerf_grid"`, see [kerf_grid_spec()].
#' @export
#' @examples
#' if (requireNamespace("terra", quietly = TRUE)) {
#'   f <- system.file("ex/elev.tif", package = "terra")
#'   kerf_raster(f)
#' }
kerf_raster <- function(x, ...) {
  UseMethod("kerf_raster")
}

#' @export
kerf_raster.kerf_grid <- function(x, ...) {
  x
}

#' @export
kerf_raster.SpatRaster <- function(x, ...) {
  dm <- c(terra::ncol(x), terra::nrow(x))
  ex <- as.vector(terra::ext(x))
  bs <- terra_blocksize(x)
  src <- terra::sources(x)
  dsn <- if (length(src) && nzchar(src[1L])) src[1L] else NA_character_
  kerf_grid_spec(dm, ex, crs = terra::crs(x), blocksize = bs, dsn = dsn)
}

#' @export
kerf_raster.character <- function(x, ...) {
  if (length(x) != 1L || is.na(x)) {
    stop("'x' must be a single data source name", call. = FALSE)
  }
  if (requireNamespace("gdalraster", quietly = TRUE)) {
    ds <- methods::new(gdalraster::GDALRaster, x, read_only = TRUE)
    on.exit(ds$close(), add = TRUE)
    bb <- ds$bbox()
    return(kerf_grid_spec(
      c(ds$getRasterXSize(), ds$getRasterYSize()),
      c(bb[1L], bb[3L], bb[2L], bb[4L]),
      crs = ds$getProjection(),
      blocksize = ds$getBlockSize(1L),
      dsn = x
    ))
  }
  if (requireNamespace("terra", quietly = TRUE)) {
    g <- kerf_raster(terra::rast(x))
    g$dsn <- x
    return(g)
  }
  stop("reading raster metadata needs 'gdalraster' or 'terra'", call. = FALSE)
}

#' @export
kerf_raster.default <- function(x, ...) {
  stop(sprintf("no kerf_raster() method for class '%s'", class(x)[1L]),
       call. = FALSE)
}

## native block size of a file-backed SpatRaster as c(ncol, nrow), or NULL
terra_blocksize <- function(x) {
  b <- tryCatch(terra::fileBlocksize(x), error = function(e) NULL)
  if (is.null(b) || !length(b)) return(NULL)
  bs <- c(b[1L, "cols"], b[1L, "rows"])
  ## in-memory rasters report no usable block size
  if (anyNA(bs) || any(bs < 1)) return(NULL)
  bs
}
