#' Make a plan
#'
#' Resolves a key generator, and optionally a context, against a source and
#' returns the plan: a plain data frame with one row per batch.
#'
#' For a raster the columns are `batch`, the owned window `xoff, yoff,
#' xsize, ysize` (0-based pixel offsets and sizes, as GDAL uses them), the
#' owned extent `xmin, ymin, xmax, ymax`, the same groups prefixed `ctx_`
#' when a context is given, `rule` (always `"centre"`) and `n_owned` (cells).
#' See [kerf_contract()] for the full contract.
#'
#' `batch` is derived from the owned window's top-left offset, so it is
#' stable for a given source and generator and readable as a branch label.
#'
#' Because a plan is just a data frame, it can be filtered, sorted, written to
#' Parquet or built by hand; [kerf_check()] validates it either way.
#'
#' @param x a source: anything [kerf_raster()] accepts.
#' @param keys a key generator, see [kerf_keys].
#' @param context `NULL` for no context, or a context such as [kerf_halo()].
#'
#' @return a data frame.
#' @export
#' @examples
#' g <- kerf_grid_spec(c(100, 80), c(0, 1000, 0, 800), blocksize = c(32, 16))
#' p <- kerf_plan(g, kerf_blocks(2, 2), context = kerf_halo(4))
#' head(p)
#' kerf_check(p, g)
kerf_plan <- function(x, keys, context = NULL) {
  if (!inherits(keys, "kerf_keys")) {
    stop("'keys' must be a key generator such as kerf_blocks()",
         call. = FALSE)
  }
  if (!is.null(context) && !inherits(context, "kerf_context")) {
    stop("'context' must be NULL or a context such as kerf_halo()",
         call. = FALSE)
  }
  src <- kerf_raster(x)
  w <- resolve_windows(src, keys)

  plan <- data.frame(
    batch = sprintf("x%d_y%d", w$xoff, w$yoff),
    xoff = w$xoff, yoff = w$yoff, xsize = w$xsize, ysize = w$ysize,
    stringsAsFactors = FALSE
  )
  plan <- cbind(plan, window_extent(src, w$xoff, w$yoff, w$xsize, w$ysize))

  if (!is.null(context)) {
    cw <- resolve_context(src, w, context)
    ctx <- cbind(
      data.frame(xoff = cw$xoff, yoff = cw$yoff,
                 xsize = cw$xsize, ysize = cw$ysize),
      window_extent(src, cw$xoff, cw$yoff, cw$xsize, cw$ysize)
    )
    names(ctx) <- paste0("ctx_", names(ctx))
    plan <- cbind(plan, ctx)
  }
  plan$rule <- "centre"
  plan$n_owned <- as.numeric(w$xsize) * w$ysize
  rownames(plan) <- NULL
  plan
}

## ---- internals -------------------------------------------------------------

## owned windows, row-major from the top-left
resolve_windows <- function(src, keys) {
  dm <- src$dimension
  p <- keys$params
  switch(keys$type,
    blocks = {
      if (is.null(src$blocksize)) {
        stop("the source has no native block size; use kerf_tile_size() ",
             "or kerf_grid() instead", call. = FALSE)
      }
      size <- src$blocksize * c(p$n_col, p$n_row)
      grout_windows(src, size)
    },
    tile_size = grout_windows(src, c(p$ncol, p$nrow)),
    grid = split_windows(dm, p$ncol, p$nrow),
    stop(sprintf("unknown key generator type '%s'", keys$type),
         call. = FALSE)
  )
}

## fixed-size tiles via grout, clipped at the grid edge
grout_windows <- function(src, size) {
  size <- as.integer(pmin(size, src$dimension))
  g <- grout::grout(src$dimension, src$extent, blocksize = size)
  ti <- grout::tile_index(g)
  list(xoff = as.integer(ti$offset_x), yoff = as.integer(ti$offset_y),
       xsize = as.integer(ti$ncol), ysize = as.integer(ti$nrow))
}

## an ncol x nrow split with edges snapped to the nearest cell edge
split_windows <- function(dm, ncol, nrow) {
  xe <- split_edges(dm[1L], ncol)
  ye <- split_edges(dm[2L], nrow)
  windows_from_edges(xe, ye)
}

split_edges <- function(n, k) {
  e <- floor(seq(0, k) * n / k + 0.5)
  unique(as.integer(e))
}

windows_from_edges <- function(xe, ye) {
  nx <- length(xe) - 1L
  ny <- length(ye) - 1L
  ci <- rep(seq_len(nx), times = ny)
  ri <- rep(seq_len(ny), each = nx)
  list(xoff = xe[ci], yoff = ye[ri],
       xsize = xe[ci + 1L] - xe[ci], ysize = ye[ri + 1L] - ye[ri])
}

## context windows: the owned window widened and clipped to the grid
resolve_context <- function(src, w, context) {
  dm <- src$dimension
  switch(context$type,
    halo = {
      hx <- context$params$ncol
      hy <- context$params$nrow
      x0 <- pmax(0L, w$xoff - hx)
      y0 <- pmax(0L, w$yoff - hy)
      x1 <- pmin(dm[1L], w$xoff + w$xsize + hx)
      y1 <- pmin(dm[2L], w$yoff + w$ysize + hy)
      list(xoff = as.integer(x0), yoff = as.integer(y0),
           xsize = as.integer(x1 - x0), ysize = as.integer(y1 - y0))
    },
    stop(sprintf("unknown context type '%s'", context$type), call. = FALSE)
  )
}

## extent of pixel windows on a grid, in contract column order
window_extent <- function(src, xoff, yoff, xsize, ysize) {
  ex <- src$extent
  dm <- src$dimension
  rx <- (ex[2L] - ex[1L]) / dm[1L]
  ry <- (ex[4L] - ex[3L]) / dm[2L]
  data.frame(
    xmin = ex[1L] + xoff * rx,
    ymin = ex[4L] - (yoff + ysize) * ry,
    xmax = ex[1L] + (xoff + xsize) * rx,
    ymax = ex[4L] - yoff * ry
  )
}
