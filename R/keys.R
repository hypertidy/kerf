## Key generators and context forms. Each is a small description; nothing
## is computed until kerf_plan() resolves it against a source.

new_keys <- function(type, ...) {
  structure(list(type = type, params = list(...)), class = "kerf_keys")
}

new_context <- function(type, ...) {
  structure(list(type = type, params = list(...)), class = "kerf_context")
}

check_count <- function(x, name) {
  if (length(x) != 1L || is.na(x) || !is.numeric(x) || x < 1 ||
      x != round(x)) {
    stop(sprintf("'%s' must be a single positive whole number", name),
         call. = FALSE)
  }
  as.integer(x)
}

#' Raster key generators
#'
#' Key generators describe how a raster is cut into batches. They are
#' resolved against a source by [kerf_plan()].
#'
#' * `kerf_blocks()` cuts along the source's native blocks, `n_row` blocks
#'   high and `n_col` blocks wide per batch. Edge batches are clipped to the
#'   grid. This matches `geotargets::tile_blocksize()`.
#' * `kerf_tile_size()` cuts into tiles of a fixed pixel size, clipped at the
#'   edges.
#' * `kerf_grid()` splits the grid into `ncol` by `nrow` batches whose edges
#'   are snapped to the nearest cell. This matches `geotargets::tile_grid()`.
#' * `kerf_n_tiles()` picks a near-square `ncol` by `nrow` split giving about
#'   `n` batches. This matches `geotargets::tile_n()`.
#'
#' Batches are ordered row-major from the top-left. Every raster plan owns
#' cells by the `centre` rule, so owned windows partition the grid.
#'
#' `kerf_grid()` places the edge between batches `i` and `i + 1` at the cell
#' boundary nearest to `i * n / k`, rounding exact half-cell ties up (right
#' and down). terra settles those exact ties by floating-point rounding, so
#' `geotargets::tile_grid()` and `geotargets::tile_n()` may put such an edge
#' one cell the other way; every other edge is identical.
#'
#' @param n_row,n_col number of native blocks per batch, vertically and
#'   horizontally.
#' @param ncol,nrow for `kerf_tile_size()`, the tile size in pixels; for
#'   `kerf_grid()`, the number of batch columns and rows.
#' @param n approximate total number of batches.
#'
#' @return a `"kerf_keys"` object.
#' @name kerf_keys
#' @examples
#' g <- kerf_grid_spec(c(100, 80), blocksize = c(32, 16))
#' kerf_plan(g, kerf_blocks(2, 1))
#' kerf_plan(g, kerf_grid(3, 2))
NULL

#' @rdname kerf_keys
#' @export
kerf_blocks <- function(n_row = 1, n_col = 1) {
  new_keys("blocks", n_row = check_count(n_row, "n_row"),
           n_col = check_count(n_col, "n_col"))
}

#' @rdname kerf_keys
#' @export
kerf_tile_size <- function(ncol, nrow = ncol) {
  new_keys("tile_size", ncol = check_count(ncol, "ncol"),
           nrow = check_count(nrow, "nrow"))
}

#' @rdname kerf_keys
#' @export
kerf_grid <- function(ncol, nrow) {
  new_keys("grid", ncol = check_count(ncol, "ncol"),
           nrow = check_count(nrow, "nrow"))
}

#' @rdname kerf_keys
#' @export
kerf_n_tiles <- function(n) {
  n <- check_count(n, "n")
  ## the same near-square choice as geotargets::tile_n()
  nrow <- floor(sqrt(n))
  q <- n / nrow
  ncol <- if (q == round(q)) q else ceiling(q)
  new_keys("grid", ncol = as.integer(ncol), nrow = as.integer(nrow))
}

#' @export
print.kerf_keys <- function(x, ...) {
  p <- x$params
  cat(sprintf("<kerf_keys> %s(%s)\n", x$type,
              paste(names(p), unlist(p), sep = " = ", collapse = ", ")))
  invisible(x)
}

#' Raster context: a halo of cells
#'
#' A halo widens each batch's read window by a number of cells on every side,
#' clipped to the grid. The halo is context only: [kerf_trim()] and
#' [kerf_reduce()] cut results back to the owned window, so focal and other
#' neighbourhood operations give the same answer as on the whole grid.
#'
#' @param cells number of cells, one value for both directions or
#'   `c(ncol, nrow)`.
#'
#' @return a `"kerf_context"` object.
#' @export
#' @examples
#' g <- kerf_grid_spec(c(100, 80))
#' kerf_plan(g, kerf_grid(2, 2), context = kerf_halo(3))
kerf_halo <- function(cells) {
  cells <- rep_len(as.numeric(cells), 2L)
  if (anyNA(cells) || any(cells < 0) || any(cells != round(cells))) {
    stop("'cells' must be whole numbers >= 0", call. = FALSE)
  }
  new_context("halo", ncol = as.integer(cells[1L]),
              nrow = as.integer(cells[2L]))
}

#' @export
print.kerf_context <- function(x, ...) {
  p <- x$params
  cat(sprintf("<kerf_context> %s(%s)\n", x$type,
              paste(names(p), unlist(p), sep = " = ", collapse = ", ")))
  invisible(x)
}
