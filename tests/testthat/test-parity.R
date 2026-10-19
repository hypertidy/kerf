## kerf's raster generators reproduce terra::getTileExtents(), which is what
## the geotargets tile helpers wrap

as_matrix <- function(e) {
  m <- do.call(rbind, e)
  dimnames(m) <- NULL
  m
}

## distinct tile edges in pixel units, from extents
edges_x <- function(m, r) {
  round((sort(unique(c(m[, 1], m[, 2]))) - terra::xmin(r)) / terra::xres(r))
}
edges_y <- function(m, r) {
  round((terra::ymax(r) - sort(unique(c(m[, 3], m[, 4])),
                                decreasing = TRUE)) / terra::yres(r))
}

## each edge equals terra's, except an exact half-cell tie, which terra may
## round either way
expect_edges <- function(got, want, n, k) {
  p <- (0:k) * n / k
  tie <- abs(p - floor(p) - 0.5) < 1e-9
  expect_equal(got[!tie], want[!tie])
  expect_true(all(abs(got[tie] - want[tie]) <= 1))
}

test_that("kerf_grid matches terra across awkward dimensions", {
  skip_if_no_terra()
  set.seed(42)
  dims <- rbind(c(75, 53), c(10, 10), c(7, 300), cbind(
    sample(5:400, 40, replace = TRUE), sample(5:400, 40, replace = TRUE)))
  for (i in seq_len(nrow(dims))) {
    r <- terra::rast(ncols = dims[i, 1], nrows = dims[i, 2],
                     xmin = -1234.5, xmax = 987.25, ymin = -77, ymax = 412)
    for (k in list(c(2, 2), c(3, 1), c(1, 4), c(5, 3), c(4, 7))) {
      if (k[1] > dims[i, 1] || k[2] > dims[i, 2]) next
      p <- kerf_plan(r, kerf_grid(k[1], k[2]))
      template <- terra::rast(terra::ext(r), ncols = k[1], nrows = k[2])
      m <- terra::getTileExtents(r, template)
      expect_equal(nrow(p), nrow(m))
      g <- as_matrix(kerf_extents(p))
      expect_edges(edges_x(g, r), edges_x(m, r), dims[i, 1], k[1])
      expect_edges(edges_y(g, r), edges_y(m, r), dims[i, 2], k[2])
    }
  }
})

test_that("kerf_grid matches terra exactly when no edge is a tie", {
  skip_if_no_terra()
  r <- terra::rast(ncols = 72, nrows = 60, xmin = 0.1, xmax = 7.3,
                   ymin = -3, ymax = 3)
  for (k in list(c(2, 2), c(3, 4), c(4, 5), c(6, 6), c(1, 3))) {
    template <- terra::rast(terra::ext(r), ncols = k[1], nrows = k[2])
    expect_equal(as_matrix(kerf_extents(kerf_plan(r, kerf_grid(k[1], k[2])))),
                 unname(terra::getTileExtents(r, template)),
                 tolerance = 1e-9)
  }
})

test_that("kerf_blocks matches terra's block tiling", {
  skip_if_no_terra()
  f <- make_tif(75L, 53L, block = 16L)
  r <- terra::rast(f)
  for (k in list(c(1, 1), c(2, 1), c(1, 3), c(3, 2))) {
    got <- as_matrix(kerf_extents(kerf_plan(r, kerf_blocks(k[1], k[2]))))
    want <- unname(terra::getTileExtents(r, c(16, 16) * k))
    expect_equal(got, want, tolerance = 1e-9)
  }
})

test_that("kerf matches the geotargets tile helpers directly", {
  skip_if_no_terra()
  skip_if_not_installed("geotargets")
  ## 72 x 60 avoids half-cell ties for every split used here
  f <- make_tif(72L, 60L, block = 16L)
  r <- terra::rast(f)
  same <- function(p, tiles) {
    expect_equal(as_matrix(kerf_extents(p)),
                 unname(do.call(rbind, tiles)), tolerance = 1e-9)
  }
  same(kerf_plan(r, kerf_grid(3, 2)), geotargets::tile_grid(r, 3, 2))
  same(kerf_plan(r, kerf_blocks(2, 1)),
       geotargets::tile_blocksize(r, n_blocks_row = 2, n_blocks_col = 1))
  for (n in c(1, 4, 6, 7, 10)) {
    same(kerf_plan(r, kerf_n_tiles(n)),
         suppressMessages(geotargets::tile_n(r, n)))
  }
})

test_that("a terra source carries its block size", {
  skip_if_no_terra()
  f <- make_tif(40L, 30L, block = 16L)
  g <- kerf_raster(terra::rast(f))
  expect_equal(g$blocksize, c(16L, 16L))
  expect_equal(g$dimension, c(40L, 30L))
  expect_equal(g$dsn, f)
})
