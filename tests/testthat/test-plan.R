test_that("grid spec validates its inputs", {
  expect_error(kerf_grid_spec(c(0, 10)), "positive integers")
  expect_error(kerf_grid_spec(c(10, 10), c(0, -1, 0, 1)), "xmax > xmin")
  expect_error(kerf_grid_spec(c(10, 10), blocksize = 3), "two positive")
  g <- kerf_grid_spec(c(10, 8))
  expect_equal(g$extent, c(0, 10, 0, 8))
  expect_s3_class(kerf_raster(g), "kerf_grid")
})

test_that("tile and block plans partition the grid, row-major", {
  g <- kerf_grid_spec(c(75, 53), c(0, 75, 0, 53), blocksize = c(16, 16))
  for (k in list(kerf_blocks(), kerf_blocks(2, 3), kerf_tile_size(20, 7),
                 kerf_tile_size(100), kerf_grid(4, 3), kerf_grid(1, 1),
                 kerf_n_tiles(7))) {
    p <- kerf_plan(g, k)
    expect_silent(kerf_check(p, g))
    expect_equal(sum(p$n_owned), 75 * 53)
    expect_false(is.unsorted(p$yoff))
    expect_false(anyDuplicated(p$batch) > 0)
  }
  p <- kerf_plan(g, kerf_blocks())
  expect_equal(nrow(p), 5L * 4L)
  expect_equal(p$xsize[1:5], c(16L, 16L, 16L, 16L, 11L))
  expect_equal(p$ysize[c(1, 20)], c(16L, 5L))
  expect_equal(p$batch[1:2], c("x0_y0", "x16_y0"))
})

test_that("extents follow the windows", {
  g <- kerf_grid_spec(c(10, 4), c(100, 200, -20, 20))
  p <- kerf_plan(g, kerf_grid(2, 2))
  expect_equal(p$xmin, c(100, 150, 100, 150))
  expect_equal(p$xmax, c(150, 200, 150, 200))
  expect_equal(p$ymax, c(20, 20, 0, 0))
  expect_equal(p$ymin, c(0, 0, -20, -20))
})

test_that("halo context is clipped to the grid and contains owned", {
  g <- kerf_grid_spec(c(20, 10))
  p <- kerf_plan(g, kerf_grid(2, 2), context = kerf_halo(3))
  expect_true(all(c("ctx_xoff", "ctx_ymax") %in% names(p)))
  expect_equal(p$ctx_xoff, c(0L, 7L, 0L, 7L))
  expect_equal(p$ctx_xsize, c(13L, 13L, 13L, 13L))
  expect_equal(p$ctx_yoff, c(0L, 0L, 2L, 2L))
  expect_equal(p$ctx_ysize, c(8L, 8L, 8L, 8L))
  expect_silent(kerf_check(p, g))
  p2 <- kerf_plan(g, kerf_grid(2, 2), context = kerf_halo(c(1, 0)))
  expect_equal(p2$ctx_ysize, p2$ysize)
})

test_that("blocks need a native block size", {
  expect_error(kerf_plan(kerf_grid_spec(c(10, 10)), kerf_blocks()),
               "no native block size")
})

test_that("n_tiles makes the same choice as geotargets::tile_n", {
  shape <- function(n) {
    p <- kerf_n_tiles(n)$params
    c(p$nrow, p$ncol)
  }
  expect_equal(shape(1), c(1, 1))
  expect_equal(shape(4), c(2, 2))
  expect_equal(shape(6), c(2, 3))
  expect_equal(shape(7), c(2, 4))
  expect_equal(shape(12), c(3, 4))
})

test_that("generators reject bad arguments", {
  expect_error(kerf_grid(0, 1), "positive whole")
  expect_error(kerf_tile_size(2.5), "positive whole")
  expect_error(kerf_halo(-1), ">= 0")
  expect_error(kerf_plan(kerf_grid_spec(c(5, 5)), "grid"), "key generator")
})
