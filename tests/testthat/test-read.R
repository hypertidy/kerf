test_that("a focal result per batch, halo trimmed, equals the whole grid", {
  skip_if_no_terra()
  f <- make_tif(75L, 53L, block = 16L)
  r <- terra::rast(f)
  w <- matrix(1, 5, 5)
  whole <- terra::focal(r, w, fun = "mean")

  p <- kerf_plan(f, kerf_blocks(), context = kerf_halo(2))
  parts <- lapply(seq_len(nrow(p)), function(i) {
    terra::focal(kerf_read(f, p[i, ]), w, fun = "mean")
  })
  merged <- kerf_reduce(parts, p)
  expect_equal(as.vector(terra::ext(merged)), as.vector(terra::ext(r)))
  expect_equal(terra::values(merged), terra::values(whole))
})

test_that("without a halo the batch edges differ, as they should", {
  skip_if_no_terra()
  f <- make_tif(40L, 30L, block = 16L)
  w <- matrix(1, 3, 3)
  whole <- terra::focal(terra::rast(f), w, fun = "mean")
  p <- kerf_plan(f, kerf_blocks())
  parts <- lapply(seq_len(nrow(p)), function(i) {
    terra::focal(kerf_read(f, p[i, ]), w, fun = "mean")
  })
  merged <- kerf_reduce(parts, p)
  expect_false(isTRUE(all.equal(terra::values(merged),
                                terra::values(whole))))
})

test_that("owned reads ignore the context", {
  skip_if_no_terra()
  f <- make_tif(40L, 30L, block = 16L)
  p <- kerf_plan(f, kerf_grid(2, 2), context = kerf_halo(3))
  own <- kerf_read(f, p[4, ], role = "owned")
  ctx <- kerf_read(f, p[4, ])
  expect_equal(dim(own)[1:2], c(p$ysize[4], p$xsize[4]))
  expect_equal(dim(ctx)[1:2], c(p$ctx_ysize[4], p$ctx_xsize[4]))
  expect_equal(terra::values(kerf_trim(ctx, p[4, ])), terra::values(own))
})

test_that("gdalraster arrays reassemble to the full grid", {
  skip_if_no_terra()
  skip_if_not_installed("gdalraster")
  f <- make_tif(75L, 53L, block = 16L, nlyr = 2L)
  p <- kerf_plan(f, kerf_grid(3, 2), context = kerf_halo(c(2, 1)))
  parts <- lapply(seq_len(nrow(p)), function(i) {
    kerf_read(f, p[i, ], backend = "gdalraster")
  })
  expect_equal(attr(parts[[1L]], "kerf_window"),
               unlist(p[1L, c("ctx_xoff", "ctx_yoff", "ctx_xsize",
                              "ctx_ysize")], use.names = FALSE))
  full <- kerf_reduce(parts, p)
  r <- terra::rast(f)
  want <- terra::as.array(r)
  expect_equal(dim(full), dim(want))
  expect_equal(full, want, ignore_attr = TRUE)
  expect_equal(attr(full, "kerf_window"), c(0L, 0L, 75L, 53L))
})

test_that("kerf_raster reads metadata through gdalraster", {
  skip_if_no_terra()
  skip_if_not_installed("gdalraster")
  f <- make_tif(40L, 30L, block = 16L)
  g <- kerf_raster(f)
  gt <- kerf_raster(terra::rast(f))
  expect_equal(g$dimension, gt$dimension)
  expect_equal(g$extent, gt$extent)
  expect_equal(g$blocksize, c(16L, 16L))
})

test_that("in-memory rasters read with terra only", {
  skip_if_no_terra()
  r <- terra::rast(ncols = 10, nrows = 10, vals = 1:100)
  p <- kerf_plan(r, kerf_grid(2, 2))
  expect_equal(terra::values(kerf_read(r, p[1, ]))[1:5], 1:5)
  expect_error(kerf_read(r, p[1, ], backend = "gdalraster"), "in memory")
  expect_error(kerf_read(r, p), "exactly one row")
})
