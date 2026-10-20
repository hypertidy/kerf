test_that("extents come back in terra order", {
  g <- kerf_grid_spec(c(10, 4), c(100, 200, -20, 20))
  p <- kerf_plan(g, kerf_grid(2, 1), context = kerf_halo(1))
  e <- kerf_extents(p)
  expect_length(e, 2L)
  expect_equal(names(e[[1L]]), c("xmin", "xmax", "ymin", "ymax"))
  expect_equal(unname(e[[1L]]), c(100, 150, -20, 20))
  ec <- kerf_extents(p, "context")
  expect_equal(unname(ec[[1L]]), c(100, 160, -20, 20))
  expect_error(kerf_extents(p[, c("batch", "xoff")]), "no owned extent")
})

test_that("GDAL arguments follow the selectors present", {
  g <- kerf_grid_spec(c(100, 80), c(0, 1000, 0, 800))
  p <- kerf_plan(g, kerf_grid(2, 2), context = kerf_halo(4))
  expect_equal(kerf_gdal_args(p[2, ]), c("--window", "50,0,50,40"))
  expect_equal(kerf_gdal_args(p[2, ], "context"), c("--window", "46,0,54,44"))
  ext_only <- data.frame(batch = "a", xmin = 0, ymin = 0.5, xmax = 10,
                         ymax = 20, where = "kind = 'road'")
  expect_equal(kerf_gdal_args(ext_only),
               c("--bbox", "0,0.5,10,20", "--where", "kind = 'road'"))
  expect_error(kerf_gdal_args(p), "exactly one row")
  expect_equal(kerf_gdal_args(data.frame(batch = "a", fid_lo = 1, fid_hi = 2)),
               c("--where", "FID >= 1 AND FID <= 2"))
  expect_error(kerf_gdal_args(data.frame(batch = "a", n_owned = 2)),
               "no window")
})
