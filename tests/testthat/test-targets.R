test_that("a plan drives dynamic branching in targets", {
  skip_if_no_terra()
  skip_if_not_installed("targets")
  skip_on_cran()
  f <- make_tif(75L, 53L, block = 16L)
  want <- terra::values(terra::focal(terra::rast(f), matrix(1, 3, 3),
                                     fun = "mean"))
  targets::tar_dir({
    targets::tar_script({
      list(
        targets::tar_target(src, Sys.getenv("KERF_TEST_TIF")),
        targets::tar_target(plan, kerf::kerf_plan(src, kerf::kerf_blocks(2, 2),
                                                  context = kerf::kerf_halo(1))),
        targets::tar_target(
          part,
          terra::wrap(terra::focal(kerf::kerf_read(src, plan),
                                   matrix(1, 3, 3), fun = "mean")),
          pattern = map(plan), iteration = "list"
        ),
        targets::tar_target(
          result,
          terra::values(kerf::kerf_reduce(lapply(part, terra::unwrap), plan))
        )
      )
    })
    withr::with_envvar(c(KERF_TEST_TIF = f), {
      targets::tar_make(callr_function = NULL, reporter = "silent")
    })
    expect_equal(length(targets::tar_read(part)), 6L)
    expect_equal(targets::tar_read(result), want)
  })
})
