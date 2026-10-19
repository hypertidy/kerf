g <- kerf_grid_spec(c(10, 10))
hand <- function() {
  data.frame(batch = c("left", "right"), xoff = c(0L, 5L), yoff = 0L,
             xsize = 5L, ysize = 10L, stringsAsFactors = FALSE)
}

test_that("a hand-built plan passes and gains n_owned", {
  p <- kerf_check(hand(), g)
  expect_equal(p$n_owned, c(50, 50))
})

test_that("structural problems are reported together", {
  p <- hand()
  p$batch <- c("a", "a")
  p$xoff <- NULL
  err <- tryCatch(kerf_check(p), error = conditionMessage)
  expect_match(err, "not unique")
  expect_match(err, "group 'window' is incomplete")
})

test_that("overlaps, gaps and out-of-bounds windows fail against a source", {
  p <- hand()
  p$xsize <- c(6L, 5L)
  expect_error(kerf_check(p, g), "overlap")
  p <- hand()
  p$xsize <- c(4L, 5L)
  expect_error(kerf_check(p, g), "uncovered")
  p <- hand()
  p$xoff <- c(0L, 6L)
  expect_error(kerf_check(p, g), "extends past")
})

test_that("context must contain owned", {
  p <- hand()
  p$ctx_xoff <- c(0L, 6L)
  p$ctx_yoff <- 0L
  p$ctx_xsize <- c(6L, 4L)
  p$ctx_ysize <- 10L
  expect_error(kerf_check(p), "does not contain")
})

test_that("rules, types and selectors are validated", {
  p <- hand()
  p$rule <- "nearby"
  expect_error(kerf_check(p), "unknown 'rule'")
  expect_error(kerf_check(data.frame(batch = "a")), "no owned selector")
  p <- hand()
  p$xoff <- c(0.5, 5)
  expect_error(kerf_check(p), "whole numbers")
  p <- data.frame(batch = "a", xmin = 1, ymin = 0, xmax = 0, ymax = 1)
  expect_error(kerf_check(p), "must exceed")
})

test_that("non-raster selectors pass structural checks", {
  p <- data.frame(batch = c("a", "b"), fid_lo = c(1, 101), fid_hi = c(100, 200),
                  rule = "point_in")
  expect_silent(kerf_check(p))
  p <- data.frame(batch = c("east", "west"),
                  where = c("region = 'east'", "region = 'west'"))
  expect_silent(kerf_check(p))
})

test_that("the contract describes every group", {
  k <- kerf_contract()
  expect_true(all(c("identity", "window", "extent", "fid", "where", "cells",
                    "context", "rule", "source", "stats") %in% k$group))
})
