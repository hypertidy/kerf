pts_index <- function(n = 500L, seed = 2L) {
  set.seed(seed)
  x <- stats::runif(n)
  y <- stats::runif(n)
  data.frame(fid = seq_len(n), xmin = x, ymin = y, xmax = x, ymax = y)
}

test_that("hilbert_index is a bijection with unit steps", {
  expect_equal(kerf:::hilbert_index(c(0, 0, 1, 1), c(0, 1, 1, 0), 1),
               c(0, 1, 2, 3))
  for (level in 1:5) {
    n <- 2^level
    g <- expand.grid(x = 0:(n - 1), y = 0:(n - 1))
    d <- kerf:::hilbert_index(g$x, g$y, level)
    expect_setequal(d, 0:(n * n - 1))
    o <- order(d)
    step <- abs(diff(g$x[o])) + abs(diff(g$y[o]))
    expect_true(all(step == 1))
  }
})

test_that("a data frame is a feature index source", {
  idx <- pts_index(10)
  v <- kerf_vector(idx)
  expect_s3_class(v, "kerf_vector")
  expect_equal(v$n, 10)
  expect_error(kerf_vector(idx[, -1]), "needs columns: fid")
  idx$fid[2] <- 1
  expect_error(kerf_vector(idx), "unique")
  expect_error(kerf_vector_index(pts_index(5), fields = "grp"), "no field")
})

test_that("chunks are FID ranges that cover every feature once", {
  idx <- pts_index(53)
  p <- kerf_plan(idx, kerf_chunk(10))
  expect_equal(nrow(p), 6L)
  expect_equal(p$n_owned, c(10, 10, 10, 10, 10, 3))
  expect_equal(p$fid_lo, c(1, 11, 21, 31, 41, 51))
  expect_equal(p$batch[1], "fid1_10")
  expect_silent(kerf_check(p, idx))
  bad <- p
  bad$fid_hi[1] <- 12
  expect_error(kerf_check(bad), "overlap")
  bad <- p[-2, ]
  expect_error(kerf_check(bad, idx), "not owned by any batch")
})

test_that("Hilbert chunks partition the features and stay compact", {
  idx <- pts_index(2000)
  p <- kerf_plan(idx, kerf_hilbert_chunk(100))
  expect_silent(kerf_check(p, idx))
  expect_equal(sum(p$n_owned), 2000)
  expect_true(all(p$n_owned == 100))
  expect_false(is.unsorted(p$cell_lo))
  expect_true(all(utils::head(p$cell_hi, -1) < utils::tail(p$cell_lo, -1)))
  expect_equal(unique(p$rule), "point_in")

  box_area <- function(members) {
    sum(vapply(members, function(f) {
      i <- match(f, idx$fid)
      diff(range(idx$xmin[i])) * diff(range(idx$ymin[i]))
    }, numeric(1)))
  }
  pc <- kerf_plan(idx, kerf_chunk(100))
  chunk_members <- Map(function(a, b) a:b, pc$fid_lo, pc$fid_hi)
  expect_lt(box_area(p$fids), 0.3 * box_area(chunk_members))
})

test_that("Hilbert cuts move past tied keys", {
  idx <- pts_index(30)
  idx[1:12, c("xmin", "xmax")] <- 0.5
  idx[1:12, c("ymin", "ymax")] <- 0.5
  p <- kerf_plan(idx, kerf_hilbert_chunk(5))
  expect_silent(kerf_check(p, idx))
  expect_true(all(utils::head(p$cell_hi, -1) < utils::tail(p$cell_lo, -1)))
  tied <- vapply(p$fids, function(f) all(1:12 %in% f), logical(1))
  expect_equal(sum(tied), 1L)
})

test_that("features without geometry get their own final batch", {
  idx <- pts_index(20)
  idx[c(3, 7), c("xmin", "ymin", "xmax", "ymax")] <- NA
  p <- kerf_plan(idx, kerf_hilbert_chunk(6), context = kerf_buffer(0.1))
  expect_equal(p$batch[nrow(p)], "h16_empty")
  expect_equal(p$fids[[nrow(p)]], c(3, 7))
  expect_true(is.na(p$ctx_xmin[nrow(p)]))
  expect_silent(kerf_check(p, idx))
})

test_that("a fixed key extent fixes each feature's key", {
  idx <- pts_index(200)
  k1 <- kerf:::centre_keys(idx, c(0, 1, 0, 1), 16)
  more <- rbind(idx, data.frame(fid = 201, xmin = 0.99, ymin = 0.99,
                                xmax = 0.99, ymax = 0.99))
  k2 <- kerf:::centre_keys(more, c(0, 1, 0, 1), 16)
  expect_equal(k2[1:200], k1)
  expect_error(kerf_hilbert_chunk(10, level = 27), "at most 26")
})

test_that("attribute plans write predicates, quoting values", {
  idx <- pts_index(6)
  idx$kind <- c("road", "road", "o'brien", NA, "track", "track")
  p <- kerf_plan(idx, kerf_attribute("kind"))
  expect_setequal(p$where, c("\"kind\" IS NULL", "\"kind\" = 'o''brien'",
                             "\"kind\" = 'road'", "\"kind\" = 'track'"))
  expect_equal(sum(p$n_owned), 6)
  expect_silent(kerf_check(p))
  idx$num <- c(1, 1, 2.5, 2.5, 2.5, NA)
  q <- kerf_plan(idx, kerf_attribute("num"))
  expect_setequal(q$where, c("\"num\" = 1", "\"num\" = 2.5", "\"num\" IS NULL"))
})

test_that("buffer context boxes the owned features", {
  idx <- pts_index(300)
  p <- kerf_plan(idx, kerf_hilbert_chunk(50), context = kerf_buffer(0.05))
  for (i in seq_len(nrow(p))) {
    j <- match(p$fids[[i]], idx$fid)
    expect_equal(p$ctx_xmin[i], min(idx$xmin[j]) - 0.05)
    expect_equal(p$ctx_ymax[i], max(idx$ymax[j]) + 0.05)
  }
  expect_error(kerf_plan(idx, kerf_chunk(10), context = kerf_halo(1)),
               "kerf_buffer")
  expect_error(kerf_plan(kerf_grid_spec(c(5, 5)), kerf_grid(1, 1),
                         context = kerf_buffer(1)), "kerf_halo")
  expect_error(kerf_buffer(-1), ">= 0")
})

test_that("checks catch shared FIDs and half-missing context", {
  p <- data.frame(batch = c("a", "b"))
  p$fids <- list(c(1, 2), c(2, 3))
  expect_error(kerf_check(p), "more than one batch")
  p$fids <- list(c(1, 2), c(3, 4))
  p$ctx_xmin <- c(0, NA)
  p$ctx_ymin <- c(0, 0)
  p$ctx_xmax <- c(1, 1)
  p$ctx_ymax <- c(1, 1)
  expect_error(kerf_check(p), "partly missing")
  p$ctx_ymin[2] <- p$ctx_xmax[2] <- p$ctx_ymax[2] <- NA
  expect_silent(kerf_check(p))
})

test_that("GDAL arguments for FID sets and predicates", {
  p <- data.frame(batch = "a", where = "grp = 'a'")
  p$fids <- list(c(4, 9))
  expect_equal(kerf_gdal_args(p), c("--where", "(FID IN (4, 9)) AND (grp = 'a')"))
})

## ---- reading from a data source ------------------------------------------

test_that("the planning read gives FIDs and boxes", {
  skip_if_no_vector()
  f <- make_points(50L)
  src <- kerf_vector(f)
  expect_equal(src$n, 50)
  idx <- kerf_vector_index(f, fields = "grp")
  expect_equal(idx$fid, as.numeric(1:50))
  v <- terra::vect(f)
  expect_equal(idx$xmin, terra::crds(v)[, 1])
  expect_equal(idx$xmin, idx$xmax)
  expect_true(all(idx$grp %in% c("a", "b", "c")))
  expect_error(kerf_vector_index(f, fields = "nope"), "no field")
})

test_that("owned reads match across backends and reassemble the layer", {
  skip_if_no_vector()
  f <- make_points(300L)
  p <- kerf_plan(f, kerf_hilbert_chunk(70), context = kerf_buffer(0.02))
  kerf_check(p, f)
  for (i in seq_len(nrow(p))) {
    g <- kerf_read(f, p[i, ], role = "owned", backend = "gdalraster")
    t <- kerf_read(f, p[i, ], role = "owned")
    expect_setequal(as.numeric(g$FID), p$fids[[i]])
    expect_setequal(as.numeric(t$kerf_fid), p$fids[[i]])
    expect_true(all(g$kerf_owned))
  }
  parts <- lapply(seq_len(nrow(p)), function(i) {
    kerf_read(f, p[i, ], backend = "gdalraster")
  })
  expect_true(any(vapply(parts, function(d) any(!d$kerf_owned), logical(1))))
  all_rows <- kerf_reduce(parts, p)
  expect_equal(sort(as.numeric(all_rows$FID)), as.numeric(1:300))
  expect_false("kerf_owned" %in% names(all_rows))
})

test_that("neighbour counts per batch with a buffer equal the whole layer", {
  skip_if_no_vector()
  f <- make_points(400L, seed = 7L)
  d <- 0.06
  xy <- terra::crds(terra::vect(f))
  dist2 <- outer(xy[, 1], xy[, 1], "-")^2 + outer(xy[, 2], xy[, 2], "-")^2
  want <- rowSums(dist2 <= d^2) - 1

  count_owned <- function(v) {
    c <- terra::crds(v)
    own <- which(v$kerf_owned)
    n <- vapply(own, function(i) {
      sum((c[, 1] - c[i, 1])^2 + (c[, 2] - c[i, 2])^2 <= d^2) - 1
    }, numeric(1))
    data.frame(fid = v$kerf_fid[own], n = n)
  }
  p <- kerf_plan(f, kerf_hilbert_chunk(60), context = kerf_buffer(d))
  got <- do.call(rbind, lapply(seq_len(nrow(p)), function(i) {
    count_owned(kerf_read(f, p[i, ]))
  }))
  got <- got[order(got$fid), ]
  expect_equal(got$fid, 1:400)
  expect_equal(got$n, unname(want))

  ## and without the buffer, edge features undercount
  p0 <- kerf_plan(f, kerf_hilbert_chunk(60))
  got0 <- do.call(rbind, lapply(seq_len(nrow(p0)), function(i) {
    count_owned(kerf_read(f, p0[i, ]))
  }))
  got0 <- got0[order(got0$fid), ]
  expect_true(any(got0$n < want))
})

test_that("attribute plans read their groups, with context", {
  skip_if_no_vector()
  f <- system.file("ex/lux.shp", package = "terra")
  p <- kerf_plan(f, kerf_attribute("NAME_1"), context = kerf_buffer(0))
  expect_equal(nrow(p), 3L)
  expect_equal(sum(p$n_owned), 12)
  for (i in seq_len(nrow(p))) {
    r <- kerf_read(f, p[i, ], backend = "gdalraster")
    expect_equal(sum(r$kerf_owned), p$n_owned[i])
    expect_true(all(r$NAME_1[r$kerf_owned] == p$batch[i]))
  }
  lux <- kerf_plan(f, kerf_chunk(5))
  expect_silent(kerf_check(lux, f))
  expect_equal(lux$fid_lo, c(0, 5, 10))
})

test_that("an index source cannot be read", {
  idx <- pts_index(10)
  p <- kerf_plan(idx, kerf_chunk(5))
  expect_error(kerf_read(idx, p[1, ]), "holds no geometry")
})
