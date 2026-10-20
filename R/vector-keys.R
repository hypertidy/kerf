#' Vector key generators
#'
#' Key generators for vector sources, resolved by [kerf_plan()] against the
#' feature index from [kerf_vector_index()]. Every feature is owned by
#' exactly one batch.
#'
#' * `kerf_chunk()` cuts the features, in FID order, into runs of `n`. Each
#'   batch owns an FID range (`fid_lo, fid_hi`), which GDAL can filter on
#'   directly. Nothing spatial is implied.
#' * `kerf_hilbert_chunk()` orders features along a Hilbert curve by the
#'   centre of their bounding box and cuts that order into runs of about
#'   `n`, giving spatially compact batches without reading full geometry.
#'   Each batch owns an explicit FID set (the list column `fids`) and records
#'   its interval of Hilbert keys (`cell_lo, cell_hi, cell_level`). Cuts are
#'   moved to the next change of key, so key intervals never overlap; batch
#'   sizes can differ from `n` by the number of features that share a key.
#'   Centres are taken in the source CRS, so a projected CRS gives the most
#'   compact batches.
#' * `kerf_attribute()` makes one batch per distinct value of a field. Each
#'   batch owns a `where` predicate (missing values become `IS NULL`).
#'
#' Features with no geometry have no bounding-box centre; `kerf_hilbert_chunk()`
#' puts them together in a final batch.
#'
#' @param n number of features per batch.
#' @param level Hilbert curve level, 1 to 26: the key grid is `2^level`
#'   cells on each side.
#' @param extent the key extent `c(xmin, xmax, ymin, ymax)`. Defaults to the
#'   extent of the features. Fixing it keeps every feature's key the same when
#'   features are added or removed elsewhere. Batches that come before the
#'   change along the curve then keep their key intervals and names; cuts
#'   after it can move, because batches are cut by count.
#' @param field the attribute field to group by.
#'
#' @return a `"kerf_keys"` object.
#' @name kerf_vector_keys
#' @examples
#' set.seed(1)
#' x <- runif(500)
#' y <- runif(500)
#' idx <- data.frame(fid = 1:500, xmin = x, ymin = y, xmax = x, ymax = y)
#' p <- kerf_plan(idx, kerf_hilbert_chunk(100))
#' p[, c("batch", "cell_lo", "cell_hi", "n_owned")]
NULL

#' @rdname kerf_vector_keys
#' @export
kerf_chunk <- function(n) {
  new_keys("chunk", n = check_count(n, "n"))
}

#' @rdname kerf_vector_keys
#' @export
kerf_hilbert_chunk <- function(n, level = 16, extent = NULL) {
  level <- check_count(level, "level")
  if (level > 26L) {
    stop("'level' must be at most 26, so keys stay exact as doubles",
         call. = FALSE)
  }
  if (!is.null(extent)) {
    extent <- as.numeric(extent)
    if (length(extent) != 4L || anyNA(extent) ||
        extent[2L] < extent[1L] || extent[4L] < extent[3L]) {
      stop("'extent' must be c(xmin, xmax, ymin, ymax)", call. = FALSE)
    }
  }
  new_keys("hilbert_chunk", n = check_count(n, "n"), level = level,
           extent = extent)
}

#' @rdname kerf_vector_keys
#' @export
kerf_attribute <- function(field) {
  if (!is.character(field) || length(field) != 1L || is.na(field) ||
      !nzchar(field)) {
    stop("'field' must be a single field name", call. = FALSE)
  }
  new_keys("attribute", field = field)
}

#' Vector context: a buffer around each batch
#'
#' The context of a vector batch is every feature whose bounding box meets the
#' box around the batch's owned features, widened by `distance`. It is read
#' through a spatial filter and marked so [kerf_trim()] can drop it. Use it for
#' neighbourhood work, such as counting features within `distance` of each
#' owned feature.
#'
#' @param distance buffer distance in source CRS units, `>= 0`.
#'
#' @return a `"kerf_context"` object.
#' @export
#' @examples
#' idx <- data.frame(fid = 1:4, xmin = c(0, 1, 5, 6), ymin = 0,
#'                   xmax = c(0, 1, 5, 6), ymax = 0)
#' kerf_plan(idx, kerf_chunk(2), context = kerf_buffer(0.5))
kerf_buffer <- function(distance) {
  distance <- as.numeric(distance)
  if (length(distance) != 1L || is.na(distance) || distance < 0) {
    stop("'distance' must be a single number >= 0", call. = FALSE)
  }
  new_context("buffer", distance = distance)
}

## ---- internals -------------------------------------------------------------

vector_key_types <- c("chunk", "hilbert_chunk", "attribute")

vector_plan <- function(src, keys, context) {
  if (!is.null(context) && context$type != "buffer") {
    stop("vector plans take kerf_buffer() context", call. = FALSE)
  }
  p <- keys$params
  fields <- if (keys$type == "attribute") p$field else NULL
  idx <- kerf_vector_index(src, fields = fields)
  if (!nrow(idx)) stop("the source has no features", call. = FALSE)

  plan <- switch(keys$type,
    chunk = plan_chunk(idx, p$n),
    hilbert_chunk = plan_hilbert(idx, p$n, p$level,
                                 p$extent %||% src_extent(idx)),
    attribute = plan_attribute(idx, p$field)
  )
  if (!is.null(context)) {
    plan <- cbind(plan, buffer_context(idx, plan$.members,
                                       context$params$distance))
  }
  plan$.members <- NULL
  rownames(plan) <- NULL
  plan
}

`%||%` <- function(a, b) if (is.null(a)) b else a

src_extent <- function(idx) {
  ok <- !is.na(idx$xmin)
  if (!any(ok)) return(c(0, 0, 0, 0))
  c(min(idx$xmin[ok]), max(idx$xmax[ok]), min(idx$ymin[ok]),
    max(idx$ymax[ok]))
}

## each plan builder returns a data frame with a temporary list column
## `.members` of row indices into idx, used for context and then dropped

plan_chunk <- function(idx, n) {
  grp <- (seq_len(nrow(idx)) - 1L) %/% n
  members <- unname(split(seq_len(nrow(idx)), grp))
  lo <- vapply(members, function(i) idx$fid[i[1L]], numeric(1))
  hi <- vapply(members, function(i) idx$fid[i[length(i)]], numeric(1))
  out <- data.frame(batch = sprintf("fid%s_%s", fmt_num(lo), fmt_num(hi)),
                    fid_lo = lo, fid_hi = hi,
                    n_owned = lengths(members), stringsAsFactors = FALSE)
  out$.members <- members
  out
}

plan_hilbert <- function(idx, n, level, extent) {
  key <- centre_keys(idx, extent, level)
  ord <- order(key, idx$fid, na.last = TRUE)
  k <- key[ord]
  has <- !is.na(k)
  kk <- k[has]
  m <- length(kk)
  ## cut every n, then move each cut forward to the next change of key
  starts <- 1L
  while (m > 0L) {
    cut <- starts[length(starts)] + n
    if (cut > m) break
    while (cut <= m && kk[cut] == kk[cut - 1L]) cut <- cut + 1L
    if (cut > m) break
    starts <- c(starts, cut)
  }
  members <- list()
  if (m > 0L) {
    ends <- c(starts[-1L] - 1L, m)
    members <- Map(function(a, b) ord[has][a:b], starts, ends)
  }
  lo <- vapply(members, function(i) min(key[i]), numeric(1))
  hi <- vapply(members, function(i) max(key[i]), numeric(1))
  batch <- sprintf("h%i_%s_%s", level, fmt_num(lo), fmt_num(hi))
  if (any(!has)) {
    members <- c(members, list(ord[!has]))
    lo <- c(lo, NA_real_)
    hi <- c(hi, NA_real_)
    batch <- c(batch, sprintf("h%i_empty", level))
  }
  out <- data.frame(batch = batch, stringsAsFactors = FALSE)
  out$fids <- lapply(members, function(i) sort(idx$fid[i]))
  out$cell_lo <- lo
  out$cell_hi <- hi
  out$cell_level <- ifelse(is.na(lo), NA_integer_, level)
  out$rule <- "point_in"
  out$n_owned <- lengths(members)
  out$.members <- members
  out
}

plan_attribute <- function(idx, field) {
  v <- idx[[field]]
  if (is.factor(v)) v <- as.character(v)
  if (inherits(v, "integer64")) v <- as.numeric(v)
  key <- ifelse(is.na(v), "\001NA", as.character(v))
  members <- split(seq_len(nrow(idx)), key)
  vals <- lapply(members, function(i) v[i[1L]])
  where <- vapply(vals, function(val) sql_equals(field, val), character(1))
  label <- vapply(vals, function(val) {
    if (is.na(val)) "NA" else as.character(val)
  }, character(1))
  label <- make.unique(gsub("[^A-Za-z0-9._-]+", "_", label), sep = "_")
  out <- data.frame(batch = unname(label), where = unname(where),
                    n_owned = unname(lengths(members)),
                    stringsAsFactors = FALSE)
  out$.members <- unname(members)
  out
}

## the box around each batch's owned features, widened by d
buffer_context <- function(idx, members, d) {
  box <- t(vapply(members, function(i) {
    ok <- i[!is.na(idx$xmin[i])]
    if (!length(ok)) return(rep(NA_real_, 4L))
    c(min(idx$xmin[ok]) - d, min(idx$ymin[ok]) - d,
      max(idx$xmax[ok]) + d, max(idx$ymax[ok]) + d)
  }, numeric(4)))
  data.frame(ctx_xmin = box[, 1L], ctx_ymin = box[, 2L],
             ctx_xmax = box[, 3L], ctx_ymax = box[, 4L])
}

fmt_num <- function(x) {
  ifelse(is.na(x), "NA", formatC(x, format = "f", digits = 0, big.mark = ""))
}

## an OGR SQL equality test for one value
sql_equals <- function(field, val) {
  id <- sprintf("\"%s\"", gsub("\"", "\"\"", field))
  if (length(val) != 1L || is.na(val)) return(sprintf("%s IS NULL", id))
  if (is.numeric(val)) return(sprintf("%s = %s", id, format_num_sql(val)))
  if (is.logical(val)) return(sprintf("%s = %i", id, as.integer(val)))
  sprintf("%s = '%s'", id, gsub("'", "''", as.character(val)))
}

format_num_sql <- function(x) sprintf("%.17g", x)
