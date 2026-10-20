## column groups of the plan contract; the context groups repeat the
## selector groups with a "ctx_" prefix
selector_groups <- list(
  window = c("xoff", "yoff", "xsize", "ysize"),
  extent = c("xmin", "ymin", "xmax", "ymax"),
  fid = c("fid_lo", "fid_hi"),
  where = "where",
  cells = c("cell_lo", "cell_hi", "cell_level")
)

valid_rules <- c("centre", "point_in", "intersects", "clip")

#' The plan column contract
#'
#' A plan is any data frame that meets this contract. Every selector group
#' present in a row is a filter, and a row owns what all of its groups select
#' together. A group is either absent or complete: all of its columns, with
#' no missing values. Context groups repeat the selector groups with a
#' `ctx_` prefix; when no context group is present, context equals owned.
#'
#' @return a data frame describing each column: its group, type, whether it
#'   is required, and its meaning.
#' @export
#' @examples
#' kerf_contract()
kerf_contract <- function() {
  rows <- list(
    c("identity", "batch", "character", "yes",
      "unique batch key, stable across reruns"),
    c("window", "xoff, yoff, xsize, ysize", "whole number", "one selector",
      "owned raster window, 0-based pixel offsets and sizes"),
    c("extent", "xmin, ymin, xmax, ymax", "numeric", "one selector",
      "owned extent in the source CRS, half-open on the max edges"),
    c("fid", "fid_lo, fid_hi", "numeric", "one selector",
      "owned feature id range, inclusive"),
    c("where", "where", "character", "one selector",
      "owned attribute predicate, an OGR SQL fragment"),
    c("cells", "cell_lo, cell_hi, cell_level", "numeric", "one selector",
      "owned interval of finest-level cell indices"),
    c("context", "ctx_ + any selector group", "as above", "no",
      "read-only context; must contain the owned selection"),
    c("rule", "rule", "character", "no",
      "ownership rule: centre, point_in, intersects or clip"),
    c("source", "dsn", "character", "no",
      "per-row data source, overriding the source handle"),
    c("stats", "n_owned, est_bytes", "numeric", "no",
      "owned member count and estimated bytes; filled in by kerf_check()")
  )
  out <- as.data.frame(do.call(rbind, rows), stringsAsFactors = FALSE)
  names(out) <- c("group", "columns", "type", "required", "meaning")
  out
}

#' Check a plan against the contract
#'
#' Validates a plan, generated or hand-built. Without a source, the checks
#' are structural: column groups complete, types and values valid, `batch`
#' unique, context containing owned. With a raster source, owned windows must
#' also lie on the grid and partition it exactly (no overlaps, no gaps), and
#' context windows must lie on the grid.
#'
#' @param plan a data frame.
#' @param x optional source, anything [kerf_raster()] accepts.
#'
#' @return the plan, invisibly, with `n_owned` filled in when it was missing
#'   and can be computed. Problems are reported together in one error.
#' @export
#' @examples
#' g <- kerf_grid_spec(c(10, 10))
#' p <- data.frame(batch = c("left", "right"), xoff = c(0L, 5L), yoff = 0L,
#'                 xsize = 5L, ysize = 10L)
#' kerf_check(p, g)
#'
#' ## an overlapping plan fails
#' p$xsize <- c(6L, 5L)
#' try(kerf_check(p, g))
kerf_check <- function(plan, x = NULL) {
  if (!is.data.frame(plan)) stop("'plan' must be a data frame", call. = FALSE)
  if (nrow(plan) < 1L) stop("'plan' has no rows", call. = FALSE)
  problems <- character()
  add <- function(...) problems <<- c(problems, sprintf(...))

  ## identity
  if (!"batch" %in% names(plan)) {
    add("missing the required column 'batch'")
  } else if (!is.character(plan$batch)) {
    add("'batch' must be character")
  } else {
    if (anyNA(plan$batch)) add("'batch' has missing values")
    dup <- unique(plan$batch[duplicated(plan$batch)])
    if (length(dup)) {
      add("'batch' is not unique: %s", paste(utils::head(dup, 5L),
                                             collapse = ", "))
    }
  }

  ## selector groups, owned and context
  present <- list()
  for (prefix in c("", "ctx_")) {
    for (g in names(selector_groups)) {
      cols <- paste0(prefix, selector_groups[[g]])
      has <- cols %in% names(plan)
      if (any(has) && !all(has)) {
        add("group '%s%s' is incomplete, missing: %s", prefix, g,
            paste(cols[!has], collapse = ", "))
      } else if (all(has)) {
        na <- vapply(plan[cols], anyNA, logical(1))
        if (any(na)) {
          add("group '%s%s' has missing values in: %s", prefix, g,
              paste(cols[na], collapse = ", "))
        }
        present[[paste0(prefix, g)]] <- cols
      }
    }
  }
  owned <- intersect(names(present), names(selector_groups))
  if (!length(owned)) {
    add("no owned selector group; need one of: %s",
        paste(names(selector_groups), collapse = ", "))
  }

  for (prefix in c("", "ctx_")) {
    if (paste0(prefix, "window") %in% names(present)) {
      problems <- c(problems, check_window(plan, prefix))
    }
    if (paste0(prefix, "extent") %in% names(present)) {
      problems <- c(problems, check_extent(plan, prefix))
    }
  }

  if ("rule" %in% names(plan)) {
    bad <- setdiff(unique(plan$rule), valid_rules)
    if (length(bad)) {
      add("unknown 'rule' value(s): %s; expected one of: %s",
          paste(bad, collapse = ", "), paste(valid_rules, collapse = ", "))
    }
  }

  ## context contains owned
  if (all(c("window", "ctx_window") %in% names(present)) &&
      !length(problems)) {
    inside <- plan$ctx_xoff <= plan$xoff & plan$ctx_yoff <= plan$yoff &
      plan$ctx_xoff + plan$ctx_xsize >= plan$xoff + plan$xsize &
      plan$ctx_yoff + plan$ctx_ysize >= plan$yoff + plan$ysize
    if (!all(inside)) {
      add("context window does not contain the owned window for batch: %s",
          paste(utils::head(plan$batch[!inside], 5L), collapse = ", "))
    }
  }

  ## against a raster source
  if (!is.null(x) && !length(problems)) {
    src <- kerf_raster(x)
    for (prefix in c("", "ctx_")) {
      if (paste0(prefix, "window") %in% names(present)) {
        problems <- c(problems, check_bounds(plan, prefix, src$dimension))
      }
    }
    if ("window" %in% names(present) && !length(problems)) {
      problems <- c(problems, check_partition(plan, src$dimension))
    }
  }

  if (length(problems)) {
    stop(paste0("invalid plan:\n", paste0("* ", problems, collapse = "\n")),
         call. = FALSE)
  }

  if (!"n_owned" %in% names(plan) && "window" %in% names(present)) {
    plan$n_owned <- as.numeric(plan$xsize) * plan$ysize
  }
  invisible(plan)
}

## ---- internals -------------------------------------------------------------

is_whole <- function(v) is.numeric(v) && all(v == round(v))

check_window <- function(plan, prefix) {
  cols <- paste0(prefix, selector_groups$window)
  out <- character()
  whole <- vapply(plan[cols], is_whole, logical(1))
  if (!all(whole)) {
    return(sprintf("window columns must be whole numbers: %s",
                   paste(cols[!whole], collapse = ", ")))
  }
  if (any(plan[[cols[1L]]] < 0) || any(plan[[cols[2L]]] < 0)) {
    out <- c(out, sprintf("'%s' and '%s' must be >= 0", cols[1L], cols[2L]))
  }
  if (any(plan[[cols[3L]]] < 1) || any(plan[[cols[4L]]] < 1)) {
    out <- c(out, sprintf("'%s' and '%s' must be >= 1", cols[3L], cols[4L]))
  }
  out
}

check_extent <- function(plan, prefix) {
  cols <- paste0(prefix, selector_groups$extent)
  num <- vapply(plan[cols], is.numeric, logical(1))
  if (!all(num)) {
    return(sprintf("extent columns must be numeric: %s",
                   paste(cols[!num], collapse = ", ")))
  }
  bad <- plan[[cols[3L]]] <= plan[[cols[1L]]] |
    plan[[cols[4L]]] <= plan[[cols[2L]]]
  if (any(bad)) {
    return(sprintf("%sxmax/%symax must exceed %sxmin/%symin for batch: %s",
                   prefix, prefix, prefix, prefix,
                   paste(utils::head(plan$batch[bad], 5L), collapse = ", ")))
  }
  character()
}

check_bounds <- function(plan, prefix, dm) {
  xo <- plan[[paste0(prefix, "xoff")]]
  yo <- plan[[paste0(prefix, "yoff")]]
  xs <- plan[[paste0(prefix, "xsize")]]
  ys <- plan[[paste0(prefix, "ysize")]]
  bad <- xo + xs > dm[1L] | yo + ys > dm[2L]
  if (any(bad)) {
    return(sprintf("%swindow extends past the %i x %i grid for batch: %s",
                   prefix, dm[1L], dm[2L],
                   paste(utils::head(plan$batch[bad], 5L), collapse = ", ")))
  }
  character()
}

## owned windows must tile the grid exactly: equal total area and no overlap
## (with every window in bounds, that is an exact partition)
check_partition <- function(plan, dm) {
  xo <- plan$xoff
  yo <- plan$yoff
  xs <- plan$xsize
  ys <- plan$ysize
  out <- character()
  area <- sum(as.numeric(xs) * ys)
  total <- as.numeric(dm[1L]) * dm[2L]

  ux <- sort(unique(c(xo, xo + xs)))
  uy <- sort(unique(c(yo, yo + ys)))
  if ((length(ux) - 1) * (length(uy) - 1) <= 4e6) {
    count <- matrix(0L, length(uy) - 1L, length(ux) - 1L)
    ix0 <- match(xo, ux)
    ix1 <- match(xo + xs, ux) - 1L
    iy0 <- match(yo, uy)
    iy1 <- match(yo + ys, uy) - 1L
    for (i in seq_along(xo)) {
      ii <- iy0[i]:iy1[i]
      jj <- ix0[i]:ix1[i]
      count[ii, jj] <- count[ii, jj] + 1L
    }
    if (any(count > 1L)) {
      hit <- vapply(seq_along(xo), function(i) {
        any(count[iy0[i]:iy1[i], ix0[i]:ix1[i]] > 1L)
      }, logical(1))
      out <- c(out, sprintf("owned windows overlap for batch: %s",
                            paste(utils::head(plan$batch[hit], 5L),
                                  collapse = ", ")))
    }
  } else {
    warning("too many distinct window edges to test for overlap; ",
            "checked total area only", call. = FALSE)
  }
  if (area < total && !length(out)) {
    out <- c(out, sprintf("owned windows leave %s cells of the grid uncovered",
                          format(total - area, big.mark = ",")))
  }
  if (area > total && !length(out)) {
    out <- c(out, "owned windows cover more cells than the grid holds")
  }
  out
}
