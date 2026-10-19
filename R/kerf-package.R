#' kerf: plan spatial batches
#'
#' A plan is a plain data frame with one row per batch. Each row says which
#' members the batch *owns* and, optionally, which larger *context* it may
#' read. Owned sets partition the input exactly, so results reduced over
#' owned members never double count, while context (a halo of cells, for
#' example) is read-only.
#'
#' The main functions are:
#'
#' * [kerf_raster()] and [kerf_grid_spec()] describe a source from metadata.
#' * [kerf_plan()] builds a plan from a key generator ([kerf_blocks()],
#'   [kerf_grid()], [kerf_n_tiles()], [kerf_tile_size()]) and an optional
#'   context ([kerf_halo()]).
#' * [kerf_check()] validates any plan, generated or hand-built, against the
#'   column contract described by [kerf_contract()].
#' * [kerf_read()], [kerf_trim()] and [kerf_reduce()] materialise and
#'   recombine batches; [kerf_extents()] and [kerf_gdal_args()] hand a plan
#'   row to other tools.
#'
#' @keywords internal
"_PACKAGE"
