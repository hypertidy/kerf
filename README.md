# kerf

<!-- badges: start -->
[![R-CMD-check](https://github.com/hypertidy/kerf/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/hypertidy/kerf/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

kerf plans how to split spatial data into batches. A plan is a plain data
frame with one row per batch, saying which members the batch *owns* and,
optionally, which larger *context* it may read. Plans are made from metadata
only, they can be built by hand, and any backend can read them, so they fit
dynamic branching in [targets](https://github.com/ropensci/targets) and
windowed reads with [terra](https://github.com/rspatial/terra) or
[gdalraster](https://github.com/firelab/gdalraster).

A kerf is the cut a saw makes, and the width of that cut.

This is an early, experimental package. It covers raster plans (windows,
with a halo of cells as context) and vector plans (FID ranges, Hilbert-ordered
FID sets and attribute predicates, with a buffer as context).

## Owned and context

Every batch has two member sets.

* **Owned** members partition the input exactly. Results are reduced over
  owned members only.
* **Context** members are a read-only superset the computation needs, such
  as a halo of cells for a focal operation.

So a neighbourhood operation run batch by batch, each with a halo, then
trimmed to its owned window, gives exactly the whole-grid answer. The same
holds for vectors: count each owned feature's neighbours within a distance,
using a buffer of that distance as context, and the batches together give
the whole-layer counts.

## Example

```r
library(kerf)
f <- system.file("ex/elev.tif", package = "terra")

## a plan: 2 x 2 batches, each with a 2-cell halo
p <- kerf_plan(f, kerf_grid(2, 2), context = kerf_halo(2))
p[, c("batch", "xoff", "yoff", "xsize", "ysize", "ctx_xoff", "ctx_xsize")]

## any data frame meeting the contract is a plan
kerf_check(p, f)

## read, process, reduce
w <- matrix(1, 5, 5)
parts <- lapply(seq_len(nrow(p)), function(i) {
  terra::focal(kerf_read(f, p[i, ]), w, fun = "mean")
})
out <- kerf_reduce(parts, p)
```

## Key generators

For rasters:

| Generator | Batches |
| --- | --- |
| `kerf_blocks(n_row, n_col)` | multiples of the source's native block size |
| `kerf_tile_size(ncol, nrow)` | fixed pixel size, clipped at the edges |
| `kerf_grid(ncol, nrow)` | an `ncol` by `nrow` split, edges snapped to cells |
| `kerf_n_tiles(n)` | a near-square split into about `n` batches |

`kerf_blocks()`, `kerf_grid()` and `kerf_n_tiles()` reproduce the
`tile_blocksize()`, `tile_grid()` and `tile_n()` helpers in
[geotargets](https://github.com/ropensci/geotargets), and `kerf_extents()`
returns plan rows in the same form those helpers do.

For vectors:

| Generator | Batches | Owned selection |
| --- | --- | --- |
| `kerf_chunk(n)` | runs of `n` features in FID order | `fid_lo, fid_hi` |
| `kerf_hilbert_chunk(n)` | about `n` features, compact along a Hilbert curve | `fids`, with key interval `cell_lo, cell_hi` |
| `kerf_attribute(field)` | one per distinct value | `where` |

Vector planning reads only FIDs and bounding boxes (`kerf_vector_index()`).
A data frame with `fid, xmin, ymin, xmax, ymax` columns is also a source, so
a plan can be made from a feature index produced by any tool, such as a slim
`gdal vector sql` projection.

```r
f <- system.file("ex/lux.shp", package = "terra")
p <- kerf_plan(f, kerf_hilbert_chunk(4), context = kerf_buffer(0.05))
v <- kerf_read(f, p[1, ])   # owned features plus context, marked kerf_owned
kerf_trim(v, p[1, ])        # just the owned features
```

## The plan contract

`kerf_contract()` prints the full contract. In short: `batch` is required,
and each row needs at least one owned selector group. Every selector group
present in a row is a filter, and the row owns what all of its groups
select together. A row may leave a group entirely missing, and then that
group does not constrain it.

| Group | Columns |
| --- | --- |
| identity | `batch` |
| window | `xoff, yoff, xsize, ysize` (0-based pixels, as GDAL uses) |
| extent | `xmin, ymin, xmax, ymax` (half-open on the max edges) |
| fid | `fid_lo, fid_hi` |
| fids | `fids` (a list column of FID vectors) |
| where | `where` (an OGR SQL fragment) |
| cells | `cell_lo, cell_hi, cell_level` |
| context | `ctx_` plus any selector group |
| optional | `rule`, `dsn`, `n_owned`, `est_bytes` |

## In a targets pipeline

A plan is an ordinary data frame target, and `pattern = map(plan)` branches
over its rows.

```r
library(targets)
list(
  tar_target(src, "big.tif", format = "file"),
  tar_target(plan, kerf::kerf_plan(src, kerf::kerf_blocks(4, 4),
                                   context = kerf::kerf_halo(16))),
  tar_target(part, terra::wrap(process(kerf::kerf_read(src, plan))),
             pattern = map(plan), iteration = "list"),
  tar_target(result, kerf::kerf_reduce(lapply(part, terra::unwrap), plan))
)
```

With `backend = "gdalraster"`, `kerf_read()` returns plain arrays for rasters
and plain data frames (geometry as WKB) for vectors, which targets stores
without any wrapping.

## Handing rows to GDAL

`kerf_gdal_args()` turns a plan row into arguments for the GDAL command line
interface, such as `--window` for `gdal raster clip`, or `--bbox` and
`--where` for `gdal vector filter`. FID ranges and sets become `--where`
clauses (`FID >= ... AND FID <= ...`, `FID IN (...)`).
