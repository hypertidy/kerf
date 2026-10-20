# kerf 0.0.0.9000

* Vector plans. `kerf_vector()` and `kerf_vector_index()` read FIDs and
  bounding boxes only; a data frame of `fid, xmin, ymin, xmax, ymax` is also
  a source. Generators `kerf_chunk()`, `kerf_hilbert_chunk()` (bounding-box
  centres along a Hilbert curve) and `kerf_attribute()`; `kerf_buffer()`
  context. `kerf_read()`, `kerf_trim()`, `kerf_reduce()` and
  `kerf_gdal_args()` handle vector plans with terra and gdalraster.
* The contract gains the `fids` list-column group, and a row may leave a
  group entirely missing when another group selects for it.
* First milestone: raster plans. `kerf_plan()` with `kerf_blocks()`,
  `kerf_tile_size()`, `kerf_grid()`, `kerf_n_tiles()` and `kerf_halo()`;
  the plan contract in `kerf_contract()` and `kerf_check()`; `kerf_read()`,
  `kerf_trim()` and `kerf_reduce()` for terra and gdalraster; `kerf_extents()`
  and `kerf_gdal_args()` for handing rows to geotargets and the GDAL CLI.
