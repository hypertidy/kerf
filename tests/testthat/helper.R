## a small tiled GeoTIFF with 16 x 16 blocks and non-divisible dimensions
make_tif <- function(ncol = 75L, nrow = 53L, block = 16L, nlyr = 1L) {
  r <- terra::rast(ncols = ncol, nrows = nrow, nlyrs = nlyr,
                   xmin = 100, xmax = 100 + ncol * 30,
                   ymin = -5000, ymax = -5000 + nrow * 30,
                   crs = "EPSG:3031")
  set.seed(1)
  terra::values(r) <- stats::runif(terra::ncell(r) * nlyr)
  f <- tempfile(fileext = ".tif")
  terra::writeRaster(r, f, gdal = c("TILED=YES",
                                    sprintf("BLOCKXSIZE=%i", block),
                                    sprintf("BLOCKYSIZE=%i", block)))
  f
}

skip_if_no_terra <- function() {
  testthat::skip_if_not_installed("terra", minimum_version = "1.8-10")
}
