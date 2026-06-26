#' Count neighbouring presence cells within a radius
#'
#' Computes the spatial corroboration variable used by DetectConf as a
#' predictor: for each target cell, the number of \emph{other} presence
#' cells within \code{radius_km}. 
#'
#' @param target_cells Integer vector of \code{cell_id} values for which
#'   to compute the count. May be the presence cells, pseudo-absence
#'   cells, or projection-region cells.
#' @param presence_cells_ref Integer vector of reference presence
#'   \code{cell_id} values to count against. Typically the global set of
#'   cleaned presence cells.
#' @param raster_template The \code{terra::SpatRaster} used to obtain
#'   cell coordinates. Must be the same raster used to assign
#'   \code{cell_id} values throughout the pipeline.
#' @param radius_km Search radius in kilometres. Default \code{50}, which
#'   roughly matches one neighbouring 5-arcminute grid cell at mid
#'   latitudes.
#' @param block_size Integer block size for memory-safe distance
#'   computation. Default \code{2000}. Larger values are faster but use
#'   more RAM.
#'
#' @return Integer vector of length \code{length(target_cells)}, giving
#'   the count of \code{presence_cells_ref} within \code{radius_km} of
#'   each target cell.
#'
#' @importFrom terra xyFromCell
#' @importFrom fields rdist.earth
#' @export
dc_compute_corroboration <- function(target_cells,
                                     presence_cells_ref,
                                     raster_template,
                                     radius_km  = 50,
                                     block_size = 2000L) {

  pres_xy   <- terra::xyFromCell(raster_template, presence_cells_ref)
  target_xy <- terra::xyFromCell(raster_template, target_cells)

  n    <- nrow(target_xy)
  corr <- integer(n)

  for (i in seq(1L, n, by = block_size)) {
    idx <- i:min(i + block_size - 1L, n)
    d_mat <- fields::rdist.earth(target_xy[idx, , drop = FALSE],
                                 pres_xy, miles = FALSE)
    corr[idx] <- rowSums(d_mat <= radius_km, na.rm = TRUE)
    rm(d_mat)
  }
  corr
}
