#' Build a national detection-confidence projection surface
#'
#' Computes all effort and
#' prior variables for every grid cell within a target country, ready
#' for prediction by \code{dc_predict}. Implements adaptive temporal
#' binning: a single GBIF download is requested per year bin, and any
#' bin whose download exceeds \code{max_size_bytes} is split at its
#' midpoint and re-submitted until all bins are within the memory
#' threshold.
#'
#' @param country_code ISO 3166-1 alpha-3 country code (e.g.
#'   \code{"BEL"} for Belgium, \code{"DEU"} for Germany).
#' @param year_bins List of initial year bins (e.g.
#'   \code{list(c(1800, 2005), c(2006, 2020), c(2021, 2026))}). The
#'   function may refine these adaptively.
#' @param raster_template Classified climate envelope raster.
#' @param presence_cells Integer vector of presence cell IDs (used for
#'   the spatial corroboration predictor).
#' @param sinas_priors List of priors from \code{dc_sinas_prior}.
#' @param focal_species Scientific name of the focal species (with
#'   space, e.g. \code{"Asclepias syriaca"}). Records matching this
#'   name are removed from effort calculations to prevent circularity.
#' @param focal_exclusion_keys Character vector of cell-observer
#'   exclusion keys (from \code{dc_focal_exclusion_keys}).
#' @param cfg Taxonomic relevance configuration (see
#'   \code{?dc_expertise_config}).
#' @param keys_file Optional file path for caching GBIF download keys
#'   across runs. If supplied, existing keys are reused and new keys
#'   are saved as downloads complete, making the function resumable.
#' @param scratch_dir Directory for temporary zip downloads. Default
#'   \code{tempdir()}.
#' @param max_size_bytes Maximum GBIF download size in bytes before a
#'   bin is split. Default \code{8e8} (~800 MB).
#' @param seven_zip_path Path to 7-Zip executable.
#' @param baseline_prior Confirmation-prior baseline for cells with no
#'   SInAS record. Default \code{0.1} (matches \code{dc_sinas_prior}).
#' @param vars_sel Character vector of predictors to validate before
#'   returning. Defaults to \code{dc_default_vars()}.
#' @param corroboration_radius_km Radius in km for spatial corroboration.
#'   Default \code{50}.
#'
#' @return A \code{data.table} with one row per country grid cell and
#'   all columns required by \code{dc_predict}, plus the auxiliary
#'   \code{zero_effort} flag indicating cells with no GBIF records at
#'   all in the queried bins.
#'
#' @details
#' GBIF downloads can be large for well-surveyed countries. Reading a
#' multi-gigabyte file into memory for in-place year-bin filtering may
#' exceed available RAM. The adaptive binning strategy avoids this:
#' GBIF reports each download's compressed size in its metadata, and
#' bins whose download exceeds \code{max_size_bytes} are split into
#' halves and re-requested. The function iterates until all bins are
#' stable.
#'
#' For Belgium (the demonstration species' projection country), the
#' initial 8 bins typically remain unsplit. For larger countries with
#' denser sampling, the function may produce 12-20 final bins.
#'
#' All effort variables computed here use the identical aggregation
#' logic as \code{dc_extract_effort} and \code{dc_collapse_effort},
#' ensuring that projection predictions are comparable to training
#' data.
#'
#' @importFrom data.table data.table as.data.table setnames rbind
#' @importFrom terra rast mask vect values xyFromCell
#' @importFrom sf st_as_sf st_make_valid st_as_text st_union
#' @export
dc_project_country <- function(country_code,
                               year_bins,
                               raster_template,
                               presence_cells,
                               sinas_priors,
                               focal_species,
                               focal_exclusion_keys,
                               cfg,
                               keys_file              = NULL,
                               scratch_dir            = tempdir(),
                               max_size_bytes         = 8e8,
                               seven_zip_path         = "C:/Program Files/7-Zip/7z.exe",
                               baseline_prior         = 0.1,
                               vars_sel               = dc_default_vars(),
                               corroboration_radius_km = 50) {

  if (!requireNamespace("rgbif",   quietly = TRUE) ||
      !requireNamespace("geodata", quietly = TRUE))
    stop("Packages 'rgbif' and 'geodata' are required.")

  decimalLongitude <- decimalLatitude <- cell_id <- NULL
  national_presence <- confirmation_prior <- spatial_corroboration <- NULL
  zero_effort <- climate_class <- NULL

  # ── 1. Country grid cells ────────────────────────────────────
  country_sf <- sf::st_make_valid(sf::st_as_sf(
    geodata::gadm(country_code, level = 0, path = tempdir())))
  country_mask  <- terra::mask(raster_template, terra::vect(country_sf))
  country_cells <- which(!is.na(terra::values(country_mask)))

  country_xy <- data.table::as.data.table(
    terra::xyFromCell(raster_template, country_cells))
  data.table::setnames(country_xy,
                       c("decimalLongitude", "decimalLatitude"))
  country_xy[, cell_id := country_cells]
  message("Country grid cells: ", nrow(country_xy))

  country_wkt <- sf::st_as_text(sf::st_union(country_sf))

  # ── 2. Load or initialise key store ──────────────────────────
  keys_map <- if (!is.null(keys_file) && file.exists(keys_file))
    readRDS(keys_file) else list()

  # ── 3. Adaptive binning ──────────────────────────────────────
  refine_and_download <- function(bins, keys_map) {
    out <- list()
    for (bin in bins) {
      bin_label <- paste(bin, collapse = "-")
      key       <- keys_map[[bin_label]]

      if (is.null(key)) {
        message("Submitting NEW bin: ", bin_label)
        dl <- rgbif::occ_download(
          rgbif::pred("hasCoordinate", TRUE),
          rgbif::pred("geometry",      country_wkt),
          rgbif::pred_gte("year",      bin[1]),
          rgbif::pred_lte("year",      bin[2]),
          format = "SIMPLE_CSV")
        key <- if (is.character(dl)) dl else dl$key
        dc_wait_for_gbif(key)
        keys_map[[bin_label]] <- key
        if (!is.null(keys_file)) saveRDS(keys_map, keys_file)
      }

      meta <- tryCatch(rgbif::occ_download_meta(key),
                       error = function(e) NULL)

      if (!is.null(meta$size) && meta$size <= max_size_bytes) {
        out[[length(out) + 1]] <- bin
        next
      }
      message("Splitting bin: ", bin_label,
              " (", round(meta$size / 1e6, 1), " MB)")
      years <- seq(bin[1], bin[2])
      if (length(years) == 1L) {
        out[[length(out) + 1]] <- bin
      } else {
        mid <- floor(mean(years))
        out[[length(out) + 1]] <- c(bin[1], mid)
        out[[length(out) + 1]] <- c(mid + 1L, bin[2])
      }
    }
    list(bins = out, keys = keys_map)
  }

  projection_bins <- year_bins
  repeat {
    res          <- refine_and_download(projection_bins, keys_map)
    new_bins     <- res$bins
    keys_map     <- res$keys
    if (identical(new_bins, projection_bins)) break
    projection_bins <- new_bins
  }
  message("Final projection bins: ", length(projection_bins))

  # ── 4. Extract effort per bin ────────────────────────────────
  all_bin_aggs <- list()

  for (b in seq_along(projection_bins)) {
    bin       <- projection_bins[[b]]
    bin_label <- paste(bin, collapse = "-")
    key       <- keys_map[[bin_label]]
    message("\n--- ", country_code, " | ", bin_label, " ---")

    zip_path <- normalizePath(
      rgbif::occ_download_get(key, path = scratch_dir,
                              overwrite = TRUE),
      winslash = "/")

    agg_list <- dc_extract_effort(
      zip_path             = zip_path,
      year_bins            = list(bin),
      raster_template      = raster_template,
      focal_species        = focal_species,
      focal_exclusion_keys = focal_exclusion_keys,
      cfg                  = cfg,
      cell_filter          = country_cells,
      seven_zip_path       = seven_zip_path,
      tmpdir               = scratch_dir)

    file.remove(zip_path)
    all_bin_aggs[[bin_label]] <- agg_list[[1]]
  }

  # ── 5. Collapse and derive ───────────────────────────────────
  country_effort_grid <- dc_collapse_effort(all_bin_aggs)

  # ── 6. Merge onto all country cells ──────────────────────────
  country_surface <- merge(country_xy, country_effort_grid,
                           by = "cell_id", all.x = TRUE)
  country_surface[, zero_effort := is.na(get("n_records_all"))]

  effort_zero_cols <- intersect(
    c("log_total_records", "n_events",
      "relevance_ratio", "recorder_relevance_ratio",
      "recorder_dominance", "n_years_with_records",
      "recording_span", "recent_activity",
      "mean_coord_uncertainty"),
    names(country_surface))
  for (col in effort_zero_cols)
    country_surface[is.na(get(col)), (col) := 0]

  # ── 7. Add priors and climate ────────────────────────────────
  country_surface[, climate_class := as.integer(
    terra::values(raster_template)[cell_id])]
  country_surface[, national_presence := dc_extract_by_cell(
    sinas_priors$national_presence, cell_id)]
  country_surface[, confirmation_prior := dc_extract_by_cell(
    sinas_priors$confirmation_prior, cell_id)]
  country_surface[is.na(national_presence),  national_presence  := 0]
  country_surface[is.na(confirmation_prior), confirmation_prior := baseline_prior]

  # ── 8. Spatial corroboration ─────────────────────────────────
  message("Computing spatial corroboration for projection surface...")
  country_surface[, spatial_corroboration := dc_compute_corroboration(
    cell_id, presence_cells, raster_template,
    radius_km = corroboration_radius_km)]

  # ── 9. Sanity check ──────────────────────────────────────────
  missing <- setdiff(vars_sel, names(country_surface))
  if (length(missing))
    stop("Projection surface missing vars: ",
         paste(missing, collapse = ", "))

  message("\nProjection surface ready: ",
          nrow(country_surface), " cells")
  message("Zero-effort cells: ", sum(country_surface$zero_effort),
          " (", round(100 * mean(country_surface$zero_effort), 1), "%)")

  country_surface
}


#' Poll GBIF until a download is ready
#'
#' Convenience helper that blocks until a GBIF download key reaches
#' \code{SUCCEEDED} status, or errors out if the request fails or times
#' out. Used internally by \code{dc_project_country}.
#'
#' @param key GBIF download key (character).
#' @param poll_interval Seconds between status checks. Default 300.
#' @param max_wait Maximum total wait in seconds. Default 8 hours.
#'
#' @return \code{TRUE} (invisibly) on success; stops on failure or timeout.
#' @export
dc_wait_for_gbif <- function(key,
                             poll_interval = 300,
                             max_wait      = 3600 * 8) {

  if (!requireNamespace("rgbif", quietly = TRUE))
    stop("Package 'rgbif' is required.")

  start_time <- Sys.time()
  repeat {
    status <- rgbif::occ_download_meta(key)$status
    message("GBIF status: ", status, " | ", Sys.time())
    if (status == "SUCCEEDED") {
      message("Download complete: ", key); return(invisible(TRUE))
    }
    if (status %in% c("FAILED", "KILLED"))
      stop("GBIF download failed (status = ", status, "): ", key)
    if (as.numeric(difftime(Sys.time(), start_time, units = "secs")) > max_wait)
      stop("Timed out waiting for GBIF download: ", key)
    Sys.sleep(poll_interval)
  }
}
