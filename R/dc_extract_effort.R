#' Extract per-cell effort statistics from a GBIF download zip
#'
#' Core function of the DetectConf pipeline. Streams a GBIF
#' SIMPLE_CSV zip through 7-Zip without unzipping to disk, applies the
#' standard cleaning rules (focal-species removal, focal-observer
#' exclusion, coordinate validation, grid-cell assignment), tags each
#' record as taxonomically relevant or not, and aggregates per-cell
#' effort statistics within one or more year bins.
#'
#' Used identically in three pipeline contexts: (1) accumulating effort
#' for cells in the presence buffer;
#' (2) extracting effort for pseudo-absence cells per cluster
#' (section 23); and (3) computing effort for projection-region cells
#' for the country/location of interest. The aggregation logic is fixed so that
#' all three contexts produce mutually compatible variables.
#'
#' @param zip_path Path to a GBIF SIMPLE_CSV zip download.
#' @param year_bins A list of length-2 integer vectors of the form
#'   \code{list(c(1800, 2005), c(2006, 2009), ...)} specifying the
#'   temporal bins to aggregate. If a single bin is supplied, the
#'   function returns the aggregate for that bin only.
#' @param raster_template A \code{terra::SpatRaster} (the classified
#'   climate envelope) used to assign \code{cell_id} via
#'   \code{terra::cellFromXY}. All three call sites must share the same
#'   raster so that \code{cell_id} is consistent across the pipeline.
#' @param focal_species Scientific name of the target invasive species
#'   (e.g. \code{"Asclepias syriaca"}). Records matching this name are
#'   removed to prevent circularity: the focal species cannot
#'   simultaneously be the signal of interest and a component of the
#'   effort proxy.
#' @param focal_exclusion_keys Character vector of
#'   \code{"cell_id|recordedBy"} keys identifying observer-cell pairs
#'   that recorded the focal species. Records from these observers in
#'   these cells are excluded to remove observer-level circularity.
#' @param cfg A taxonomic relevance specification (one element of
#'   \code{DETECTCONF_EXPERTISE}, see \code{?dc_expertise_config}). May
#'   contain \code{phylum}, \code{orders}, and/or \code{exclude_phyla}
#'   character vectors. Records meeting the criteria are tagged
#'   \code{is_taxon_relevant = TRUE}.
#' @param cell_filter Optional integer vector of \code{cell_id} values
#'   to restrict aggregation to. If \code{NULL} (default), all cells
#'   encountered in the download are returned. Used by pseudo-absence
#'   and projection extraction to avoid wasted computation on cells
#'   outside the area of interest.
#' @param seven_zip_path Path to the 7-Zip executable. Default is the
#'   Windows install location; on macOS or Linux, supply the path to
#'   \code{7z} or \code{p7zip}.
#' @param tmpdir Temporary directory for \code{data.table::fread}.
#'   Default is \code{tempdir()}.
#'
#' @return A named list of \code{data.table} objects, one per bin. Each
#'   data.table has columns:
#'   \describe{
#'     \item{\code{cell_id}}{integer, the grid cell index in
#'       \code{raster_template}}
#'     \item{\code{n_records_all}}{total records in the cell}
#'     \item{\code{n_records_relevant}}{records meeting taxonomic
#'       relevance criteria}
#'     \item{\code{n_recorders}, \code{n_recorders_relevant}}{unique
#'       recorder counts}
#'     \item{\code{n_events}}{unique event proxies (year + recorder +
#'       coordinates rounded to 3 dp)}
#'     \item{\code{n_years_with_records}, \code{year_min},
#'       \code{year_max}}{temporal coverage}
#'     \item{\code{sum_uncertainty}, \code{n_uncertainty}}{components for
#'       computing \code{mean_coord_uncertainty}}
#'     \item{\code{recorder_max_bin}}{maximum records by a single
#'       recorder (for \code{recorder_dominance})}
#'   }
#'   Names of the returned list are the bin labels
#'   (e.g. \code{"2006-2009"}).
#'
#' @details
#' \strong{Taxonomic relevance} describes \emph{what} was recorded in
#' the cell, not recorder competence. A cell with many beetle records
#' will have high \code{relevance_ratio} if the focal species is a
#' beetle, regardless of whether the recorders were entomologists.
#'
#' The function trusts that records in \code{focal_exclusion_keys} have
#' already been computed from the cleaned presence layer (see
#' \code{?dc_focal_exclusion_keys}). It does not re-derive them.
#'
#' @examples
#' \dontrun{
#' cfg <- dc_expertise_config("terrestrial_vascular_plants")
#' effort_by_bin <- dc_extract_effort(
#'   zip_path             = "0001234-260101000000000.zip",
#'   year_bins            = list(c(1800, 2005), c(2006, 2020),
#'                               c(2021, 2026)),
#'   raster_template      = mean_r,
#'   focal_species        = "Asclepias syriaca",
#'   focal_exclusion_keys = exclusion_keys,
#'   cfg                  = cfg
#' )
#' }
#'
#' @importFrom data.table fread uniqueN
#' @importFrom terra cellFromXY
#' @export
dc_extract_effort <- function(zip_path,
                              year_bins,
                              raster_template,
                              focal_species,
                              focal_exclusion_keys,
                              cfg,
                              cell_filter    = NULL,
                              seven_zip_path = "C:/Program Files/7-Zip/7z.exe",
                              tmpdir         = tempdir()) {

  # Declare data.table column variables to satisfy R CMD check
  decimalLongitude <- decimalLatitude <- countryCode <- NULL
  scientificName <- phylum <- order <- recordedBy <- year <- NULL
  cell_id <- exclusion_key <- coordinateUncertaintyInMeters <- NULL
  coord_uncert <- event_id_proxy <- is_taxon_relevant <- NULL
  basisOfRecord <- NULL

  if (!file.exists(zip_path))
    stop("Zip file not found: ", zip_path)

  message("Extracting from: ", basename(zip_path),
          " (", round(file.size(zip_path) / 1e9, 2), " GB)")

  # ── Stream zip through 7-Zip ─────────────────────────────────
  dt <- data.table::fread(
    cmd = sprintf('"%s" e -so "%s"', seven_zip_path, zip_path),
    sep = "\t", quote = "", fill = TRUE,
    select = c("decimalLongitude", "decimalLatitude", "year",
               "countryCode", "scientificName", "phylum", "order",
               "recordedBy", "basisOfRecord",
               "coordinateUncertaintyInMeters"),
    showProgress = TRUE, tmpdir = tmpdir)

  message("Rows read: ", formatC(nrow(dt), format = "d", big.mark = ","))

  # ── Cleaning step 1: coordinate validity ─────────────────────
  dt <- dt[is.finite(decimalLongitude) & is.finite(decimalLatitude)]

  # ── Cleaning step 2: remove focal species records ────────────
  n_before <- nrow(dt)
  dt <- dt[!grepl(focal_species, scientificName, ignore.case = TRUE)]
  message("Focal species records removed: ", n_before - nrow(dt))

  # ── Cleaning step 3: tag taxonomic relevance ─────────────────
  dt[, is_taxon_relevant := TRUE]
  if (!is.null(cfg$phylum))
    dt[!phylum %in% cfg$phylum,        is_taxon_relevant := FALSE]
  if (!is.null(cfg$orders))
    dt[!order  %in% cfg$orders,        is_taxon_relevant := FALSE]
  if (!is.null(cfg$exclude_phyla))
    dt[ phylum %in% cfg$exclude_phyla, is_taxon_relevant := FALSE]

  # ── Cleaning step 4: assign cell_id ──────────────────────────
  dt[, cell_id := terra::cellFromXY(
    raster_template, cbind(decimalLongitude, decimalLatitude))]
  dt <- dt[!is.na(cell_id)]

  # ── Cleaning step 5: optional cell filter ────────────────────
  if (!is.null(cell_filter))
    dt <- dt[cell_id %in% cell_filter]

  # ── Cleaning step 6: remove focal observers from focal cells ─
  dt[, exclusion_key := paste(cell_id, recordedBy, sep = "|")]
  n_before <- nrow(dt)
  dt <- dt[!exclusion_key %in% focal_exclusion_keys]
  dt[, exclusion_key := NULL]
  message("Focal observer records removed: ", n_before - nrow(dt))

  if (nrow(dt) == 0) {
    message("No records remain after cleaning.")
    return(setNames(
      replicate(length(year_bins), data.table::data.table(), simplify = FALSE),
      vapply(year_bins, function(b) paste(b, collapse = "-"), character(1))))
  }

  # ── Derived per-record fields ─────────────────────────────────
  dt[, coord_uncert := as.numeric(coordinateUncertaintyInMeters)]
  dt[!is.finite(coord_uncert), coord_uncert := NA_real_]

  dt[, event_id_proxy := paste(year, recordedBy,
                               round(decimalLongitude, 3),
                               round(decimalLatitude,  3))]

  # ── Aggregate per bin ─────────────────────────────────────────
  out <- vector("list", length(year_bins))
  names(out) <- vapply(year_bins,
                       function(b) paste(b, collapse = "-"), character(1))

  for (b in seq_along(year_bins)) {
    bin       <- year_bins[[b]]
    bin_label <- paste(bin, collapse = "-")
    dt_b      <- dt[year >= bin[1] & year <= bin[2]]

    if (nrow(dt_b) == 0) {
      message("  Bin ", bin_label, ": no records")
      out[[b]] <- data.table::data.table()
      next
    }

    message("  Bin ", bin_label, ": ",
            formatC(nrow(dt_b), format = "d", big.mark = ","),
            " records | cells: ",
            formatC(data.table::uniqueN(dt_b$cell_id),
                    format = "d", big.mark = ","))

    out[[b]] <- dt_b[, .(
      n_records_all        = .N,
      n_records_relevant   = sum(is_taxon_relevant),
      n_recorders          = data.table::uniqueN(recordedBy),
      n_recorders_relevant = data.table::uniqueN(
                               recordedBy[is_taxon_relevant == TRUE]),
      n_events             = data.table::uniqueN(event_id_proxy),
      n_years_with_records = data.table::uniqueN(year),
      year_min             = min(year, na.rm = TRUE),
      year_max             = max(year, na.rm = TRUE),
      sum_uncertainty      = sum(coord_uncert, na.rm = TRUE),
      n_uncertainty        = sum(!is.na(coord_uncert)),
      recorder_max_bin     = {
        rc <- table(recordedBy)
        if (length(rc) > 0) max(rc) else 0L }
    ), by = cell_id]
  }

  rm(dt); gc()
  out
}


#' Collapse multiple per-bin effort tables into a single grid-level table
#'
#' Aggregates a list of per-bin effort tables (e.g. the output of
#' \code{dc_extract_effort} applied across multiple zips and bins) into
#' one row per \code{cell_id}, summing record counts and taking the
#' maximum for recorder counts. Then derives the modelling variables
#' (\code{log_total_records}, \code{relevance_ratio},
#' \code{recorder_relevance_ratio}, \code{recorder_dominance},
#' \code{mean_coord_uncertainty}, \code{recording_span},
#' \code{recent_activity}).
#'
#' @param effort_list A list of \code{data.table} objects with the
#'   schema returned by \code{dc_extract_effort}. Empty data.tables
#'   are silently ignored.
#' @param current_year Integer year used for the \code{recent_activity}
#'   flag. Default \code{as.integer(format(Sys.Date(), "\%Y"))}.
#'
#' @return A single \code{data.table} with one row per \code{cell_id}
#'   and the full set of derived effort variables ready to merge into
#'   a presence, pseudo-absence, or projection data frame.
#'
#' @details Recorder counts are aggregated with \code{max} rather than
#'   \code{sum} because the same recorder may appear in multiple bins;
#'   summing would double-count them. This is a conservative
#'   approximation — the true unique count across bins could be slightly
#'   higher than the per-bin maximum but slightly lower than the sum.
#'
#' @importFrom data.table rbindlist fifelse
#' @export
dc_collapse_effort <- function(effort_list,
                               current_year = as.integer(
                                 format(Sys.Date(), "%Y"))) {

  cell_id <- n_records_all <- n_records_relevant <- NULL
  n_recorders <- n_recorders_relevant <- n_events <- NULL
  n_years_with_records <- year_min <- year_max <- NULL
  sum_uncertainty <- n_uncertainty <- recorder_max_bin <- NULL

  # Drop empty entries
  keep <- vapply(effort_list,
                 function(x) is.data.frame(x) && nrow(x) > 0, logical(1))
  if (!any(keep))
    return(data.table::data.table())

  acc <- data.table::rbindlist(effort_list[keep], fill = TRUE)

  out <- acc[, .(
    n_records_all        = sum(n_records_all,        na.rm = TRUE),
    n_records_relevant   = sum(n_records_relevant,   na.rm = TRUE),
    n_recorders          = max(n_recorders,          na.rm = TRUE),
    n_recorders_relevant = max(n_recorders_relevant, na.rm = TRUE),
    n_events             = sum(n_events,             na.rm = TRUE),
    n_years_with_records = sum(n_years_with_records, na.rm = TRUE),
    year_min             = min(year_min,             na.rm = TRUE),
    year_max             = max(year_max,             na.rm = TRUE),
    sum_uncertainty      = sum(sum_uncertainty,      na.rm = TRUE),
    n_uncertainty        = sum(n_uncertainty,        na.rm = TRUE),
    recorder_max         = max(recorder_max_bin,     na.rm = TRUE)
  ), by = cell_id]

  # Derived variables (section 22 logic)
  recorder_max <- relevance_ratio <- recorder_relevance_ratio <- NULL
  recorder_dominance <- mean_coord_uncertainty <- NULL
  recording_span <- recent_activity <- log_total_records <- NULL

  out[, `:=`(
    log_total_records        = log1p(n_records_all),
    relevance_ratio          = data.table::fifelse(
      n_records_all > 0, n_records_relevant / n_records_all, 0),
    recorder_relevance_ratio = data.table::fifelse(
      n_recorders > 0, n_recorders_relevant / n_recorders, 0),
    recorder_dominance       = data.table::fifelse(
      n_records_all > 0, recorder_max / n_records_all, 0),
    mean_coord_uncertainty   = data.table::fifelse(
      n_uncertainty > 0, sum_uncertainty / n_uncertainty, NA_real_),
    recording_span           = data.table::fifelse(
      n_years_with_records > 1L, year_max - year_min, 0L),
    recent_activity          = data.table::fifelse(
      year_max >= (current_year - 5L), 1L, 0L)
  )]

  out[, c("sum_uncertainty", "n_uncertainty", "recorder_max") := NULL]
  out
}


#' Derive the focal-observer exclusion keys
#'
#' Builds the character vector of \code{"cell_id|recordedBy"} keys used
#' by \code{dc_extract_effort} to remove records contributed by the same
#' observer who reported the focal species in the same cell. This
#' breaks the circularity that would otherwise inflate effort estimates
#' in cells where the focal species was found.
#'
#' @param presence_records A data.table with columns \code{cell_id} and
#'   \code{recordedBy}. Typically the cleaned presence records produced
#'   by \code{dc_climate_envelope} downstream processing.
#'
#' @return Character vector of unique \code{"cell_id|recordedBy"} keys.
#'
#' @export
dc_focal_exclusion_keys <- function(presence_records) {

  cell_id <- recordedBy <- exclusion_key <- NULL

  if (!data.table::is.data.table(presence_records))
    presence_records <- data.table::as.data.table(presence_records)

  focal_obs <- unique(presence_records[
    !is.na(recordedBy) & nchar(trimws(recordedBy)) > 0,
    .(cell_id, recordedBy)])
  focal_obs[, exclusion_key := paste(cell_id, recordedBy, sep = "|")]
  focal_obs$exclusion_key
}
