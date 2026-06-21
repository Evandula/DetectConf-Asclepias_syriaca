#' Assemble the BART model frame
#'
#' Wraps section 24 of the DetectConf pipeline. Combines presence cells
#' and pseudo-absence cells into a single data.table ready for BART
#' fitting, attaching climate, effort, knowledge-prior, and spatial
#' corroboration variables to each row. Imputes missing values, drops
#' zero-effort pseudo-absences, and validates that all required
#' predictors are present.
#'
#' @param presence_cells Integer vector of \code{cell_id} values for
#'   presence cells (typically from the non-native introduced range).
#' @param pseudo_dt The \code{pseudo_dt} component of the output of
#'   \code{dc_generate_pseudo_absences}.
#' @param effort_grid The collapsed effort table for presence-buffer
#'   downloads (output of \code{dc_collapse_effort} applied to section-20
#'   downloads).
#' @param pseudo_effort_grid The collapsed effort table for
#'   pseudo-absence downloads (output of \code{dc_collapse_effort}
#'   applied to section-23 PA-cluster downloads).
#' @param raster_template The classified climate envelope raster.
#' @param sinas_priors Named list returned by \code{dc_sinas_prior},
#'   containing \code{confirmation_prior} and \code{national_presence}
#'   rasters.
#' @param vars_sel Character vector of predictor variable names to
#'   validate before returning. Defaults to the full DetectConf
#'   variable set.
#' @param corroboration_radius_km Radius in kilometres for the spatial
#'   corroboration predictor. Default \code{50}.
#'
#' @return A \code{data.table} with one row per observation (presence
#'   or pseudo-absence), columns including \code{cell_id},
#'   \code{decimalLongitude}, \code{decimalLatitude},
#'   \code{climate_class}, \code{detected} (0/1 response), and all
#'   predictors in \code{vars_sel}.
#'
#' @details
#' \strong{Imputation}: Effort variables, taxonomic relevance ratios,
#' temporal coverage variables, spatial corroboration, and
#' \code{national_presence} default to \code{0} when missing.
#' \code{mean_coord_uncertainty} defaults to \code{0}.
#' \code{confirmation_prior} defaults to \code{0.1} (the baseline used
#' by \code{dc_sinas_prior}).
#'
#' \strong{Zero-effort pseudo-absences}: Pseudo-absences with
#' \code{log_total_records == 0} after merging are dropped. A
#' pseudo-absence with no effort contributes no information about
#' detection probability — only about candidate location — and would
#' otherwise dilute the signal.
#'
#' @importFrom data.table data.table as.data.table rbind setnames
#' @importFrom data.table copy merge.data.table
#' @importFrom terra xyFromCell values
#' @export
dc_assemble_model_frame <- function(presence_cells,
                                    pseudo_dt,
                                    effort_grid,
                                    pseudo_effort_grid,
                                    raster_template,
                                    sinas_priors,
                                    vars_sel = dc_default_vars(),
                                    corroboration_radius_km = 50) {

  cell_id <- decimalLongitude <- decimalLatitude <- climate_class <- NULL
  detected <- obs_type <- log_total_records <- NULL
  national_presence <- confirmation_prior <- spatial_corroboration <- NULL

  # ── Presence dataframe ────────────────────────────────────────
  pres_xy_dt <- data.table::as.data.table(
    terra::xyFromCell(raster_template, presence_cells))
  data.table::setnames(pres_xy_dt,
                       c("decimalLongitude", "decimalLatitude"))

  presence_df <- data.table::data.table(
    cell_id          = presence_cells,
    climate_class    = as.integer(
      dc_extract_by_cell(raster_template, presence_cells)),
    decimalLongitude = pres_xy_dt$decimalLongitude,
    decimalLatitude  = pres_xy_dt$decimalLatitude)

  presence_df <- merge(presence_df, effort_grid, by = "cell_id", all.x = TRUE)
  presence_df[, national_presence := dc_extract_by_cell(
    sinas_priors$national_presence, cell_id)]
  presence_df[, confirmation_prior := dc_extract_by_cell(
    sinas_priors$confirmation_prior, cell_id)]

  message("Computing spatial corroboration for presence cells...")
  presence_df[, spatial_corroboration := dc_compute_corroboration(
    cell_id, presence_cells, raster_template,
    radius_km = corroboration_radius_km)]

  presence_df[, detected := 1L]
  presence_df[, obs_type := "presence"]

  # ── Pseudo-absence dataframe ──────────────────────────────────
  pseudo_df <- data.table::copy(pseudo_dt)
  pseudo_df <- merge(pseudo_df, pseudo_effort_grid,
                     by = "cell_id", all.x = TRUE)
  pseudo_df[, detected := 0L]
  pseudo_df[, obs_type := "pseudo_absence"]
  pseudo_df[, national_presence := dc_extract_by_cell(
    sinas_priors$national_presence, cell_id)]
  pseudo_df[, confirmation_prior := dc_extract_by_cell(
    sinas_priors$confirmation_prior, cell_id)]

  message("Computing spatial corroboration for PA cells...")
  pseudo_df[, spatial_corroboration := dc_compute_corroboration(
    cell_id, presence_cells, raster_template,
    radius_km = corroboration_radius_km)]

  # Drop pseudo-absence-only auxiliary columns
  drop_cols <- intersect(c("continent", "cluster_id", "zero_effort",
                           "region", "cluster_key", "db_cluster"),
                         names(pseudo_df))
  if (length(drop_cols)) pseudo_df[, (drop_cols) := NULL]

  # ── Stack and validate ────────────────────────────────────────
  bart_cols <- c("cell_id", "decimalLongitude", "decimalLatitude",
                 "climate_class", "detected", "obs_type", vars_sel)

  missing_p <- setdiff(bart_cols, names(presence_df))
  missing_a <- setdiff(bart_cols, names(pseudo_df))
  if (length(missing_p))
    stop("Missing from presence_df: ", paste(missing_p, collapse = ", "))
  if (length(missing_a))
    stop("Missing from pseudo_df: ",   paste(missing_a, collapse = ", "))

  model_df <- rbind(presence_df[, bart_cols, with = FALSE],
                    pseudo_df[, bart_cols, with = FALSE])

  # ── Impute NAs ────────────────────────────────────────────────
  impute_zero <- intersect(c("log_total_records", "n_events",
                             "relevance_ratio", "recorder_relevance_ratio",
                             "recorder_dominance", "n_years_with_records",
                             "recording_span", "recent_activity",
                             "spatial_corroboration", "national_presence"),
                           names(model_df))
  for (col in impute_zero)
    model_df[is.na(get(col)), (col) := 0]

  if ("mean_coord_uncertainty" %in% names(model_df))
    model_df[is.na(get("mean_coord_uncertainty")),
             ("mean_coord_uncertainty") := 0]
  if ("confirmation_prior" %in% names(model_df))
    model_df[is.na(get("confirmation_prior")),
             ("confirmation_prior") := 0.1]

  # ── Drop zero-effort pseudo-absences ──────────────────────────
  n_before <- nrow(model_df)
  model_df <- model_df[!(detected == 0 & log_total_records == 0)]
  message("Zero-effort PAs removed: ", n_before - nrow(model_df))
  message("Final model frame: ", nrow(model_df), " rows")

  model_df
}


#' Build the projection-ready prediction frame
#'
#' Wraps section 25 of the DetectConf pipeline. Assembles a prediction
#' data.table for all cells covered by either presence-buffer effort or
#' pseudo-absence-buffer effort. Used internally by \code{dc_predict}
#' and to construct cross-validation comparisons.
#'
#' @param presence_cells Integer vector of presence cell IDs.
#' @param effort_grid Collapsed presence-buffer effort table.
#' @param pseudo_effort_grid Collapsed pseudo-absence-buffer effort
#'   table.
#' @param raster_template Classified climate envelope raster.
#' @param sinas_priors List of priors from \code{dc_sinas_prior}.
#' @param vars_sel Character vector of predictor variables to validate.
#'   Defaults to \code{dc_default_vars()}.
#' @param corroboration_radius_km Radius in kilometres for the spatial
#'   corroboration predictor. Default \code{50}.
#'
#' @return A \code{data.table} with one row per surveyed cell and all
#'   columns required for \code{dc_predict}.
#'
#' @importFrom data.table data.table as.data.table rbind setnames
#' @importFrom terra xyFromCell values
#' @export
dc_build_prediction_frame <- function(presence_cells,
                                      effort_grid,
                                      pseudo_effort_grid,
                                      raster_template,
                                      sinas_priors,
                                      vars_sel = dc_default_vars(),
                                      corroboration_radius_km = 50) {

  cell_id <- decimalLongitude <- decimalLatitude <- climate_class <- NULL
  national_presence <- confirmation_prior <- spatial_corroboration <- NULL

  pred_cells <- sort(unique(c(effort_grid$cell_id,
                              pseudo_effort_grid$cell_id,
                              presence_cells)))
  message("Total prediction cells: ", length(pred_cells))

  pred_xy <- data.table::as.data.table(
    terra::xyFromCell(raster_template, pred_cells))
  data.table::setnames(pred_xy, c("decimalLongitude", "decimalLatitude"))

  pred_frame <- data.table::data.table(
    cell_id          = pred_cells,
    decimalLongitude = pred_xy$decimalLongitude,
    decimalLatitude  = pred_xy$decimalLatitude,
    climate_class    = as.integer(terra::values(raster_template)[pred_cells]))

  effort_vars_join <- intersect(
    c("cell_id", "log_total_records", "n_events",
      "relevance_ratio", "recorder_relevance_ratio",
      "recorder_dominance", "n_years_with_records",
      "recording_span", "recent_activity",
      "mean_coord_uncertainty"),
    union(names(effort_grid), names(pseudo_effort_grid)))

  eff_pres   <- effort_grid[, intersect(effort_vars_join, names(effort_grid)),
                             with = FALSE]
  eff_pseudo <- pseudo_effort_grid[
    cell_id %in% setdiff(pseudo_effort_grid$cell_id, effort_grid$cell_id),
    intersect(effort_vars_join, names(pseudo_effort_grid)),
    with = FALSE]

  pred_frame <- merge(pred_frame,
                      rbind(eff_pres, eff_pseudo, fill = TRUE),
                      by = "cell_id", all.x = TRUE)

  for (col in setdiff(effort_vars_join, "cell_id"))
    pred_frame[is.na(get(col)), (col) := 0]

  pred_frame[, national_presence := dc_extract_by_cell(
    sinas_priors$national_presence, cell_id)]
  pred_frame[, confirmation_prior := dc_extract_by_cell(
    sinas_priors$confirmation_prior, cell_id)]
  pred_frame[is.na(national_presence),  national_presence  := 0]
  pred_frame[is.na(confirmation_prior), confirmation_prior := 0.1]

  message("Computing spatial corroboration for prediction frame...")
  pred_frame[, spatial_corroboration := dc_compute_corroboration(
    cell_id, presence_cells, raster_template,
    radius_km = corroboration_radius_km)]

  missing_pf <- setdiff(vars_sel, names(pred_frame))
  if (length(missing_pf))
    stop("pred_frame missing: ", paste(missing_pf, collapse = ", "))

  message("Prediction frame ready: ", nrow(pred_frame), " cells")
  pred_frame
}


#' Default DetectConf predictor variable set
#'
#' Returns the canonical set of 13 predictor variables used by
#' DetectConf models. Provided as a function so the canonical list lives
#' in exactly one place.
#'
#' @return Character vector of predictor variable names.
#'
#' @export
dc_default_vars <- function() {
  c("climate_class",
    "log_total_records",
    "n_events",
    "relevance_ratio",
    "recorder_relevance_ratio",
    "recorder_dominance",
    "n_years_with_records",
    "recording_span",
    "recent_activity",
    "mean_coord_uncertainty",
    "spatial_corroboration",
    "national_presence",
    "confirmation_prior")
}
