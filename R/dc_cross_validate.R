#' Spatial buffer cross-validation for DetectConf models
#'
#' Implements leave-one-buffer-out spatial
#' cross-validation: for each fold, a
#' randomly selected presence cell becomes the focal point, all
#' observations within \code{buffer_km} are held out as the test set,
#' and the remaining observations train the model. Quality is measured
#' with AUC, Boyce continuous index, TSS, and Type I / Type II error
#' rates.
#'
#' The buffer radius is computed as a fraction of the convex-hull
#' radius of all presence points, clamped between sensible minimum and
#' maximum values, so the procedure adapts to species ranges of any
#' size.
#'
#' @param model_df Model frame produced by
#'   \code{dc_assemble_model_frame}.
#' @param raster_template Classified climate envelope raster.
#' @param presence_coords Either a data.frame/data.table with columns
#'   \code{decimalLongitude} and \code{decimalLatitude}, or a
#'   \code{SpatVector} of presence points. Used to size the buffer
#'   automatically.
#' @param vars_sel Character vector of predictor variables.
#' @param n_folds Maximum number of folds to evaluate. Default \code{20}.
#'   Fewer may be returned if not enough valid focal cells exist.
#' @param mcp_buffer_frac Fraction of the convex-hull radius to use as
#'   the buffer. Default \code{0.10}.
#' @param buffer_min_km Minimum buffer in km. Default \code{100}.
#' @param buffer_max_km Maximum buffer in km. Default \code{2000}.
#' @param min_test_pres Minimum presences required in the test buffer
#'   for a fold to be considered valid. Default \code{5}.
#' @param min_test_abs Minimum pseudo-absences required in the test
#'   buffer. Default \code{5}.
#' @param test_prev_range Tuple of \code{c(min, max)} test-buffer
#'   prevalences considered valid. Default \code{c(0.2, 0.8)}.
#' @param max_prev_diff Maximum allowed difference between test and
#'   train prevalences. Default \code{0.25}.
#' @param save_final_model_to Optional file path. If supplied, the
#'   model fitted in the last fold is saved (with \code{storeState()}
#'   already called) for downstream projection use.
#' @param seed Random seed for focal-cell sampling. Default \code{45}.
#'
#' @return A \code{data.table} with one row per fold containing
#'   \code{fold}, \code{focal_cell}, \code{buffer_km}, \code{n_train},
#'   \code{n_test}, \code{auc}, \code{boyce}, \code{tss},
#'   \code{type1_error}, \code{type2_error}.
#'
#' @importFrom data.table data.table as.data.table rbindlist fifelse
#' @importFrom terra vect convHull crs expanse xyFromCell
#' @importFrom fields rdist.earth
#' @export
dc_cross_validate <- function(model_df,
                              raster_template,
                              presence_coords,
                              vars_sel        = dc_default_vars(),
                              n_folds         = 20L,
                              mcp_buffer_frac = 0.10,
                              buffer_min_km   = 100,
                              buffer_max_km   = 2000,
                              min_test_pres   = 5L,
                              min_test_abs    = 5L,
                              test_prev_range = c(0.2, 0.8),
                              max_prev_diff   = 0.25,
                              save_final_model_to = NULL,
                              seed                = 45L) {

  if (!requireNamespace("dbarts", quietly = TRUE) ||
      !requireNamespace("pROC",   quietly = TRUE) ||
      !requireNamespace("modEvA", quietly = TRUE))
    stop("Packages 'dbarts', 'pROC', and 'modEvA' are required.")

  cell_id <- detected <- decimalLongitude <- decimalLatitude <- NULL
  dist_to_focal <- NULL

  if (!data.table::is.data.table(model_df))
    model_df <- data.table::as.data.table(model_df)

  # ── Size buffer from MCP radius ──────────────────────────────
  if (inherits(presence_coords, "SpatVector")) {
    pts_vc_cv <- presence_coords
  } else {
    pts_vc_cv <- terra::vect(presence_coords,
                             geom = c("decimalLongitude", "decimalLatitude"),
                             crs  = "EPSG:4326")
  }
  mcp_global <- terra::convHull(pts_vc_cv)
  terra::crs(mcp_global) <- "EPSG:4326"

  mcp_area_km2  <- terra::expanse(mcp_global, unit = "km")
  mcp_radius_km <- sqrt(mcp_area_km2 / pi)
  buffer_km     <- round(mcp_radius_km * mcp_buffer_frac)
  buffer_km     <- min(max(buffer_km, buffer_min_km), buffer_max_km)
  message("MCP radius (km): ", round(mcp_radius_km),
          " | Buffer (km): ", buffer_km)

  # ── Pre-screen focal cells ───────────────────────────────────
  set.seed(seed)
  candidate_cells_cv <- sample(model_df[detected == 1, cell_id],
                               min(100L, model_df[detected == 1, .N]))
  valid_focals <- integer()
  message("Pre-screening focal cells...")

  for (focal_cell in candidate_cells_cv) {
    focal_xy <- terra::xyFromCell(raster_template, focal_cell)
    all_xy   <- model_df[, .(decimalLongitude, decimalLatitude)]
    dists_km <- fields::rdist.earth(
      matrix(c(focal_xy[1], focal_xy[2]), nrow = 1L),
      as.matrix(all_xy), miles = FALSE)[1, ]
    model_df[, dist_to_focal := dists_km]
    test_df  <- model_df[dist_to_focal <= buffer_km]
    train_df <- model_df[dist_to_focal >  buffer_km]
    tp <- sum(test_df$detected == 1); ta <- sum(test_df$detected == 0)
    if ((tp + ta) == 0) next
    test_prev  <- tp / (tp + ta)
    train_prev <- sum(train_df$detected == 1) / nrow(train_df)
    if (tp >= min_test_pres && ta >= min_test_abs &&
        test_prev >= test_prev_range[1] &&
        test_prev <= test_prev_range[2] &&
        abs(test_prev - train_prev) <= max_prev_diff &&
        nrow(train_df) >= 50) {
      valid_focals <- c(valid_focals, focal_cell)
      if (length(valid_focals) >= n_folds) break
    }
  }
  model_df[, dist_to_focal := NULL]

  if (length(valid_focals) == 0L)
    stop("No valid focal cells. Reduce buffer or relax thresholds.")
  if (length(valid_focals) < n_folds)
    message("Note: only ", length(valid_focals),
            " valid folds found (requested ", n_folds, ").")

  # ── Per-fold evaluation ──────────────────────────────────────
  buffer_eval_results <- list()

  for (i in seq_along(valid_focals)) {
    focal_cell <- valid_focals[i]
    focal_xy   <- terra::xyFromCell(raster_template, focal_cell)
    message("\nFold ", i, "/", length(valid_focals))

    dists_km <- fields::rdist.earth(
      matrix(c(focal_xy[1], focal_xy[2]), nrow = 1L),
      as.matrix(model_df[, .(decimalLongitude, decimalLatitude)]),
      miles = FALSE)[1, ]
    model_df[, dist_to_focal := dists_km]
    test_df  <- model_df[dist_to_focal <= buffer_km]
    train_df <- model_df[dist_to_focal >  buffer_km]

    message("Train: ", nrow(train_df),
            " | Test: ", nrow(test_df),
            " | pres: ", sum(test_df$detected == 1),
            " | abs: ",  sum(test_df$detected == 0))

    mod_buf <- dbarts::bart(x.train   = train_df[, vars_sel, with = FALSE],
                            y.train   = train_df[["detected"]],
                            keeptrees = TRUE)

    test_preds <- dc_predict(mod_buf, test_df[, vars_sel, with = FALSE])

    roc_i   <- pROC::roc(test_df$detected, test_preds, quiet = TRUE)
    auc_i   <- as.numeric(pROC::auc(roc_i))
    boyce_i <- suppressWarnings(tryCatch(
      modEvA::Boyce(obs = test_df$detected, pred = test_preds)[["Boyce"]],
      error = function(e) NA_real_))
    coords_i <- pROC::coords(roc_i, "best",
                             best.method = "closest.topleft",
                             ret = c("threshold", "sensitivity",
                                     "specificity"))
    if (is.data.frame(coords_i) && nrow(coords_i) > 1)
      coords_i <- coords_i[1, ]
    tss_i         <- coords_i$sensitivity + coords_i$specificity - 1
    pred_binary_i <- as.integer(test_preds >= coords_i$threshold)
    tp_i <- sum(pred_binary_i == 1 & test_df$detected == 1)
    tn_i <- sum(pred_binary_i == 0 & test_df$detected == 0)
    fp_i <- sum(pred_binary_i == 1 & test_df$detected == 0)
    fn_i <- sum(pred_binary_i == 0 & test_df$detected == 1)
    type1_i <- data.table::fifelse(fp_i + tn_i > 0,
                                   fp_i / (fp_i + tn_i), NA_real_)
    type2_i <- data.table::fifelse(fn_i + tp_i > 0,
                                   fn_i / (fn_i + tp_i), NA_real_)

    message("AUC=",     round(auc_i,    3),
            " | Boyce=", round(boyce_i,  3),
            " | TSS=",   round(tss_i,    3),
            " | Type I=",  round(type1_i, 3),
            " | Type II=", round(type2_i, 3))

    buffer_eval_results[[i]] <- data.table::data.table(
      fold = i, focal_cell = focal_cell, buffer_km = buffer_km,
      n_train = nrow(train_df), n_test = nrow(test_df),
      auc = auc_i, boyce = boyce_i, tss = tss_i,
      type1_error = type1_i, type2_error = type2_i)

    # Save final model
    if (i == length(valid_focals) && !is.null(save_final_model_to)) {
      mod_buf$fit$storeState()
      dir.create(dirname(save_final_model_to),
                 recursive = TRUE, showWarnings = FALSE)
      saveRDS(mod_buf, save_final_model_to)
      message("Final model saved with trees stored: ",
              save_final_model_to)
    }

    rm(mod_buf, test_preds, roc_i); gc()
  }
  model_df[, dist_to_focal := NULL]

  out <- data.table::rbindlist(buffer_eval_results, fill = TRUE)
  message("\n--- CV SUMMARY ---")
  message("Folds      : ", nrow(out))
  message("Mean AUC   : ", round(mean(out$auc,   na.rm = TRUE), 3))
  message("Mean Boyce : ", round(mean(out$boyce, na.rm = TRUE), 3))
  message("Mean TSS   : ", round(mean(out$tss,   na.rm = TRUE), 3))
  out
}
