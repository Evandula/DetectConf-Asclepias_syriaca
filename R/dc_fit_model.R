#' Fit a BART detection-confidence model
#'
#' Fits a Bayesian Additive Regression Trees (BART) model using
#' \code{dbarts::bart} with \code{keeptrees = TRUE}, then immediately
#' calls \code{mod$fit$storeState()} so the model can be serialised
#' with \code{saveRDS} without losing its C++ tree pointer.
#'
#' @param model_df Model frame produced by \code{dc_assemble_model_frame}.
#'   Must contain a \code{detected} column (0/1 response) and all
#'   columns in \code{vars_sel}.
#' @param vars_sel Character vector of predictor variable names.
#'   Defaults to \code{dc_default_vars()}.
#' @param out_file Optional file path to save the fitted model as an
#'   \code{.rds}. If \code{NULL}, the model is returned only in memory.
#' @param ... Additional arguments passed to \code{dbarts::bart}.
#'
#' @return A fitted BART model object with C++ trees stored ready for
#'   serialisation.
#'
#' @section BART serialisation:
#' \code{dbarts} BART objects rely on a C++ pointer to the fitted tree
#' state. When the object is serialised with \code{saveRDS}, that
#' pointer becomes invalid and \code{predict()} on the reloaded model
#' returns a uniform \code{0.5} for every input. This is silent — there
#' is no warning, no error, just useless predictions.
#'
#' \code{fit$storeState()} flattens the C++ state into R-side data
#' before serialisation. \code{dc_fit_model} always calls it, so
#' models saved with this function are safe to reload.
#'
#' Note: \code{exportState()} does not exist in current \code{dbarts}
#' versions and should not be used. Note also that
#' xz-compressed \code{.rds} files break the pointer irreversibly even
#' with \code{storeState()} — keep models gzip-compressed (the default).
#'
#' @examples
#' \dontrun{
#' mod <- dc_fit_model(model_df,
#'                     out_file = "inst/extdata/bart_model_final.rds")
#' # Later, in any R session:
#' mod_reloaded <- readRDS("inst/extdata/bart_model_final.rds")
#' preds <- dc_predict(mod_reloaded, pred_frame)
#' }
#'
#' @export
dc_fit_model <- function(model_df,
                         vars_sel = dc_default_vars(),
                         out_file = NULL,
                         ...) {

  if (!requireNamespace("dbarts", quietly = TRUE))
    stop("Package 'dbarts' is required. Install with install.packages('dbarts').")

  detected <- NULL

  if (!data.table::is.data.table(model_df))
    model_df <- data.table::as.data.table(model_df)

  missing_v <- setdiff(c(vars_sel, "detected"), names(model_df))
  if (length(missing_v))
    stop("model_df missing columns: ", paste(missing_v, collapse = ", "))

  x_train <- model_df[, vars_sel, with = FALSE]
  y_train <- model_df[["detected"]]

  message("Fitting BART model on ", nrow(model_df), " observations, ",
          ncol(x_train), " predictors...")

  mod <- dbarts::bart(x.train   = x_train,
                      y.train   = y_train,
                      keeptrees = TRUE,
                      ...)

  # CRITICAL: store C++ tree state before saveRDS
  mod$fit$storeState()
  message("BART trees stored for safe serialisation.")

  if (!is.null(out_file)) {
    dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
    saveRDS(mod, out_file)
    message("Model saved: ", out_file)
  }

  mod
}


#' Block-wise prediction from a BART model
#'
#' Generates predictions from a fitted BART model on a data frame of
#' new observations. Optionally splits the input into blocks for
#' memory-safe prediction on large surfaces (e.g. country projections).
#' Returns posterior mean predictions and optionally posterior
#' quantiles.
#'
#' @param object A BART model fitted by \code{dc_fit_model} (or reloaded
#'   via \code{readRDS}). The model must have been fitted with
#'   \code{keeptrees = TRUE} and have \code{storeState()} called before
#'   serialisation.
#' @param newdata Data frame or data.table containing the predictors
#'   used to fit the model. Extra columns are ignored.
#' @param quantiles Numeric vector of probabilities to compute from the
#'   posterior draws (e.g. \code{c(0.025, 0.975)} for a 95\% credible
#'   interval). Default \code{c()} returns only the posterior mean.
#' @param splitby Number of blocks to split prediction into. Default
#'   \code{1} (single block). For large projection surfaces, values of
#'   10-50 reduce peak memory.
#' @param quiet Suppress the progress bar shown during block-wise
#'   prediction. Default \code{FALSE}.
#'
#' @return If \code{quantiles} is empty, a numeric vector of posterior
#'   means with length \code{nrow(newdata)}. If \code{quantiles} are
#'   requested, a data frame with columns \code{pred} (posterior mean)
#'   and one column per requested quantile.
#'
#' @details Predictions are returned for rows of \code{newdata} that
#'   are complete cases on the model's predictors. Rows with any
#'   missing predictor values receive \code{NA} predictions.
#'
#'   Sanity check on reloaded models: if all predictions return
#'   exactly \code{0.5}, the model's C++ tree pointer was lost during
#'   serialisation. Refit the model and save it with
#'   \code{dc_fit_model}, which guarantees \code{storeState()} is
#'   called.
#'
#' @export
dc_predict <- function(object, newdata,
                       quantiles = c(),
                       splitby   = 1L,
                       quiet     = FALSE) {

  if (!requireNamespace("dbarts", quietly = TRUE))
    stop("Package 'dbarts' is required.")

  xnames    <- attr(object$fit$data@x, "term.labels")
  newdata_x <- as.data.frame(newdata)[, xnames, drop = FALSE]
  input_mat <- as.matrix(newdata_x)
  n_rows    <- nrow(input_mat)
  blankout  <- data.frame(matrix(NA_real_,
                                 ncol = 1 + length(quantiles),
                                 nrow = n_rows))

  which_valid <- which(stats::complete.cases(input_mat))
  input_mat   <- input_mat[which_valid, , drop = FALSE]

  if (nrow(input_mat) == 0L) {
    if (ncol(blankout) > 1)
      colnames(blankout) <- c("pred",
                              paste0("q", gsub("\\.", "", quantiles)))
    return(if (ncol(blankout) > 1) blankout else blankout[, 1])
  }

  summarise_block <- function(pred_mat) {
    if (length(quantiles) == 0L)
      colMeans(pred_mat)
    else
      cbind(data.frame(colMeans(pred_mat)),
            matrixStats::colQuantiles(pred_mat, probs = quantiles))
  }

  # Use the S3 generic: dbarts ships predict.bart but does not export
  # it from the dbarts namespace directly. Loading dbarts via
  # requireNamespace registers the S3 method, so the unqualified
  # predict() call dispatches correctly.
  loadNamespace("dbarts")

  if (splitby == 1L) {
    pred         <- predict(object, input_mat)
    pred_summary <- summarise_block(pred)
    blankout[which_valid, ] <- as.matrix(pred_summary)
  } else {
    block_sz <- ceiling(nrow(input_mat) / splitby)
    if (!quiet) pb <- utils::txtProgressBar(min = 0, max = splitby,
                                            style = 3)
    for (i in seq_len(splitby)) {
      idx_s <- (i - 1L) * block_sz + 1L
      idx_e <- min(i * block_sz, nrow(input_mat))
      if (idx_s > nrow(input_mat)) break
      blk    <- input_mat[idx_s:idx_e, , drop = FALSE]
      pred_b <- predict(object, blk)
      bs     <- summarise_block(pred_b)
      blankout[which_valid[idx_s:idx_e], ] <- as.matrix(bs)
      rm(pred_b, bs); gc()
      if (!quiet) utils::setTxtProgressBar(pb, i)
    }
    if (!quiet) close(pb)
  }

  if (ncol(blankout) > 1) {
    colnames(blankout) <- c("pred", paste0("q", gsub("\\.", "", quantiles)))
    blankout
  } else {
    blankout[, 1]
  }
}
