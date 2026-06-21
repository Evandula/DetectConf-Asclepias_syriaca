################################################################
## DetectConf demo: Asclepias syriaca in Belgium
##
## Demonstrates the full DetectConf workflow using the cached
## artefacts shipped in inst/extdata/. This runs in minutes
## (versus the multi-day full pipeline) because the heavy
## GBIF downloads and effort aggregation are already done.
##
## To run the full pipeline from scratch instead, see
## inst/scripts/full_pipeline.R.
################################################################

library(DetectConf)
library(terra)
library(data.table)

## ── 1. Load cached artefacts ────────────────────────────────────
ext <- function(x) system.file("extdata", x, package = "DetectConf")

mean_r            <- rast(ext("climate_envelope.tif"))
conf_prior_rast   <- rast(ext("sinas_confirmation_prior.tif"))
nat_presence_rast <- rast(ext("sinas_national_presence.tif"))
sinas_priors      <- list(confirmation_prior = conf_prior_rast,
                          national_presence  = nat_presence_rast)

effort_grid        <- readRDS(ext("effort_grid_final.rds"))
pseudo_effort_grid <- readRDS(ext("pseudo_effort_grid_final.rds"))
belgium_surface    <- readRDS(ext("belgium_projection_surface.rds"))
bart_model         <- readRDS(ext("bart_model_final.rds"))
cv_results         <- readRDS(ext("cv_results.rds"))

## ── 2. Cross-validation summary ─────────────────────────────────
message("\nCross-validation results:")
print(cv_results[, .(fold, auc, boyce, tss,
                     type1_error, type2_error)])
message("Mean AUC   : ", round(mean(cv_results$auc,   na.rm = TRUE), 3))
message("Mean Boyce : ", round(mean(cv_results$boyce, na.rm = TRUE), 3))
message("Mean TSS   : ", round(mean(cv_results$tss,   na.rm = TRUE), 3))

## ── 3. Predict over Belgium ─────────────────────────────────────
message("\nGenerating Belgium detection-confidence surface...")
belgium_preds <- dc_predict(bart_model, belgium_surface,
                            quantiles = c(0.025, 0.975),
                            splitby   = 5)

belgium_surface[, detect_confidence := belgium_preds$pred]
belgium_surface[, dc_lower          := belgium_preds$q0025]
belgium_surface[, dc_upper          := belgium_preds$q0975]

## ── 4. Quick map ────────────────────────────────────────────────
# A bare-bones map of detection confidence over Belgium. For
# publication-quality maps see vignettes/belgium_demo.Rmd.
plot_r <- rast(mean_r)
values(plot_r) <- NA
values(plot_r)[belgium_surface$cell_id] <- belgium_surface$detect_confidence

plot(plot_r,
     main = "DetectConf: Asclepias syriaca in Belgium",
     col  = hcl.colors(50, "viridis"))

## ── 5. Tripartition (extension beyond core claims) ──────────────
# Three operational outputs anchored on detection confidence:
#   credible detection   : high confidence + record exists
#   monitoring priority  : high confidence + no record + suitable
#   surveillance gap     : low confidence + suitable

belgium_surface[, status := fcase(
  detect_confidence >= 0.7 &  zero_effort == FALSE, "credible_detection",
  detect_confidence >= 0.7 &  zero_effort == TRUE,  "monitoring_priority",
  detect_confidence <  0.4,                          "surveillance_gap",
  default = "intermediate")]

print(belgium_surface[, .N, by = status])
