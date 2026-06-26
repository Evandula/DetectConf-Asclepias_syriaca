#' Generate pseudo-absences within invasion-front envelopes
#'
#' Performs global DBSCAN
#' clustering on presence points to identify natural invasion fronts,
#' draws a buffered convex hull around each cluster, and samples
#' pseudo-absences from candidate cells inside those hulls.
#'
#' Pseudo-absences are defined as locations within the spatial context
#' where the species could plausibly occur \emph{and} where the
#' observation system was demonstrably active. They are explicitly not
#' inferred true absences. This restricted sampling avoids the common
#' pitfall of treating every unsurveyed cell as evidence of absence.
#'
#' @param presence_xy A data.table or data.frame with columns
#'   \code{decimalLongitude} and \code{decimalLatitude} for cleaned
#'   presence points (typically the introduced-range subset).
#' @param raster_template The classified climate envelope raster (output
#'   of \code{dc_climate_envelope}). Cells with value 3 (native range)
#'   are excluded from candidate sampling.
#' @param pa_ratio Number of pseudo-absences per presence cell within
#'   each cluster. Default \code{3}.
#' @param target_cluster_size Target k-means batch size for grouping
#'   pseudo-absences into GBIF download requests. Default \code{20}.
#' @param dbscan_eps DBSCAN epsilon in degrees. If \code{NULL}, auto-tuned
#'   from the kNN-distance curve. Default \code{3} degrees.
#' @param dbscan_min_pts DBSCAN minimum cluster size. Default \code{5}.
#' @param hull_buffer_deg Buffer in degrees around each cluster's convex
#'   hull. Default \code{2}.
#' @param exclusion_radius_km Pseudo-absences are forbidden within this
#'   distance of any cluster-member presence point. Default \code{50} km.
#' @param seed Integer seed for reproducible sampling. Default \code{42}.
#'
#' @return A named list with elements:
#'   \describe{
#'     \item{\code{pseudo_dt}}{A \code{data.table} of sampled
#'       pseudo-absence cells with columns \code{cell_id},
#'       \code{db_cluster}, \code{cluster_key},
#'       \code{decimalLongitude}, \code{decimalLatitude},
#'       \code{climate_class}, \code{continent}, \code{cluster_id} (the
#'       k-means batch label used for GBIF downloads).}
#'     \item{\code{cluster_hulls}}{A named list of \code{sf} polygon
#'       objects, one per DBSCAN cluster.}
#'     \item{\code{noise_presences}}{A \code{data.table} of presence
#'       points classified as DBSCAN noise. Reserved for external
#'       validation and not used in the model frame.}
#'     \item{\code{cluster_wkt_list}}{A named list of simplified WKT
#'       polygon strings (one per k-means batch), suitable for passing
#'       to \code{rgbif::pred("geometry", ...)} in GBIF download
#'       requests for pseudo-absence effort extraction.}
#'   }
#'
#' @details
#' The function uses fixed region polygons (continent-scale bounding
#' boxes) for tagging only — clustering itself is global, so an
#' invasion that straddles two regions (e.g. France and Belgium)
#' produces one cluster, not two. Points falling outside all regions
#' are tagged as \code{"Other"} and excluded from DBSCAN.
#'
#' WKT polygons are produced via per-point buffering (10 km) followed
#' by union and adaptive simplification (1 km tolerance, increased
#' until the polygon has at most 10000 vertices and the WKT string is
#' at most 3000 characters — the practical limit for GBIF geometry
#' predicates).
#'
#' @examples
#' \dontrun{
#' # Assuming presence_xy and mean_r are already in memory:
#' pa <- dc_generate_pseudo_absences(presence_xy, mean_r,
#'                                   pa_ratio = 3)
#' # Submit one GBIF download per cluster for effort extraction:
#' for (cid in names(pa$cluster_wkt_list)) {
#'   key <- rgbif::occ_download(
#'     rgbif::pred("hasCoordinate", TRUE),
#'     rgbif::pred("geometry", pa$cluster_wkt_list[[cid]]),
#'     rgbif::pred_gte("year", 1800),
#'     rgbif::pred_lte("year", 2026),
#'     format = "SIMPLE_CSV")
#' }
#' }
#'
#' @importFrom data.table data.table as.data.table rbindlist uniqueN
#' @importFrom terra rast cellFromXY xyFromCell values vect rasterize
#' @importFrom terra project buffer aggregate crs
#' @importFrom sf st_as_sf st_as_sfc st_buffer st_convex_hull st_union
#' @importFrom sf st_join st_within st_make_valid st_simplify st_coordinates
#' @importFrom sf st_geometry st_as_text st_bbox
#' @importFrom dbscan dbscan kNNdist
#' @importFrom countrycode countrycode
#' @importFrom stats kmeans
#' @export
dc_generate_pseudo_absences <- function(presence_xy,
                                        raster_template,
                                        pa_ratio            = 3,
                                        target_cluster_size = 20,
                                        dbscan_eps          = 3,
                                        dbscan_min_pts      = 5,
                                        hull_buffer_deg     = 2,
                                        exclusion_radius_km = 50,
                                        seed                = 42) {

  decimalLongitude <- decimalLatitude <- cell_id <- climate_class <- NULL
  db_cluster <- cluster_key <- cluster_id <- continent <- NULL

  if (!data.table::is.data.table(presence_xy))
    presence_xy <- data.table::as.data.table(presence_xy)

  presence_cells_pa <- unique(
    terra::cellFromXY(raster_template,
                      as.matrix(presence_xy[, .(decimalLongitude,
                                                decimalLatitude)])))
  presence_cells_pa <- presence_cells_pa[!is.na(presence_cells_pa)]

  # ── 1. Region definitions ────────────────────────────────────
  REGIONS <- list(
    NorthAmerica    = sf::st_as_sfc(
      "POLYGON((-180 14,-52 14,-52 84,-180 84,-180 14))", crs = 4326),
    SouthAmerica    = sf::st_as_sfc(
      "POLYGON((-82 -56,-34 -56,-34 14,-82 14,-82 -56))",  crs = 4326),
    Europe          = sf::st_as_sfc(
      "POLYGON((-25 34,45 34,45 72,-25 72,-25 34))",        crs = 4326),
    RussiaNorthAsia = sf::st_as_sfc(
      "POLYGON((45 50,180 50,180 78,45 78,45 50))",          crs = 4326),
    WestCentralAsia = sf::st_as_sfc(
      "POLYGON((25 12,75 12,75 50,25 50,25 12))",            crs = 4326),
    EastAsia        = sf::st_as_sfc(
      "POLYGON((100 18,145 18,145 55,100 55,100 18))",       crs = 4326),
    SouthEastAsia   = sf::st_as_sfc(
      "POLYGON((92 -11,142 -11,142 28,92 28,92 -11))",      crs = 4326),
    Oceania         = sf::st_as_sfc(
      "POLYGON((110 -55,180 -55,180 0,110 0,110 -55))",     crs = 4326),
    AfricaMidEast   = sf::st_as_sfc(
      "POLYGON((-20 -35,62 -35,62 38,-20 38,-20 -35))",    crs = 4326)
  )
  regions_sf <- do.call(rbind, lapply(names(REGIONS), function(nm)
    sf::st_sf(region = nm, geometry = REGIONS[[nm]])))

  # ── 2. Tag presence points with region ───────────────────────
  pts_sf_pa <- sf::st_as_sf(terra::vect(presence_xy,
    geom = c("decimalLongitude", "decimalLatitude"),
    crs  = "EPSG:4326"))
  pts_sf_pa <- sf::st_join(pts_sf_pa, regions_sf,
                           join = sf::st_within, left = TRUE)
  pts_sf_pa$region[is.na(pts_sf_pa$region)] <- "Other"

  # ── 3. Global DBSCAN ─────────────────────────────────────────
  pts_sf_valid <- pts_sf_pa[pts_sf_pa$region != "Other", ]
  coords_all   <- sf::st_coordinates(pts_sf_valid)

  eps_val <- if (!is.null(dbscan_eps)) {
    dbscan_eps
  } else {
    if (nrow(coords_all) <= dbscan_min_pts) NA_real_ else {
      kd_s <- sort(dbscan::kNNdist(coords_all, k = dbscan_min_pts))
      d2   <- diff(diff(kd_s))
      kd_s[which.max(d2) + 2L]
    }
  }
  message("Global DBSCAN: eps = ", round(eps_val, 3),
          " | minPts = ", dbscan_min_pts)

  db_global <- dbscan::dbscan(coords_all, eps = eps_val,
                              minPts = dbscan_min_pts)
  pts_sf_valid$db_cluster <- db_global$cluster
  message("Clusters: ", max(db_global$cluster),
          " | Noise: ", sum(db_global$cluster == 0))

  # ── Noise → external validation ──────────────────────────────
  noise_pts <- pts_sf_valid[pts_sf_valid$db_cluster == 0L, ]
  if (nrow(noise_pts) > 0) {
    nc <- sf::st_coordinates(noise_pts)
    noise_df <- data.table::data.table(
      decimalLongitude = nc[, "X"],
      decimalLatitude  = nc[, "Y"],
      db_cluster       = 0L,
      label            = "noise_presence",
      cell_id          = terra::cellFromXY(raster_template, nc))
    noise_df[, climate_class := as.integer(
      terra::values(raster_template)[cell_id])]
    noise_df[, region := noise_pts$region]
  } else {
    noise_df <- data.table::data.table()
  }

  # ── Helper: buffered convex hull ─────────────────────────────
  make_cluster_hull <- function(pts_mat, buffer_deg) {
    pts_df <- as.data.frame(pts_mat)
    if (!all(c("X", "Y") %in% names(pts_df)))
      colnames(pts_df) <- c("X", "Y")
    pts_tmp <- sf::st_as_sf(pts_df, coords = c("X", "Y"), crs = 4326)
    if (nrow(pts_df) < 3)
      return(sf::st_buffer(sf::st_as_sfc(sf::st_bbox(pts_tmp)), buffer_deg))
    sf::st_buffer(sf::st_convex_hull(sf::st_union(pts_tmp)), buffer_deg)
  }

  # ── 4. Per-cluster sampling ──────────────────────────────────
  pseudo_list        <- list()
  cluster_hulls_list <- list()
  valid_clusters     <- sort(unique(
    pts_sf_valid$db_cluster[pts_sf_valid$db_cluster > 0]))

  set.seed(seed)

  for (cl in valid_clusters) {
    cl_key <- paste0("global_cl", cl)
    cl_pts <- pts_sf_valid[pts_sf_valid$db_cluster == cl, ]
    cl_mat <- sf::st_coordinates(cl_pts)

    cl_presence_cells <- unique(terra::cellFromXY(raster_template, cl_mat))
    cl_presence_cells <- cl_presence_cells[!is.na(cl_presence_cells)]
    n_cl_cells        <- length(cl_presence_cells)
    message("\n── Cluster ", cl, ": ", nrow(cl_pts),
            " pts | ", n_cl_cells, " cells ──")

    if (n_cl_cells < 3) { message("  Skipping."); next }

    cl_hull <- make_cluster_hull(cl_mat, hull_buffer_deg)
    cluster_hulls_list[[cl_key]] <- cl_hull

    hull_vect      <- terra::vect(cl_hull)
    terra::crs(hull_vect) <- "EPSG:4326"
    hull_rast      <- terra::rasterize(hull_vect, raster_template, field = 1)
    hull_cells     <- which(!is.na(terra::values(hull_rast)))

    candidate_cells <- setdiff(hull_cells, presence_cells_pa)
    candidate_cells <- candidate_cells[!is.na(
      terra::values(raster_template)[candidate_cells])]
    candidate_cells <- candidate_cells[
      terra::values(raster_template)[candidate_cells] != 3L]

    if (exclusion_radius_km > 0 && nrow(cl_pts) > 0) {
      excl_buf <- terra::buffer(
        terra::project(terra::vect(cl_pts), "EPSG:3857"),
        width = exclusion_radius_km * 1000)
      excl_buf  <- terra::project(terra::aggregate(excl_buf), "EPSG:4326")
      excl_rast <- terra::rasterize(excl_buf, raster_template, background = 0)
      candidate_cells <- setdiff(candidate_cells,
                                 which(terra::values(excl_rast) == 1))
    }

    message("  Candidate cells: ", length(candidate_cells))
    if (length(candidate_cells) == 0) { message("  Skipping."); next }

    n_sample      <- min(n_cl_cells * pa_ratio, length(candidate_cells))
    sampled_cells <- sample(candidate_cells, n_sample, replace = FALSE)
    xy_cl         <- data.table::as.data.table(
      terra::xyFromCell(raster_template, sampled_cells))
    data.table::setnames(xy_cl, c("decimalLongitude", "decimalLatitude"))

    pseudo_list[[cl_key]] <- data.table::data.table(
      cell_id          = sampled_cells,
      db_cluster       = cl,
      cluster_key      = cl_key,
      decimalLongitude = xy_cl$decimalLongitude,
      decimalLatitude  = xy_cl$decimalLatitude)
  }

  # ── 5. Combine + sanity checks ───────────────────────────────
  pseudo_dt <- data.table::rbindlist(pseudo_list, fill = TRUE)
  stopifnot(length(intersect(pseudo_dt$cell_id, presence_cells_pa)) == 0)
  pseudo_dt[, climate_class := as.integer(
    terra::values(raster_template)[cell_id])]
  stopifnot(nrow(pseudo_dt[climate_class == 3L]) == 0)
  message("\nTotal pseudo-absences: ", nrow(pseudo_dt))

  # ── 6. Continent assignment ──────────────────────────────────
  world_v <- geodata::world(resolution = 1, path = tempdir())
  world_sf_cont <- sf::st_as_sf(world_v)
  world_sf_cont$continent <- suppressWarnings(
    countrycode::countrycode(world_sf_cont$GID_0, "iso3c", "continent"))
  world_sf_cont$continent[is.na(world_sf_cont$continent)] <- "Unknown"
  world_sf_cont <- sf::st_make_valid(
    world_sf_cont[, c("continent", "geometry")])

  pseudo_sf_pts <- sf::st_as_sf(pseudo_dt,
    coords = c("decimalLongitude", "decimalLatitude"), crs = 4326)
  pseudo_dt[, continent := sf::st_join(pseudo_sf_pts, world_sf_cont,
                                       join = sf::st_within)$continent]
  pseudo_dt[is.na(continent), continent := "Unknown"]

  # ── 7. k-means batching for GBIF downloads ───────────────────
  set.seed(seed)
  pseudo_dt[, cluster_id := NA_character_]
  for (cont in pseudo_dt[continent != "Unknown", unique(continent)]) {
    idx <- which(pseudo_dt$continent == cont)
    k   <- max(1L, ceiling(length(idx) / target_cluster_size))
    if (k == 1L) {
      pseudo_dt[idx, cluster_id := paste0(cont, "_1")]
    } else {
      km <- stats::kmeans(
        as.matrix(pseudo_dt[idx, .(decimalLongitude, decimalLatitude)]),
        centers = k, nstart = 10, iter.max = 50)
      pseudo_dt[idx, cluster_id := paste0(cont, "_", km$cluster)]
    }
  }
  unk_idx <- which(pseudo_dt$continent == "Unknown")
  if (length(unk_idx) > 0)
    pseudo_dt[unk_idx, cluster_id := paste0("Unknown_", seq_len(.N))]

  # ── 8. Build WKT polygons per k-means batch ──────────────────
  cluster_ids      <- pseudo_dt[continent != "Unknown", unique(cluster_id)]
  cluster_wkt_list <- list()

  for (cid in cluster_ids) {
    cells <- pseudo_dt[cluster_id == cid, cell_id]
    xy    <- terra::xyFromCell(raster_template, cells)
    pts_m <- terra::project(terra::vect(xy, crs = "EPSG:4326"),
                            "EPSG:3857")
    union_wgs <- terra::project(
      terra::aggregate(terra::buffer(pts_m, width = 10000)),
      "EPSG:4326")
    cl_sf <- sf::st_make_valid(sf::st_as_sf(union_wgs))

    tol <- 1000; step <- 500; max_tol <- 50000
    repeat {
      cl_s  <- sf::st_simplify(cl_sf, dTolerance = tol,
                               preserveTopology = TRUE)
      wkt_c <- sf::st_as_text(sf::st_geometry(cl_s))
      if (nrow(sf::st_coordinates(cl_s)) <= 10000 && nchar(wkt_c) <= 3000) {
        break
      }
      tol <- tol + step
      if (tol > max_tol) {
        cl_s  <- sf::st_as_sf(sf::st_as_sfc(sf::st_bbox(cl_sf)))
        wkt_c <- sf::st_as_text(sf::st_geometry(cl_s))
        message("  ", cid, " bbox fallback"); break
      }
    }
    cluster_wkt_list[[cid]] <- wkt_c
  }

  list(pseudo_dt        = pseudo_dt,
       cluster_hulls    = cluster_hulls_list,
       noise_presences  = noise_df,
       cluster_wkt_list = cluster_wkt_list)
}
