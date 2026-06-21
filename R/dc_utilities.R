#' Taxonomic relevance configuration
#'
#' Returns the taxonomic relevance specification for one of the
#' supported recorder communities. The returned object is a list with
#' \code{phylum}, \code{orders}, and/or \code{exclude_phyla} elements
#' that is passed to \code{dc_extract_effort} via the \code{cfg}
#' argument. Records matching the configuration are tagged
#' \code{is_taxon_relevant = TRUE}; this drives the
#' \code{relevance_ratio} and \code{recorder_relevance_ratio}
#' predictors.
#'
#' "Taxonomic relevance" describes \emph{what was recorded in a cell},
#' not recorder competence. A cell containing many beetle records has
#' a high relevance ratio for a beetle focal species, regardless of the
#' formal qualifications of the recorders.
#'
#' For vertebrates the configuration is keyed by class via Linnaean
#' orders. For plants and invertebrates it is keyed by recording
#' practice (terrestrial vascular plants vs. bryophytes vs. marine,
#' for example) because plant systematics does not map cleanly to a
#' single rank that recorders work at.
#'
#' @param group Character, one of \code{"fish"}, \code{"mammal"},
#'   \code{"reptile"}, \code{"anuran"}, \code{"bird"},
#'   \code{"terrestrial_vascular_plants"}, \code{"bryophytes"},
#'   \code{"terrestrial_invertebrates"},
#'   \code{"freshwater_invertebrates"}, or \code{"unsupported"}.
#'
#' @return A named list with elements \code{domain},
#'   \code{phylum} or \code{kingdom}, and optionally \code{orders}
#'   and/or \code{exclude_phyla}.
#'
#' @examples
#' cfg <- dc_expertise_config("terrestrial_vascular_plants")
#' cfg
#' @export
dc_expertise_config <- function(group) {

  configs <- list(
    fish = list(
      domain = "vertebrates", phylum = "Chordata",
      orders = c("Cypriniformes", "Siluriformes", "Salmoniformes",
                 "Cyprinodontiformes", "Cichliformes", "Perciformes",
                 "Gobiiformes", "Clupeiformes", "Esociformes",
                 "Anguilliformes", "Osmeriformes", "Gadiformes",
                 "Atheriniformes", "Beloniformes", "Carcharhiniformes",
                 "Lamniformes", "Rajiformes", "Myliobatiformes")
    ),
    mammal = list(
      domain = "vertebrates", phylum = "Chordata",
      orders = c("Artiodactyla", "Perissodactyla", "Carnivora", "Rodentia",
                 "Lagomorpha", "Primates", "Chiroptera", "Proboscidea",
                 "Hyracoidea", "Sirenia", "Pilosa", "Cingulata",
                 "Didelphimorphia", "Dasyuromorphia", "Peramelemorphia",
                 "Diprotodontia", "Monotremata", "Eulipotyphla")
    ),
    reptile = list(
      domain = "vertebrates", phylum = "Chordata",
      orders = c("Squamata", "Testudines", "Crocodylia", "Rhynchocephalia")
    ),
    anuran = list(
      domain = "vertebrates", phylum = "Chordata",
      orders = "Anura"
    ),
    bird = list(
      domain = "vertebrates", phylum = "Chordata",
      orders = c("Passeriformes", "Anseriformes", "Galliformes",
                 "Accipitriformes", "Falconiformes", "Strigiformes",
                 "Columbiformes", "Charadriiformes", "Pelecaniformes",
                 "Ciconiiformes", "Gruiformes", "Procellariiformes",
                 "Podicipediformes", "Gaviiformes", "Apodiformes",
                 "Piciformes", "Coraciiformes", "Psittaciformes",
                 "Struthioniformes", "Tinamiformes")
    ),
    terrestrial_vascular_plants = list(
      domain = "plants", kingdom = "Plantae",
      exclude_phyla = c("Bryophyta", "Marchantiophyta",
                       "Anthocerotophyta")
    ),
    bryophytes = list(
      domain = "plants", kingdom = "Plantae",
      phylum = c("Bryophyta", "Marchantiophyta", "Anthocerotophyta")
    ),
    terrestrial_invertebrates = list(
      domain = "invertebrates",
      phylum = c("Arthropoda", "Mollusca", "Annelida"),
      orders = c("Coleoptera", "Lepidoptera", "Diptera", "Hymenoptera",
                 "Hemiptera", "Orthoptera", "Odonata", "Mantodea",
                 "Phasmatodea", "Blattodea", "Araneae", "Opiliones",
                 "Scorpiones", "Chilopoda", "Diplopoda")
    ),
    freshwater_invertebrates = list(
      domain = "invertebrates",
      phylum = c("Arthropoda", "Mollusca", "Annelida"),
      orders = c("Unionida", "Venerida", "Basommatophora",
                 "Ephemeroptera", "Trichoptera", "Plecoptera", "Odonata",
                 "Decapoda", "Amphipoda", "Isopoda"),
      freshwater_only = TRUE
    ),
    unsupported = list(domain = "unsupported")
  )

  if (!group %in% names(configs))
    stop("Unknown expertise group: '", group,
         "'. Supported groups: ",
         paste(names(configs), collapse = ", "))

  configs[[group]]
}


#' Build a buffered WKT polygon from presence points
#'
#' Buffers a set of presence points and returns a simplified WKT polygon
#' suitable for use as a GBIF geometry predicate. The simplification
#' tolerance is increased adaptively until the polygon has at most
#' \code{max_vertices} vertices and the WKT string is at most
#' \code{max_chars} characters (GBIF's practical limits).
#'
#' @param presence_xy A data.frame or data.table with columns
#'   \code{decimalLongitude} and \code{decimalLatitude}.
#' @param buffer_m Buffer in metres around each point. Default
#'   \code{15000} (15 km).
#' @param max_vertices Maximum vertices in the output polygon.
#'   Default \code{10000}.
#' @param max_chars Maximum WKT string length. Default \code{10000}.
#' @param initial_tol Initial simplification tolerance in metres.
#'   Default \code{1000}.
#' @param step Tolerance increment in metres. Default \code{500}.
#' @param max_tol Maximum allowed tolerance before giving up. Default
#'   \code{10000}.
#'
#' @return A WKT string suitable for
#'   \code{rgbif::pred("geometry", wkt)}.
#'
#' @importFrom terra vect project buffer aggregate crs
#' @importFrom sf st_as_sf st_make_valid st_simplify st_as_text
#' @importFrom sf st_geometry st_coordinates
#' @export
dc_buffered_wkt <- function(presence_xy,
                            buffer_m     = 15000,
                            max_vertices = 10000,
                            max_chars    = 10000,
                            initial_tol  = 1000,
                            step         = 500,
                            max_tol      = 10000) {

  pres_v   <- terra::vect(presence_xy, crs = "EPSG:4326")
  pres_m   <- terra::project(pres_v, "EPSG:3857")
  buffers  <- terra::buffer(pres_m, width = buffer_m)
  buf_wgs  <- terra::project(terra::aggregate(buffers), "EPSG:4326")
  buf_sf   <- sf::st_make_valid(sf::st_as_sf(buf_wgs))

  tol <- initial_tol
  repeat {
    simp <- sf::st_simplify(buf_sf, dTolerance = tol,
                            preserveTopology = TRUE)
    n_v  <- nrow(sf::st_coordinates(simp))
    wkt  <- sf::st_as_text(sf::st_geometry(simp))
    if (n_v <= max_vertices && nchar(wkt) <= max_chars) return(wkt)
    tol <- tol + step
    if (tol > max_tol)
      stop("Cannot reduce geometry below limits ",
           "(max vertices = ", max_vertices,
           ", max chars = ", max_chars, ").")
  }
}
