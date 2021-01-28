#' Retrieve the shapefile for a country at the given level.
#' GADM identifies countries with the Alpha-3 country code.
#' If the file doesn't exist, it will return a NULL.
#' @param alpha3 The three-letter country code.
#' @param admin_level 0 is the whole country. Different countries have more or
#'     fewer levels.
#' @param gadm_version When this version of GADM was downloaded.
#' @return A shapefile from the `sf` package, read with `st_read`.
#'
#' This function depends on having a working `tempfile()` command because it
#' has to unzip the shapefile into a temporary directory.
#'
#' @export
gadm_country_shapefile <- function(alpha3, admin_level = 0, gadm_version = "201104") {
    gadm_rp <- rampdata::ramp_path(
        "/inputs/gadmshape",
        version = gadm_version
        )
    pfn <- rampdata::as.path(gadm_rp)
    if (!dir.exists(pfn)) {
        flog.error(paste("The directory for the gadm file doesn't exist", pfn))
        return(NULL)
    }
    matching <- list.files(pfn, paste0("*_", toupper(alpha3), "_*"))
    if (length(matching) == 1) {
        zip_rp <- rampdata::add_path(gadm_rp, file = matching)
        inner <- utils::unzip(rampdata::as.path(zip_rp), list = TRUE)
        outline_idx <- grep(paste0("_", admin_level, "."), inner$Name)
        if (length(outline_idx) > 0) {
            inner_fn <- inner$Name[outline_idx]
            unzip_dir <- tempfile()
            country_sf <- tryCatch({
                utils::unzip(rampdata::as.path(zip_rp), files = inner_fn, exdir = unzip_dir)
                utils::capture.output({shape <- sf::st_read(unzip_dir)})
                shape
            },
            finally = {
                if (dir.exists(unzip_dir)) {
                    unlink(unzip_dir, recursive = TRUE)
                }
            })
        } else {
            flog.error(paste("The country", alpha3,
                    "doesn't have a gadm outline at level", admin_level))
            return(NULL)
        }
    } else {
        flog.error(paste("The gadm file with country code", alpha3, "not found"))
        return(NULL)
    }
    country_sf
}
