v


#' Given a bunch of filenames, find all four-number digit strings.
#' @param filenames A list of character strings to search.
#' @return A list of numeric years.
years_in_filenames <- function(filenames) {
    matches <- regexpr("[0-9]{4}", filenames)
    start <- matches[matches > 0]
    stop <- attr(matches, "match.length")[matches > 0] + start - 1
    year_files <- filenames[matches > 0]
    years <- substr(year_files, start, stop)
    names(year_files) <- years
    year_files
}


#' Given a bounding box in map coordinates (maybe lat-long), return pixel bbox.
#' @param map_raster The `raster` package raster for the pixels.
#' @param bbox A bounding box as returned by `sf::st_bbox`.
#' @return An extent, in pixels.
pixel_bounding_box <- function(map_raster, bbox) {
    x <- c(bbox["xmin"], bbox["xmax"])
    y <- c(bbox["ymin"], bbox["ymax"])
    r <- integer(2)
    c <- integer(2)
    for (i in 1:2) {
        c[i] <- raster::colFromX(map_raster, x[i])
        r[i] <- raster::rowFromY(map_raster, y[i])
    }
    # The row is y and row 1 has a greater y than row 2.
    c(rmin = min(r), rmax = max(r), cmin = min(c), cmax = max(c))
}


available_data <- function(input_version, country_code, select_years) {
    pfpr_dir <- rampdata::workflow_path("pfpr")
    pfpr_files <- list.files(rampdata::as.path(pfpr_dir), "*.tif")
    pfpr_years <- as.integer(names(years_in_filenames(pfpr_files)))
    am_dir <- rampdata::workflow_path("am")
    am_files <- list.files(rampdata::as.path(am_dir), "*.tif")
    am_years <- as.integer(names(years_in_filenames(am_files)))
    # We could choose, instead, to extend either one to cover missing years
    # in the other.
    shared_years <- sort(intersect(am_years, pfpr_years))
    if (!is.null(select_years)) {
        shared_years <- sort(intersect(shared_years, select_years))
    }

    pr_fn <- pfpr_files[1]
    pfpr_file <- rampdata::add_path(pfpr_dir, file = pr_fn)
    flog.info(paste("reading a sample slice for this work from", pfpr_file))
    pfpr <- raster::raster(rampdata::as.path(pfpr_file))
    whole_input_extent <- c(rmin = 1, rmax = raster::nrow(pfpr),
        cmin = 1, cmax = raster::ncol(pfpr)
    )
    if (!is.null(country_code)) {
        outline_sf <- gadm_country_shapefile(country_code)
        bbox <- sf::st_bbox(outline_sf)
        domain_extent <- pixel_bounding_box(pfpr, bbox)
        flog.info(paste("using country", country_code, "to limit extent to",
            paste(domain_extent, collapse = ",")))
        stopifnot(domain_extent["rmin"] < domain_extent["rmax"])
        stopifnot(domain_extent["cmin"] < domain_extent["cmax"])
    } else {
        domain_extent <- c(
            rmin = 1, rmax = raster::nrow(pfpr), cmin = 1, cmax = raster::ncol(pfpr))
    }
    xres <- raster::xres(pfpr)
    yres <- raster::yres(pfpr)
    domain_dimensions <- list(
        xmin = raster::xmin(pfpr) + (domain_extent["cmin"] - 1) * xres,
        xmax = raster::xmin(pfpr) + domain_extent["cmax"] * xres,
        xres = xres,
        ymax = raster::ymax(pfpr) - (domain_extent["rmin"] - 1) * yres,
        ymin = raster::ymax(pfpr) - domain_extent["rmax"] * yres,
        yres = yres,
        projection = raster::projection(pfpr)
    )
    flog.debug(paste("domain extent", paste(domain_extent, collapse = ",")))
    flog.debug(paste("domain dimensions", paste(domain_dimensions, collapse = ",")))
    # Get a sample so we can find NA values later.
    nrows <- domain_extent["rmax"] - domain_extent["rmin"] + 1
    ncols <- domain_extent["cmax"] - domain_extent["cmin"] + 1
    sample <- matrix(
        raster::getValuesBlock(
            pfpr,
            row = domain_extent["rmin"], nrows = nrows,
            col = domain_extent["cmin"], ncols = ncols
        ),
        nrow = nrows,
        byrow = TRUE  # This is important. Raster returns by in C order.
    )
    stopifnot(domain_dimensions$xmin < domain_dimensions$xmax)
    stopifnot(domain_dimensions$ymin < domain_dimensions$ymax)
    list(
        years = shared_years,
        whole_input_extent = whole_input_extent,
        domain_extent = domain_extent,
        domain_dimensions = domain_dimensions,
        sample_slice = sample,
        draws = 1
    )
}


#' Get all the data for a country as a convenience function.
#' @param country_alpha3 This is the three-letter country code, as in "uga".
#' @param years A vector of integer years, as in `2000:2020`.
data_for_country <- function(country_alpha3, years) {
    # The country will be our lat-long outline.
    outline_sf <- gadm_country_shapefile(country_alpha3)
    bbox <- sf::st_bbox(outline_sf)

    pfpr_dir <- rampdata::workflow_path("pfpr")
    am_dir <- rampdata::workflow_path("am")
    pfpr_yearly <- years_in_filenames(list.files(rampdata::as.path(pfpr_dir)))
    am_yearly <- years_in_filenames(list.files(rampdata::as.path(am_dir)))

    pfpr_file <- rampdata::add_path(pfpr_dir, file = pfpr_yearly[1])
    pfpr <- raster::raster(rampdata::as.path(pfpr_file))
    domain_extent <- pixel_bounding_box(pfpr, bbox)
    flog.debug("domain_extent", paste(domain_extent, collapse = ","))

    pfpr_all <- list()
    am_all <- list()
    for (year_idx in 1:length(years)) {
        year <- years[year_idx]
        pr_fn <- pfpr_yearly[as.character(year)]
        pfpr_file <- rampdata::add_path(pfpr_dir, file = pr_fn)
        if (!file.exists(rampdata::as.path(pfpr_file))) {
            msg <- paste("Cannot find pfpr at", rampdata::as.path(pfpr_file))
            flog.error(msg)
            stop(msg)
        }

        am_fn <- am_yearly[as.character(year)]
        am_file <- rampdata::add_path(am_dir, file = am_fn)
        if (!file.exists(rampdata::as.path(am_file))) {
            msg <- paste("Cannot find pfpr at", rampdata::as.path(am_file))
            flog.error(msg)
            stop(msg)
        }

        pfpr <- raster::raster(rampdata::as.path(pfpr_file))
        am <- raster::raster(rampdata::as.path(am_file))

        load_extent <- c(
            domain_extent["rmin"],
            rmax = min(domain_extent["rmax"], raster::nrow(pfpr)),
            domain_extent["cmin"],
            cmax = min(domain_extent["cmax"], raster::ncol(pfpr))
        )
        # We use byrow because raster returns values in C order, not
        # Fortran order. This gives us the same definition of row and column
        # as the one that raster uses.
        nrows <- load_extent["rmax"] - load_extent["rmin"] + 1
        ncols <- load_extent["cmax"] - load_extent["cmin"] + 1
        flog.debug("load_extent", paste(load_extent, collapse = ","))
        prblock <- matrix(
            raster::getValuesBlock(
                pfpr, row = load_extent["rmin"], nrows = nrows, col = load_extent["cmin"], ncols = ncols),
            nrow = nrows,
            byrow = TRUE
        )
        if (!all(load_extent == domain_extent)) {
            pr_outer <- matrix(
                nrow = domain_extent["rmax"] - domain_extent["rmin"] + 1,
                ncol = domain_extent["cmax"] - domain_extent["cmin"] + 1)
            pr_outer[1:nrows, 1:ncols] <- prblock
        } else {
            pr_outer <- prblock
        }

        pfpr_all[[as.character(year)]] <- pr_outer
        amblock <- matrix(
            raster::getValuesBlock(
                am, row = load_extent["rmin"], nrows = nrows, col = load_extent["cmin"], ncols = ncols),
            nrow = nrows,
            byrow = TRUE
        )
        if (!all(load_extent == domain_extent)) {
            am_outer <- matrix(
                nrow = domain_extent["rmax"] - domain_extent["rmin"] + 1,
                ncol = domain_extent["cmax"] - domain_extent["cmin"] + 1)
            am_outer[1:nrows, 1:ncols] <- amblock
        } else {
            am_outer <- amblock
        }
        am_all[[as.character(year)]] <- am_outer
    }

    combine <- function(list_of_matrices) {
        sentinel <- dim(list_of_matrices[[1]])
        array(do.call(c, list_of_matrices), dim = c(sentinel, length(list_of_matrices)))
    }
    pfpr_arr <- combine(pfpr_all)
    am_arr <- combine(am_all)
    stopifnot(length(dim(pfpr_arr)) == 3)
    stopifnot(length(dim(am_arr)) == 3)
    list(
        pfpr = pfpr_arr,
        am = am_arr
    )
}


#' Load input data, including pr2ar, the PfPR, and the AM.
#' @param config Input configuration file name.
#' @param pr2ar_version Which version of the `pr2ar_mesh.csv` file.
#' @param domain_extent Bounding box that we compute within the MAP raster.
#' @param years Which years to do.
#' @return A list with the loaded data.
load_data <- function(config, pr2ar_version, domain_extent, years) {
    parameters <- configr::read.config(config)[["parameters"]]
    pr2ar_rp <- rampdata::workflow_path("pr2ar")
    pr2ar_fn <- rampdata::as.path(pr2ar_rp)
    rampdata::prov.input.file(pr2ar_fn, "pr2ar")
    pr_to_ar_dt <- data.table::fread(pr2ar_fn)
    pfpr_dir <- rampdata::workflow_path("pfpr")
    am_dir <- rampdata::workflow_path("am")
    pfpr_yearly <- years_in_filenames(list.files(rampdata::as.path(pfpr_dir)))
    am_yearly <- years_in_filenames(list.files(rampdata::as.path(am_dir)))

    pfpr_all <- list()
    am_all <- list()
    for (year_idx in 1:length(years)) {
        year <- years[year_idx]
        year_str <- as.character(year)
        if (!year_str %in% names(pfpr_yearly)) {
          msg <- sprintf("Cannot find year %d in %s", year,
                             paste0(names(pfpr_yearly), collapse = ","))
          flog.error(msg)
          stop(msg)
        }
        pr_fn <- pfpr_yearly[year_str]
        pfpr_file <- rampdata::add_path(pfpr_dir, file = pr_fn)
        if (!file.exists(rampdata::as.path(pfpr_file))) {
            msg <- paste("Cannot find pfpr at", rampdata::as.path(pfpr_file))
            flog.error(msg)
            stop(msg)
        }

        am_fn <- am_yearly[as.character(year)]
        am_file <- rampdata::add_path(am_dir, file = am_fn)
        if (!file.exists(rampdata::as.path(am_file))) {
            msg <- paste("Cannot find pfpr at", rampdata::as.path(am_file))
            flog.error(msg)
            stop(msg)
        }

        pfpr_full_fn <- rampdata::as.path(pfpr_file)
        rampdata::prov.input.file(pfpr_full_fn, "pfpr_input")
        pfpr <- raster::raster(pfpr_full_fn)
        am_full_fn <- rampdata::as.path(am_file)
        rampdata::prov.input.file(am_full_fn, "am_input")
        am <- raster::raster(am_full_fn)

        load_extent <- c(
            domain_extent["rmin"],
            rmax = min(domain_extent["rmax"], raster::nrow(pfpr)),
            domain_extent["cmin"],
            cmax = min(domain_extent["cmax"], raster::ncol(pfpr))
        )
        stopifnot(names(load_extent) == c("rmin", "rmax", "cmin", "cmax"))
        # We use byrow because raster returns values in C order, not
        # Fortran order. This gives us the same definition of row and column
        # as the one that raster uses.
        nrows <- load_extent["rmax"] - load_extent["rmin"] + 1
        ncols <- load_extent["cmax"] - load_extent["cmin"] + 1
        prblock <- matrix(
            raster::getValuesBlock(
                pfpr, row = load_extent["rmin"], nrows = nrows, col = load_extent["cmin"], ncols = ncols),
            nrow = nrows,
            byrow = TRUE
        )
        if (all(load_extent == domain_extent)) {
            pfpr_tiles <- prblock
        } else {
            pfpr_tiles <- matrix(
                nrow = domain_extent["rmax"] - domain_extent["rmin"] + 1,
                ncol = domain_extent["cmax"] - domain_extent["cmin"] + 1)
            pfpr_tiles[1:nrows, 1:ncols] <- prblock
        }
        pfpr_all[[as.character(year)]] <- pfpr_tiles
        amblock <- matrix(
            raster::getValuesBlock(
                am, row = load_extent["rmin"], nrows = nrows, col = load_extent["cmin"], ncols = ncols),
            nrow = nrows,
            byrow = TRUE
        )
        if (all(load_extent == domain_extent)) {
            am_tiles <- amblock
        } else {
            am_tiles <- matrix(
                nrow = domain_extent["rmax"] - domain_extent["rmin"] + 1,
                ncol = domain_extent["cmax"] - domain_extent["cmin"] + 1)
            am_tiles[1:nrows, 1:ncols] <- amblock
        }
        am_all[[as.character(year)]] <- am_tiles
    }
    list(
        parameters = parameters,
        pr_to_ar_dt = pr_to_ar_dt,
        years = years,
        pfpr = pfpr_all,
        am = am_all,
        offset = c(row = unname(domain_extent["rmin"]), col = unname(domain_extent["cmin"]))
        )
}

