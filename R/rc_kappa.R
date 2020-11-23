#' This script reads PfPR and treatment rasters and produces several variables.
#'
#' The goal is to get R_c from PfPR and some measure of treatment.
#' Input data is from the Malaria Atlas Project (MAP). It comes as GeoTIFFs
#' for multiple years, for all of Africa.
#'
#' This script will work on a subset of input data, a subset of both
#' lat-long pixels and years. Most of the work in this script is meant to
#' process the data in parallel while preserving the clarity of the
#' central calculation on the time series of each pixel through the
#' years. That central calculation is in `pixel_work`.
#'
#' We are setting up to work on time series through the months and years.
#' A time series will look like a 3D tile (time, spatial rows, spatial cols).
#' There can be draws for these, too. Most of the complexity in this code
#' is handling the separation of work into tiles and reconstituting it.
#' The tools here will enable a splitting across jobs on the cluster and
#' splitting into parallel processes within a job.
#'
#' There are subdivisions of the raster, each a bounding box we call an extent.
#'   1. `whole_input_extent` - The input MAP data dimensions.
#'   2. `domain_extent` - The rectangle for the calculation (say, Uganda), within
#'          the `whole_input_extent`.
#'   3. `process_extent` - The part that is loaded into this process on the cluster.
#'          This is relative to the domain extent and calculated by looking at
#'          which tiles this process will compute.
#'
#' The domain extent is everything we will calculate in a single cluster job.
#' The process extent is what this particular cluster task will calculate.
#' We split the `domain_extent` into tiles, and each process will
#' calculate some set of tiles.
#'
#' What to read first:
#'
#' * `main` is an entrypoint into the program from the command line.
#' * `funcmain` is what you would call to run this from within R.
#'   This function is a roadmap to the code.

library(futile.logger)
# Using loadNamespace so that this fails early, if it's going to fail,
# but I want to enforce the use of package identifiers in the code
# so that we know where functions come from. If this code were in a
# package, then installing the package would ensure we have the right libraries.
loadNamespace("configr")
loadNamespace("data.table")
loadNamespace("docopt")
loadNamespace("pracma")
loadNamespace("raster")
loadNamespace("rgdal")  # Not explicit but needed for raster.
# Requires units which requires libudunits2-dev.
loadNamespace("sf")
loadNamespace("sp")
loadNamespace("tmap")  # For plotting.

# remotes::install_github("dd-harp/rampdata")
loadNamespace("rampdata")


# Don't know which directory will be our working directory.
to_source <- c("country_outline.R", "hilbert.R", "plan_io.R", "serialize.R")
if (dir.exists("gen_scaled_ar")) {
  for (s in to_source) source(file.path("gen_scaled_ar", s))
} else {
  for (s in to_source) source(s)
}


#' A stack trace that shows which function had the problem.
#' https://renkun.me/2020/03/31/a-simple-way-to-show-stack-trace-on-error-in-r/
improved_errors <- function() {
    options(error = function() {
        sink(stderr())
        on.exit(sink(NULL))
        traceback(3, max.lines = 1L)
        if (!interactive()) {
            q(status = 1)
        }
    })
}


#' These are the two main data structures. Tiles and blocksizes are measured
#' in rows and columns, and they are stored as a numeric vector with names
#' "row" and "col". If you forget those, then you get NA. The same goes for
#' extents, which all have "rmin", "rmax", "cmin", "cmax" for the mininum and
#' maximum of rows and columns. We use these names to make it easier to
#' visually ascertain that rows are matched with rows and columns with columns.
is_tile <- function(p) all(names(p) == c("row", "col"))
is_extent <- function(e) all(names(e) == c("rmin", "rmax", "cmin", "cmax"))


#' Compute a domain extent from a process extent within a domain extent.
#' @param domain_extent The larger bounding box.
#' @param process_extent A bounding box within the domain extent.
#' @return A bounding box measured from the domain extent's coordinate system.
apply_relative_extent <- function(domain_extent, process_extent) {
    stopifnot(is_extent(domain_extent))
    stopifnot(is_extent(process_extent))
    offset <- domain_extent[c("rmin", "cmin")] - 1
    process_extent + unname(offset[c("rmin", "rmin", "cmin", "cmin")])
}


#' Find where the tile is within the process extent.
#' Both are relative to the domain extent.
find_relative_extent <- function(process_extent, tile_extent) {
    stopifnot(is_extent(process_extent))
    stopifnot(is_extent(tile_extent))
    offset <- process_extent[c("rmin", "cmin")] - 1
    tile_extent - unname(offset[c("rmin", "rmin", "cmin", "cmin")])
}


#' Creates a function that interpolates an AR-to-PR dataset.
#' The data calculates attack rate as a function of PR and recovery.
#' @param pr_to_ar_dt A data.table with columns AR, PR, and rho.
#' @returns A function with signature function(pr, rho) -> AR.
ar_of_pr_rho <- function(pr_to_ar_dt) {
  if (length(unique(pr_to_ar_dt$rho[1:4])) < 4) {
    # This is the route taken. flog.debug("PR moves faster in pr_to_ar mesh file.")
    dtrow <- sort(unique(pr_to_ar_dt$PR))
    dtcol <- sort(unique(pr_to_ar_dt$rho))
    stopifnot(min(dtrow) == 0)
    stopifnot(max(dtrow) == 1)
    stopifnot(min(dtcol) == 0)
    stopifnot(max(dtcol) == 1)
    stopifnot(nrow(pr_to_ar_dt) == length(dtrow) * length(dtcol))
    dtz <- array(pr_to_ar_dt$AR, dim = c(length(dtrow), length(dtcol)))
    function(pr, rho) {
      stopifnot(all(pr >= 0))
      stopifnot(all(pr <= 1))
      stopifnot(all(rho >= 0))
      stopifnot(all(rho <= 1))
      stopifnot(all(is.finite(pr)))
      stopifnot(all(is.finite(rho)))
      stopifnot(length(rho) > 0)  # This was the failure.
      stopifnot(length(pr) > 0)
      stopifnot(length(pr) == length(rho))
      # This call looks like x and y are mixed up, but check the tests below.
      tryCatch(
          pracma::interp2(dtcol, dtrow, dtz, xp = rho, yp = pr, method = "linear"),
	  error = function(e) {
              cat(paste("interp2 fail for rho", paste0(rho, collapse=",")))
              cat(paste("interp2 fail for pr", paste0(pr, collapse=",")))
	  }
      )
    }
  } else {
    flog.debug("rho moves faster in pr_to_ar mesh file.")
    dtrow <- sort(unique(pr_to_ar_dt$rho))
    dtcol <- sort(unique(pr_to_ar_dt$PR))
    stopifnot(min(dtrow) == 0)
    stopifnot(max(dtrow) == 1)
    stopifnot(min(dtcol) == 0)
    stopifnot(max(dtcol) == 1)
    stopifnot(nrow(pr_to_ar_dt) == length(dtrow) * length(dtcol))
    dtz <- array(pr_to_ar_dt$AR, dim = c(length(dtrow), length(dtcol)))
    function(pr, rho) {
      stopifnot(all(pr >= 0))
      stopifnot(all(pr <= 1))
      stopifnot(all(rho >= 0))
      stopifnot(all(rho <= 1))
      pracma::interp2(dtcol, dtrow, dtz, xp = pr, yp = rho, method = "linear")
    }
  }
}


# Keep this for testing the other one. This works for a single value.
ar_of_pr_rho2 <- function(pr_to_ar_dt) {
    dt <- pr_to_ar_dt
    function(pr, rho) {
        akima::interp(x = dt$rho, y = dt$PR, z = dt$AR, xo = rho, yo = pr,
        extrap = TRUE)[[3]]
    }
}


#' This goes the other way, from AR to PR, with rho=0.
build_ar2pr <- function(pr_ar_data) {
  no_treatment <- pr_ar_data[rho < 1e-6, c("PR", "AR")]
  # The sample data doesn't go all the way to 0 and 1, but that's the asymptotic value.
  AR_sample <- c(0, no_treatment$AR, 1)
  PR_sample <- c(0, no_treatment$PR, 1)
  function(attack_rates) {
    pracma::interp1(AR_sample, PR_sample, attack_rates, method = "linear")
  }
}


#' Parse command-line arguments.
#' @param args Optionally pass the results of `commandArgs(TRUE)`.
#'     We have this parameter so that we can test without side effects.
#' @return A list where args that aren't on the command line are NULL.
arg_parser <- function(args = NULL) {
    doc <- "pr to Rc

Usage:
  rc_kappa.R [options]
  rc_kappa.R (-h | --help)

Options:
  -h --help              Show help.
  --config=<config>      A configuration file.
  --country=<alpha3>     The three-letter country code for a country's outline.
  --outvars=<outversion> Version of output to write.
  --overwrite            Whether to overwrite outputs.
  --years=<year_range>   A range of years to do, as an R range, 2000:2010.
  --cores=<core_cnt>     Tell it how many cores to use in parallel.
  --draws=<draw_cnt>     How many draws to use.
  --task=<task_id>       If this is a task, which task.
  --tasks=<task_cnt>     Total number of tasks. You should set this for workers.
"
    if (is.null(args)) {
        args <- commandArgs(TRUE)
    }
    parsed_args <- docopt::docopt(doc, version = "rc_kappa 1.0", args = args)
    if (is.null(parsed_args$config)) {
        parsed_args$config <- "rc_kappa.toml"
    }
    parsed_args
}


#' Check whether the input arguments make sense.
#' @param args A list of command-line arguments resulting from the arg parser.
#' @return The same list of command-line arguments.
check_args <- function(args) {
    OK <- TRUE
    if (!file.exists(args$config)) {
        flog.error(paste("Cannot find file", args$config, "from the config argument.",
            "from the directory", getwd()))
        OK <- FALSE
    } else {
        if (!is.list(configr::read.config(args$config))) {
            flog.error(paste("configr says the config file has a syntax error", args$config))
            OK <- FALSE
        }
    }
    if (!is.null(args$country) && nchar(args$country) != 3) {
        flog.error(paste("The country should be the three-letter Alpha-3 code:", args$country))
        OK <- FALSE
    }
    if (!is.null(args$years)) {
        start_year <- suppressWarnings(as.integer(substr(args$years, 1, 4)))
        end_year <- suppressWarnings(as.integer(substr(args$years, 6, 9)))
        if (is.na(start_year) | is.na(end_year)) {
            flog.error(paste("Expecting the year range to be 2000:2010. Found", args$years))
            OK <- FALSE
        } else {
            args$years <- start_year:end_year
        }
    }
    # The task ID is usually not set. If it is, then we take it.
    # If it isn't, we look at the SGE_TASK_ID environment variable.
    if (!is.null(args$task)) {
      args$task <- tryCatch(as.integer(args$task), warning = function(x) {
        stop(paste("Could not convert task id to integer", args$task))
      })
    } else {
      task_id <- Sys.getenv("SGE_TASK_ID", unset = NA)
      if (!is.na(task_id)) {
        args$task <- as.integer(task_id)
      }
    }
    if (!is.null(args$tasks)) {
      args$tasks <- tryCatch(as.integer(args$tasks), warning = function(x) {
        stop(paste("Could not convert tasks cnt to integer", args$tasks))
      })
    }
    args$cores <- parallel_core_cnt(args)
    stopifnot(OK)
    args
}


#' Given a bunch of filenames, find all four-number digit strings.
#' @param filenames A list of character strings to search.
#' @return A list of numeric years.
years_in_filenames <- function(filenames) {
    matches <- regexpr("[0-9]{4}", filenames)
    start <- as.numeric(matches)
    stop <- as.numeric(attr(matches, "match.length")) + start - 1
    as.numeric(substr(filenames, start, stop))
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
    pfpr_years <- years_in_filenames(pfpr_files)
    am_dir <- rampdata::workflow_path("am")
    am_files <- list.files(rampdata::as.path(am_dir), "*.tif")
    am_years <- years_in_filenames(am_files)
    # We could choose, instead, to extend either one to cover missing years
    # in the other.
    shared_years <- sort(intersect(am_years, pfpr_years))
    if (!is.null(select_years)) {
        shared_years <- sort(intersect(shared_years, select_years))
    }

    pr_fn <- paste0("PfPR_median_Africa_admin0_", shared_years[1], ".tif")
    pfpr_file <- rampdata::add_path(pfpr_dir, file = pr_fn)
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
data_for_country <- function(country_alpha3, years) {
    pfpr_dir <- rampdata::workflow_path("pfpr")
    am_dir <- rampdata::workflow_path("am")

    outline_sf <- gadm_country_shapefile(country_alpha3)
    bbox <- sf::st_bbox(outline_sf)
    pr_fn <- paste0("PfPR_median_Africa_admin0_", years[1], ".tif")
    pfpr_file <- rampdata::add_path(pfpr_dir, file = pr_fn)
    pfpr <- raster::raster(rampdata::as.path(pfpr_file))
    domain_extent <- pixel_bounding_box(pfpr, bbox)
    flog.debug("domain_extent", paste(domain_extent, collapse = ","))

    pfpr_all <- list()
    am_all <- list()
    for (year_idx in 1:length(years)) {
        year <- years[year_idx]
        pr_fn <- paste0("PfPR_median_Africa_admin0_", year, ".tif")
        pfpr_file <- rampdata::add_path(pfpr_dir, file = pr_fn)
        if (!file.exists(rampdata::as.path(pfpr_file))) {
            msg <- paste("Cannot find pfpr at", rampdata::as.path(pfpr_file))
            flog.error(msg)
            stop(msg)
        }

        am_file <- rampdata::add_path(
            am_dir,
            file = paste0(year, ".effective.treatment.tif"))
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
    list(
        pfpr = combine(pfpr_all),
        am = combine(am_all)
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
    pr_to_ar_dt <- data.table::fread(rampdata::as.path(pr2ar_rp))
    pfpr_dir <- rampdata::workflow_path("pfpr")
    am_dir <- rampdata::workflow_path("am")

    pfpr_all <- list()
    am_all <- list()
    for (year_idx in 1:length(years)) {
        year <- years[year_idx]
        pr_fn <- paste0("PfPR_median_Africa_admin0_", year, ".tif")
        pfpr_file <- rampdata::add_path(pfpr_dir, file = pr_fn)
        if (!file.exists(rampdata::as.path(pfpr_file))) {
            msg <- paste("Cannot find pfpr at", rampdata::as.path(pfpr_file))
            flog.error(msg)
            stop(msg)
        }

        am_file <- rampdata::add_path(
            am_dir,
            file = paste0(year, ".effective.treatment.tif"))
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


#' Given a list of matrices, construct a single three-dimensional array.
#' @param list_of_matrices Each item in the list is a two-dimensional array of the same size.
#' @return A single three-dimensional array where the first dimension is the
#'     number of items in the list and the next two are the dimensions of each array.
collect_and_permute <- function(list_of_matrices) {
    sentinel <- dim(list_of_matrices[[1]])
    together <- array(do.call(c, list_of_matrices), dim = c(sentinel, length(list_of_matrices)))
    aperm(together, c(3, 1, 2))
}


#' Given the raw input data, prepare it for computation.
#' @param data A list where each item is input data.
#' @param work The matrix of which tiles to do.
#' @param process_extent Bounding box of loaded data within the computational domain.
#' @param blocksize The size of each tile.
#' @return A list with transformed data items. This list is broken
#'     up into the work tiles, so that we don't transfer the whole
#'     image to each processor when we do calculations.
#'
#' This is data movement to set up for the next algorithm.
prepare_timeseries <- function(data, work, process_extent, blocksize) {
    stopifnot(is_tile(blocksize))
    pfpr <- collect_and_permute(data$pfpr)
    am <- collect_and_permute(data$am)

    work_cnt <- dim(work)[2]
    pieces <- lapply(1:work_cnt, FUN = function(col) {
        tile <- work[, col]
        tile_domain_extent <- tile_extent(tile, blocksize)
        # relative extent to the subset of data that is loaded.
        rel <- find_relative_extent(process_extent, tile_domain_extent)
        # flog.debug(paste("tile", paste(tile, collapse = ","),
        #     "relative extent", paste(rel, collapse = ","), "\n"
        # ))
        pfpr_chunk <- pfpr[, rel["rmin"]:rel["rmax"], rel["cmin"]:rel["cmax"]]
        am_chunk <- am[, rel["rmin"]:rel["rmax"], rel["cmin"]:rel["cmax"]]
        list(tile = tile, pfpr = pfpr_chunk, am = am_chunk)
    })

    list(
        parameters = data$parameters,
        pr_to_ar_dt = data$pr_to_ar_dt,
        years = data$years,
        chunks = pieces
        )
}


#' Find nonzeros in each chunk of tile data to see they loaded.
chunk_nonzero <- function(chunks) {
    vapply(
        chunks,
        function(chunk) {
            c(sum(!is.na(chunk$pfpr)), sum(!is.na(chunk$am)))
        },
        integer(2)
    )
}


#' Find the sum of each chunk of tile data to see they differ.
chunk_sums <- function(chunks) {
    vapply(
        chunks,
        function(chunk) {
            c(sum(chunk$pfpr, na.rm = TRUE), sum(chunk$am, na.rm = TRUE))
        },
        numeric(2)
    )
}

kappa_rm <- function(pfpr, c) { c * pfpr }
alpha_from_eir <- function(eir, pfpr, rho) {
  alpha <- 1 - exp(-log(1 + eir * k * tau) / k)
  alpha + ar2pr(pfpr, rho) - ar2pr(pfpr, 0)
}
rc_basic <- function(b, c, r, V, am) { b * V * c / r }

#' Compute the results of analyzing a time series of PfPR and AM.
#' @param pfpr A numeric vector of PR. Length is number of time points.
#' @param am A numeric vector of case management measure. Length is number of time points.
#' @param pr2ar A function that gives AR from PR and rho.
#' @param params Other parameters for the calculation in a list.
#' @return A list where each variable has results for the time series
#'     as a vector of numeric. The entries are `alpha`, `kappa`, `eir`, `vc`, `rc`.
#'
#' The parameters passed into this function are defined in the config
#' file, whose default name is `rc_kappa.toml`, in the `[parameters]` section.
pixel_work <- function(pfpr, am, params, strategies) {
  # PfPR and AM come in with the same set of NA patterns, where there is no land.
  with(c(params, strategies), {
    alpha <- pr_to_ar(pfpr, kam * am)
    kappa <- kappaf(pfpr, c)
    h <- -log(1 - alpha) / tau
    eir <- (exp(h * k * tau) - 1) / (b * k * tau)
    V <- eir / kappa
    R <- rcf(b, c, r, V, am)
    list(
      alpha = alpha,
      kappa = kappa,
      eir = eir * 365,
      vc = V,
      rc = R
    )
  })
}

pr2eir=function(x, b=0.55, r=1/200, k=4.2){
  ((1-x)^-k-1)*(r/k/b)
}


#' Check that the results are all finite.
#' @param results A list of arrays.
#' @param inputs The input list of arrays, used for debugging failure.
#' @return The input results. This will stop on failure.
#' We shouldn't have any NANs here. Good data in, good data out.
#' This is what let's us ignore NA later. This function is expensive
#' because it scans all points, but it's important.
check_linearized_outputs <- function(results, inputs) {
    for (check_idx in seq_along(results)) {
        check_arr <- results[[check_idx]]
        if (!all(is.finite(check_arr))) {
            bad_ones <- !is.finite(check_arr)
            res_bad <- data.frame(c(inputs, results))[bad_ones,]
            debug_file <- "debug.csv"
            data.table::fwrite(res_bad, file = debug_file)
            inf_cnt <- sum(bad_ones)
            nan_cnt <- sum(is.na(check_arr))
            var_name <- names(results)[check_idx]
            flog.warn(sprintf(
                "Found %d nan %d inf in %s. Writing to %s\n",
                nan_cnt, inf_cnt, var_name,
                normalizePath(debug_file)
                ))
            stopifnot(nan_cnt == 0 && inf_cnt == 0)
        }
    }
    results
}


#' Takes an input array in multiple dimensions and linearizes and removes NA.
#' This replaces the outputs back into multiple dimensions after making the call.
linearized_work <- function(input_list, run_func) {
    # There are NA values in the arrays in the input_list.
    # The arrays are probably three-dimensional.
    not_available <- lapply(input_list, function(check) is.na(check))
    input_dims <- dim(input_list[[1]])
    dim_keep <- length(input_dims)
    # If any array has an NA in a position in the array, then that position is
    # thrown out for all input arrays.
    keep <- !rowSums(array(do.call(c, not_available), dim = c(input_dims, 2)), dims = dim_keep)
    linear <- lapply(input_list, function(data) {
        as.numeric(data[keep])
    })

    # How could we have a tile with no data, when they were excluded?
    # They were excluded based on PR, not AM.
    if (sum(keep) > 0) {
        results <- check_linearized_outputs(run_func(linear), linear)
    } else {
        # Don't call the function with no data, but how do we know
	# variable names? Call the function with small, fake data.
        just_names <- run_func(list(pfpr=c(0.1), am = c(0.02)))
	results <- lapply(just_names, function(x) numeric(0))
    }

    # The returned values will have NA for the same pattern in all outputs.
    lapply(results, function(result) {
        complete <- array(dim = input_dims)
        complete[keep] <- result
        array(complete, dim = input_dims)
    })
}


#' Accumulate the individual time series from a computation that doesn't care about time.
#' @param chunk A list with pfpr, am, and the tile.
#' @param parameters A list of parameters for the calculation.
#' @param pr_to_ar_dt The pr-ar data from the mechanistic model.
#' @return a list with an entry for each of alpha, kappa, eir, vc, rc.
#'     It also identifies the source block by block_id.
#'
#' If the pixel function doesn't care about time, then there is no
#' reason to loop over the matrix. Turn the matrix into a single vector
#' and then do the computation.
over_block <- function(chunk, parameters, pr_to_ar_dt) {
    parameters$pr_to_ar <- ar_of_pr_rho(pr_to_ar_dt)
    parameters$ar2pr <- build_ar2pr(pr_to_ar_dt)
    strategies <- list(kappaf = kappa_rm, rcf = rc_basic)
    # The first dimension is the timeseries, so row=2, col=3.
    run_func <- function(plaquette) {
        pixel_work(plaquette$pfpr, plaquette$am, parameters, strategies)
    }
    only_data <- chunk[names(chunk)[!names(chunk) %in% "tile"]]
    results <- linearized_work(only_data, run_func)
    flog.debug(paste("chunk", paste(chunk$tile, collapse = ","),
        "has", sum(is.na(results$alpha)), "na alpha values and",
        sum(results$alpha > 0 & results$alpha < 1, na.rm = TRUE),
        "in 0 < x < 1."
    ))
    results[["array_names"]] <- names(results)
    results[["block"]] <- chunk$tile
    results
}


#' Takes a parameter that is an expression (of drawing a distribution) and calls it.
draw_parameters <- function(parameters, N) {
    if (N > 1) {
        draw_params <- with(parameters, {
            data.frame(
                kam = rep(kam, N),
                r = rnorm(N, r, r_sd),
                k = k,
                ku = runif(N),  # Use as quantile within calculation.
                b = rbeta(N, b, b_shape1, b_shape2),
                c = rep(c, N),
                tau = rep(tau, N),
                D_low = rep(D_low, N),
                D_high = rep(D_high, N)
            )
        })
    } else {
        draw_params <- data.frame(
            kam = kam,
            r = r,
            k = k,  # mean for k calculation.
            ku = 0.5,  # Use as quantile within calculation.
            b = b,
            c = c,
            tau = tau,
            D_low = D_low,
            D_high = D_high
        )
    }
    # If any parameters are less than zero, we should redraw.
    # It's a rejection method to get truncated distributions.
    stopifnot(all(as.matrix(draw_params) > 0))
    draw_params
}


pixel_two <- function(pfpr, am, params, strategies) {
  # PfPR and AM come in with the same set of NA patterns, where there is no land.
  with(c(params, strategies), {
    rho <- kam * am
    # draw of q determined by previous uniform draw.
    kd <- qgamma(params$ku, k * (1 - pfpr)^.6 * 5, 5 * (1 - pfpr)^.6)
    stopifnot(length(kd) == length(pfpr))
    # We aren't using the mechanistic model to get an absolute value of alpha
    # because that alpha isn't close enough. We are using that model to
    # estimate how much treatment would shift PfPR, given an alpha.
    alpha_cronus <- pr_to_ar(pfpr, rho)
    # nt = no_treatment
    pfpr_nt <- ar2pr(alpha_cronus)
    # This EIR estimate is from a fit to data.
    max_deir <- 1500 / 365
    eir_nt <- pmin(
      pr2eir(pfpr_nt, b, r, kd),
      rep(Inf, length(pfpr_nt))
      )
    alpha_nt <- 1 - (eir_nt * b * kd * tau + 1)^(-1/kd)
    # Then go back and estimate how much treatment would mean alpha was higher.
    lower_ar <- pr_to_ar(pfpr, numeric(length(pfpr)))
    upper_ar <- pr_to_ar(pfpr, rho)
    # Move it the same relative distance to the upper bound.
    alpha <- 1 - (1 - alpha_nt) * (1 - upper_ar) / (1 - lower_ar)
    # alpha <- min(alpha, 1 - 1e-10)

    kappa <- c * pfpr_nt
    h <- -log(1 - alpha) / tau
    # alpha can be 1 after the shift, so cap FOI.
    # These will be infinite, and that's OK. We don't
    # get rid of NAN. Those are failures.
    max_finite <- 1e10
    h[h > max_finite] <- max_finite
    eir <- (exp(h * kd * tau) - 1) / (b * kd * tau)
    eir[eir > max_finite] <- max_finite
    V <- ifelse(
        kappa > 0,
        eir / kappa,
        0
    )
    V[V > max_finite] <- max_finite
    R <- b * V * c * ((1 - rho) * D_high + rho * D_low)
    R[R > max_finite] <- max_finite
    list(
      alpha = alpha,
      kappa = kappa,
      foi = h,
      # deir = daily eir, aeir = annual eir = 365 * deir.
      aeir = eir * 365,
      vc = V,
      rc = R
    )
  })
}


#' The draws are hierarchical by variable and then quantile. This flattens them.
#' @param draw_list The hierarchical list with names.
#' @return A single list with the two names concatenated.
flatten_draws <- function(draw_list) {
    ll <- draw_list
    sl_cnt <- length(ll[[1]])
    draw_names <- lapply(1:(sl_cnt * length(ll)), function(idx) {
        list_idx <- 1 + (idx - 1) %/% sl_cnt
        sub_list <- ll[[list_idx]]
        sublist_idx <- 1 + (idx - 1) %% sl_cnt
        name <- names(ll)[list_idx]
        sub_name <- names(sub_list)[sublist_idx]
        paste(name, sub_name, sep = "_")
    })
    draw_vals <- lapply(1:(sl_cnt * length(ll)), function(idx) {
        list_idx <- 1 + (idx - 1) %/% sl_cnt
        sub_list <- ll[[list_idx]]
        sublist_idx <- 1 + (idx - 1) %% sl_cnt
        sub_list[[sublist_idx]]
    })
    names(draw_vals) <- draw_names
    draw_vals
}


#' Given a list of draws, summarize mean and confidence intervals.
#' @param draws A list of draws. Each draw is a list of variables,
#'     where each variable is a three-dimensional array.
#' @param confidence_percent As in, 95.
#' @return a list with one name for each variable. Each list
#'     item has three parts, lower, median, and upper, for that variable.
#'     Each of lower, median, and upper has dimensions year, row, col.
summarize_draws <- function(draws, confidence_percent) {
    ci <- 0.5 * (1 - 0.01 * confidence_percent)
    quantiles <- c(ci, 0.5, 1 - ci)
    draw_cnt <- length(draws)
    array_names <- names(draws[[1]])
    array_dim <- dim(draws[[1]][[1]])
    summarized <- lapply(array_names, function(name) {
        var_data <- array(
            do.call(
                c,
                lapply(1:draw_cnt, function(combine_idx) {
                    draws[[combine_idx]][[name]]
                })
                ),
            dim = c(array_dim, draw_cnt)
        )
        draws_first <- aperm(var_data, c(4, 1, 2, 3))
        quantile_first <-apply(
            draws_first,
            MARGIN = c(2, 3, 4),
            # We can ignore NA here b/c we checked earlier, and some pixels
            # _should be all NA_ in water.
            FUN = function(x) quantile(x, quantiles, na.rm = TRUE)
            )
        quantile_last <- aperm(quantile_first, c(2, 3, 4, 1))
        list(
            lower = quantile_last[, , , 1],
            median = quantile_last[, , , 2],
            upper = quantile_last[, , , 3]
        )
    })
    names(summarized) <- array_names
    flattened <- flatten_draws(summarized)
    lapply(flattened, function(x) {
      dimnames(x) <- list(year=NULL, row=NULL, col=NULL)
      x
    })
}


#' Run the pixel function on a chunk of data.
#' @param chunk a list with the pfpr and am datasets.
#' @param parameters a data.frame, at this point, with a row per draw.
#' @param pr_to_ar_dt The pr2ar mesh data.
#' @param confidence_percent Confidence interval, out of 100.
over_block_draw <- function(chunk, parameters, pr_to_ar_dt, confidence_percent) {
    stopifnot(is.numeric(confidence_percent))
    strategies <- list(
        pr_to_ar = ar_of_pr_rho(pr_to_ar_dt),
        ar2pr = build_ar2pr(pr_to_ar_dt)
    )
    # The function we call assumes that every entry in the chunk
    # is a dataset, so remove the tile.
    only_data <- chunk[names(chunk)[!names(chunk) %in% "tile"]]
    # The first dimension is the timeseries, so row=2, col=3.
    draw_cnt <- nrow(parameters)
    draws <- lapply(1:draw_cnt, function(draw_idx) {
        run_func <- function(plaquette) {
            pixel_two(plaquette$pfpr, plaquette$am, parameters[draw_idx, ], strategies)
        }
        linearized_work(only_data, run_func)
    })
    summarized <- summarize_draws(draws, confidence_percent)
    flog.debug(paste("summarized names", paste0(names(summarized), collapse = ",")))
    summarized[["block"]] <- chunk$tile
    summarized
}


#' Given everything from lapply, combine it.
#' @param tile_output is a list where each entry is a all data for a single
#'     tile. Each entry has one member named `block`, which identifies the
#'     tile. The ther entries are the arrays in (year, row, col).
#' @param blocksize An named array with "row" and "col".
#' @param only_name The name of a single array to combine from all tiles.
#' @return The data combined, where each piece is offset by the same offset.
combine_output <- function(tile_output, blocksize, only_name = NULL) {
    work <- vapply(
        tile_output,
        function(one_result) {
            one_result$block
        },
        FUN.VALUE = numeric(2)
    )
    loaded_extent <- raster_extent_from_work(work, blocksize)
    row_cnt <- unname(loaded_extent["rmax"] - loaded_extent["rmin"] + 1)
    col_cnt <- unname(loaded_extent["cmax"] - loaded_extent["cmin"] + 1)
    offset <- c(row = unname(loaded_extent["rmin"]), col = unname(loaded_extent["cmin"]))

    array_names <- names(tile_output[[1]])
    array_names <- array_names[!array_names %in% c("block", "array_names")]
    if (!is.null(only_name)) {
      if (only_name %in% array_names) {
        array_names <- only_name
      } else {
        msg <- sprintf(
          "Cannot find %s in array names %s in combine_output",
          only_name, paste0(array_names, collapse = ","))
        flog.error(msg)
        stop()
      }
    }  # else use all names, as planned.
    example <- tile_output[[1]][[array_names[1]]]
    # The order of axes is year, row, col.
    time_cnt <- dim(example)[1]
    if (any(dim(example)[2:3] != blocksize)) {
        flog.error(paste("expect blocksize to match array size",
            paste(dim(example), collapse = ",")))
    }

    # allocate first so we aren't churning memory.
    united <- lapply(array_names, function(var) {
      array(dim = c(time_cnt, row_cnt, col_cnt))
      })
    names(united) <- array_names

    for (work_idx in 1:length(tile_output)) {
        tile <- tile_output[[work_idx]]$block
        abs_extent <- tile_extent(tile, blocksize)
        stopifnot(sum(is.na(abs_extent)) == 0L)
        rel_extent <- c(
            abs_extent["rmin"] - unname(offset["row"]) + 1L,
            abs_extent["rmax"] - unname(offset["row"]) + 1L,
            abs_extent["cmin"] - unname(offset["col"]) + 1L,
            abs_extent["cmax"] - unname(offset["col"]) + 1L
        )
        re <- rel_extent
        for (name in array_names) {
            united[[name]][, re["rmin"]:re["rmax"], re["cmin"]:re["cmax"]] <- tile_output[[work_idx]][[name]]
        }
    }
    flog.debug(paste("combined has", sum(is.na(united$alpha)),
        "na alpha values and",
        sum(united$alpha > 0 & united$alpha < 1, na.rm = TRUE),
        "in 0 < x < 1."
    ))
    united
}


build_outvars_dir <- function(outarg) {
    dest_dir <- rampdata::workflow_path("outvars")
    if (!is.null(outarg)) {
        dest_dir <- rampdata::add_path(dest_dir, version = outarg)
    }
    dest_fn <- rampdata::as.path(dest_dir)
    flog.debug(paste("outvars dest_fn is", dest_fn))
    if (!dir.exists(dest_fn)) {
        dir.create(dest_fn, recursive = TRUE)
    }
    dest_dir
}



.plot.kinds <- list(
  rc = list(
    title = "Rc",
    breaks = c(0.0, 0.5, 1, 2, 5, 10, 100, Inf)
  ),
  vc = list(
    title = "Vectorial Capacity",
    breaks = c(0.0, 0.5, 1, 2, 5, 10, 100, Inf)
  ),
  aeir = list(
    title = "Annual EIR",
    breaks = c(0.0, 0.1, 0.5, 1, 5, 10, 100, Inf)
  )
)


#' Customized plotting for different variables.
#' @param raster_obj A `raster::raster` object.
#' @param filename A string filename for the png output.
#' @param name The variable and year, with an underscore between. Has to be this.
#' @param year The numeric year.
#' @param options A list of options for png width, height, and resolution.
#' @param admin0 A shapefile with admin0 boundaries.
plot_as_png <- function(raster_obj, filename, name, year, options, admin0) {
  parts <- strsplit(name, "_")[[1]]
  kind <- parts[1]
  quantile <- parts[2]

  if (kind %in% names(.plot.kinds)) {
    aplot <- tmap::tm_shape(vc_median) +
      tmap::tm_raster(
        breaks = .plot.kinds[[kind]]$breaks,
          title = sprintf(
	    "%s %s %d",
            .plot.kinds[[kind]]$title,
            tools::toTitleCase(quantile),
            year
        )
      ) +
       tmap::tm_layout(
        legend.position = c("left", "bottom") 
      ) +
      tmap::tm_shape(admin0lakes) +
      tmap::tm_borders()
    tmap::tmap_save(aplot, filename, asp = 0, height = options$pngheight)
  } else {
     png(
       file = filename,
       width = options$pngwidth,
       height = options$pngheight,
       res = options$pngres
     )
     raster::plot(
       log10(raster_obj),
       main = sprintf("log10(%s) %d", name, year),
       useRaster = TRUE
     )
     dev.off()
  }
}


#' Given outputs for a tile, write them to a folder.
#' @param output A list of results.
#' @param years Years that will save.
#' @param domain_extent Bounds and resolution of the raster.
#' @param args Command-line arguments that have been parsed.
#' @param options Settings from the toml file.
write_output <- function(output, years, domain_extent, args, options) {
    flog.debug(paste("writing output", sum(is.na(output[[1]])),
        "na values in first array and",
        sum(output[[1]] > 0 & output[[1]] < 1, na.rm = TRUE),
        "in 0 < x < 1."
    ))

    outline_rp <- rampdata::ramp_path("/inputs/country_outlines/201122")
    admin0 <- sf::st_read(rampdata::as.path(rampdata::add_path(
        outline_rp, file = "ne_10m_admin_0_countries_lakes")))

    dest_dir <- build_outvars_dir(args$outvars)
    for (name in names(output)) {
        if (!name %in% names(output)) {
            msg <- paste("output doesn't have the", name, "array")
            flog.error(msg)
            stop(msg)
        }
        by_year <- aperm(output[[name]], c(2, 3, 1))
        dm <- dim(by_year)[1:2]
        for (year_idx in 1:length(years)) {
            year <- years[year_idx]
            out_rp <- rampdata::add_path(dest_dir, file = paste0(name, "_", year, ".tif"))
            # XXX I'm worried that this should be transposed.
            ready_data <- by_year[, , year_idx]
            out_fn <- rampdata::as.path(out_rp)
            exists <- file.exists(out_fn)
            if (exists & args$overwrite) {
                unlink(out_fn)
                exists <- FALSE
            }
            if (!exists) {
                flog.info(paste("writing file", out_fn))
                raster_obj <- raster::raster(
                    nrows = dm[1], ncols = dm[2],
                    xmn = domain_extent$xmin,
                    xmx = domain_extent$xmax,  # columns are x for raster.
                    ymn = domain_extent$ymin,
                    ymx = domain_extent$ymax,
                    crs = domain_extent$projection
                )
                raster_obj <- raster::setValues(raster_obj, ready_data)
                raster::writeRaster(raster_obj, filename = out_fn, format = "GTiff")

                png_rp <- rampdata::add_path(out_rp, file = sprintf("%s_%d.png", name, year))
		png_fn <- rampdata::as.path(png_rp)
                flog.info(paste("writing file", png_fn))
                plot_as_png(raster_obj, png_fn, name, year, options, admin0)
            } else {
                flog.error(paste("cannot overwrite", out_fn))
            }
        }
    }
}


#' Given command-line arguments, determine how much data to load.
#' It's deciding splits without consulting arguments for now.
#' This decomposes the domain into tiles because, when we parallelize
#' with the parallel package, each task needs at least a few seconds
#' of work to do in order to make data transfer worthwhile.
#' @param available A description of what data could be processed.
#' @param options These parameters define `blocksize` and `single_tile_max`.
#' @return The same list with an added `tiles` member that describes
#'     the dimensions of computational tiles relative to the origin
#'     of the `domain_extent`, the portion we will compute.
plan_domain_decomposition <- function(available, options) {
    extent <- available$domain_extent
    stopifnot(is_extent(extent))
    rows <- (extent["rmax"] - extent["rmin"] + 1L)
    cols <- (extent["cmax"] - extent["cmin"] + 1L)
    if (rows * cols < options$single_tile_max) {
        tile_cnt <- c(row = 1L, col = 1L)
        blocksize <- c(row = unname(rows), col = unname(cols))
    } else {
        bs <- as.integer(options$blocksize)
        blocksize <- c(row = bs, col = bs)
        tile_cnt <- c(
            row = as.integer(ceiling(rows / blocksize["row"])),
            col = as.integer(ceiling(cols / blocksize["col"]))
            )
    }
    available[["tiles"]] <- list(
        tile_cnt = tile_cnt,
        blocksize = blocksize
    )
    available
}


#' Removes tiles that have all NA for this example slice.
#' Tile is within the domain. The sample data is already cut
#' to that domain.
remove_tiles_over_water <- function(available) {
    stopifnot("sample_slice" %in% names(available))
    # This slice is cut from the domain.
    pfpr <- available$sample_slice
    tile_cnt <- available$tiles$tile_cnt
    blocksize <- available$tiles$blocksize

    ridx <- rep(1:tile_cnt["row"], tile_cnt["col"])
    cidx <- rep(1:tile_cnt["col"], each = tile_cnt["row"])
    flog.debug("examining %d tiles for NA", length(ridx))
    tiles <- rbind(ridx, cidx)
    # If there is only one tile, we need to insist this is 2D array.
    dim(tiles) <- c(2, length(ridx))
    dimnames(tiles) <- list(c("row", "col"), NULL)
    keep_tile <- vapply(1:(dim(tiles)[2]), function(tile_idx) {
        extent <- tile_extent(tiles[, tile_idx], blocksize)
        rmax <- min(extent["rmax"], dim(pfpr)[1])
        cmax <- min(extent["cmax"], dim(pfpr)[2])
        chunk <- pfpr[
            extent["rmin"]:rmax,
            extent["cmin"]:cmax
        ]
        sum(!is.na(chunk)) > 0
    }, FUN.VALUE = logical(1))
    removed_cnt <- sum(!keep_tile)
    kept_cnt <- sum(keep_tile)
    flog.info(sprintf("Removing %d tiles with no data keep %d",
        removed_cnt, kept_cnt))
    remaining <- tiles[, keep_tile]
    # Again, if there is only one to keep, insist it be 2D.
    dim(remaining) <- c(2, kept_cnt)
    dimnames(remaining) <- list(c("row", "col"), NULL)
    remaining
}


#' Takes the extent of the tiles and returns a list of which tiles to do.
#' @param tiles A list of tiles as a 2xN matrix.
#' @return An array that's 2 x the number of tiles. It's row, col.
#' Assumes that all work is done by the same process.
#' This orders the work so that it uses cache better.
task_work_from_tiles <- function(tiles) {
    # The Hilbert curve is a space-filling curve, so it orders the 2D tiles.
    tidx <- vapply(1:dim(tiles)[2], FUN = function(tile_idx) {
            t <- tiles[, tile_idx]
            encode_hilbert(t["row"], t["col"])
        }, FUN.VALUE = integer(1)
        )
    by_curve <- order(tidx)
    tiles[, by_curve]
    dim(tiles) <- c(2, length(by_curve))
    dimnames(tiles) <- list(c("row", "col"), NULL)
    tiles
}


#' Takes the extent of the tiles and returns a list of which tiles to do.
#' @param tile_cnt The tile_cnt from the plan from `plan_domain_decomposition`.
#' @return An array that's 2 x the number of tiles. It's row, col.
#' Assumes that all work is done by the same process.
#' This orders the work so that it uses cache better.
task_work <- function(tile_cnt) {
    stopifnot(is_tile(tile_cnt))
    ridx <- rep(1:tile_cnt["row"], tile_cnt["col"])
    cidx <- rep(1:tile_cnt["col"], each = tile_cnt["row"])
    # The Hilbert curve is a space-filling curve, so it orders the 2D tiles.
    tidx <- vapply(1:length(cidx), FUN = function(tile_idx) {
            encode_hilbert(ridx[tile_idx], cidx[tile_idx])
        }, FUN.VALUE = integer(1)
        )
    by_curve <- order(tidx)
    work <- rbind(ridx[by_curve], cidx[by_curve])
    dimnames(work) <- list(c("row", "col"), NULL)
    work
}


tile_extent <- function(tile, blocksize) {
    stopifnot(is_tile(tile))
    stopifnot(is_tile(blocksize))
    c(
        rmin = unname((tile["row"] - 1L) * blocksize["row"] + 1L),
        rmax = unname(tile["row"] * blocksize["row"]),
        cmin = unname((tile["col"] - 1L) * blocksize["col"] + 1L),
        cmax = unname(tile["col"] * blocksize["col"])
    )
}


#' Given a list of tiles, figures out the bounding box of pixel values.
#' @param work A 2d matrix with a dimname on the first dimension of "row", "col".
#' @param blocksize A row and col size for blocks for the work.
raster_extent_from_work <- function(work, blocksize) {
    stopifnot(all(dimnames(work)[[1]] == c("row", "col")))
    stopifnot(is_tile(blocksize))
    small <- tile_extent(c(row = min(work["row",]), col = min(work["col", ])), blocksize)
    large <- tile_extent(c(row = max(work["row",]), col = max(work["col", ])), blocksize)
    c(small["rmin"], large["rmax"], small["cmin"], large["cmax"])
}


parallel_core_cnt <- function(args = NULL) {
    # So the user can specify in command-line arguments.
    if (is.list(args) && "cores" %in% names(args) && !is.null(args$cores)) {
        cores <- suppressWarnings(as.integer(args$cores))
        if (!is.na(cores) && length(cores) == 1L && cores > 0L) {
            return(cores)
        } else {
            flog.warn(paste("Could not use args$cores to set core count:", args$cores))
        }
    }
    # Because we run on the cluster, and the cluster sets environment variables.
    for (env_var in c("NCPUS", "fthread")) {
        env_str <- Sys.getenv(env_var, unset = NA)
        if (!is.na(env_str)) {
            env_cores <- suppressWarnings(as.integer(env_str))
            if (length(env_cores) == 1L && env_cores > 0L) {
                return(env_cores)
            }
        }
    }
    # Because we want all the cores otherwise.
    core_cnt <- parallel::detectCores()
    if (is.integer(core_cnt)) {
        return(core_cnt)
    } else {
        return(2L)
    }
}


#' If you want to run from the R command line, you can call this.
#' Use the same arguments as the commandline, but sent in as a list.
#' For example,
#' \Dontrun{
#' args <- check_args(arg_parser(c(
#'    "--config=gen_scaled_ar/rc_kappa.toml",
#'    "--country=gmb",
#'    "--years=2010:2011",
#'    "--overwrite"
#'    )))
#' args$country <- "uga"
#' funcmain(args)
#' }
funcmain <- function(args) {
    # We will want to split this work different ways for development,
    # laptops, and the cluster, so we use an explicit domain decomposition.
    rampdata::initialize_workflow(args$config)
    options <- configr::read.config(args$config)[["options"]]

    # What can be done and which part we choose to do.
    available <- available_data(args$inversion, args$country, args$years)
    # How to split that up.
    plan <- plan_domain_decomposition(available, options)
    remaining_tiles <- remove_tiles_over_water(plan)

    # The part that this particular process will do.
    work <- task_work_from_tiles(remaining_tiles)
    flog.debug(paste("tile count", dim(work)[2], "size",
        paste(plan$tiles$blocksize, collapse = ",")))

    # By this point, someone told us what to do, so load all the data.
    process_extent <- raster_extent_from_work(work, plan$tiles$blocksize)
    flog.debug(paste(
        "extent", paste(process_extent, collapse = ","), "names",
        paste(names(process_extent), collapse = ","), "\n"))
    # The load_extent is what the process needs, relative to the domain coordinates.
    load_extent <- apply_relative_extent(plan$domain_extent, process_extent)
    data <- load_data(args$config, args$pr2ar, load_extent, args$years)
    # We permute data axes for more efficient running and
    # split the data into chunks for processing.
    timeseries_data <- prepare_timeseries(data, work, process_extent, plan$tiles$blocksize)
    data <- NULL  # Help the memory cleanup. This data is now in timeseries chunks.

    params <- timeseries_data$parameters
    # Uses the parallel random generation from the `parallel` package.
    set.seed(params$random_seed, "L'Ecuyer")
    draw_params <- draw_parameters(params, args$draws)
    chunks <- timeseries_data$chunks

    core_cnt <- parallel_core_cnt(args)
    flog.debug("%d chunks to %d cores", length(chunks), core_cnt)
    # Work in parallel on this machine.
    # Thread output usually goes to /dev/null, but we can redirect
    # it to a file. Put it with our data for now. Could be a mess
    # for thread contention, but they don't write a lot.
    thread_out_dir <- build_outvars_dir(args$outvars)
    thread_out <- rampdata::add_path(thread_out_dir, file = "threads.txt")
    thread_fn <- rampdata::as.path(thread_out)
    flog.debug(sprintf("thread file out to %s", thread_fn))
    if (core_cnt > 1) {
    cluster <- parallel::makeCluster(
    	    core_cnt,
	    type = "FORK",
	    outfile = thread_fn
	    )
    output <- parallel::clusterApplyLB(
        cluster,
        chunks,
        function(chunk) {
            flog.debug(".")
            over_block_draw(
                chunk,
                draw_params,
                timeseries_data$pr_to_ar_dt,
                params$confidence_percent
                )
        }
    )
    flog.debug("End of parallel work.")
    stopCluster(cluster)
    } else {
    output <- lapply(
        chunks,
        function(chunk) {
            flog.debug(".")
            over_block_draw(
                chunk,
                draw_params,
                timeseries_data$pr_to_ar_dt,
                params$confidence_percent
                )
        }
    )
    }

    # The output chunks need to be reassembled before writing.
    ready_to_write <- combine_output(output, plan$tiles$blocksize)
    write_output(ready_to_write, args$years, plan$domain_dimensions, args, options)
}


#' If you want to run from the Bash command line, you can call this.
main <- function() {
    flog.threshold(DEBUG)
    improved_errors()
    args <- check_args(arg_parser())
    funcmain(args)
}


construct_plan <- function(args) {
  # We will want to split this work different ways for development,
  # laptops, and the cluster, so we use an explicit domain decomposition.
  rampdata::initialize_workflow(args$config)
  options <- configr::read.config(args$config)[["options"]]

  # What can be done and which part we choose to do.
  available <- available_data(args$inversion, args$country, args$years)
  # How to split that up.
  plan <- plan_domain_decomposition(available, options)
  remaining_tiles <- remove_tiles_over_water(plan)

  # The part that this particular process will do.
  work <- task_work_from_tiles(remaining_tiles)
  flog.debug(paste("tile count", dim(work)[2], "size",
                   paste(plan$tiles$blocksize, collapse = ",")))
  out_dir <- build_outvars_dir(args$outvars)
  data.table::fwrite(t(work), file = rampdata::as.path(rampdata::add_path(out_dir, file = "work.csv")))
  save_plan(plan, rampdata::as.path(rampdata::add_path(out_dir, file = "plan.json")))
}


#' Pull my work, according to my task index.
#' @param task_id integer index of task within all tasks. One-based.
#' @param task_cnt integer total number of tasks.
#' @param work Data frame with row and column of tile.
#' @return 2D array that is (row, col) by work items.
#'
#' The math here is so that, if the tasks don't split evenly, there
#' isn't a single task with only one thing to do. The first few
#' tasks get one more, and the rest get one less.
my_tasks <- function(task_id, task_cnt, work) {
  n <- nrow(work)
  m <- task_cnt
  w <- 1 + (n - 1) %/% m
  x <- n - m * (w - 1)
  if (task_id <= x) {
    rows <- (1 + (task_id - 1) * w):(task_id * w)
  } else {
    if (w > 1) {
      r <- task_id - x
      b <- x * w
      rows <- (b + 1 + (r - 1) * (w - 1)):(b + r * (w - 1))
    } else {
      return(NULL)
    }
  }
  as.matrix(t(work[rows, ]))
}


save_outputs <- function(chunks_fn, outputs) {
  for (tile_idx in seq_along(outputs)) {
    output <- outputs[[tile_idx]]
    # The outputs is a list of arrays, by name, and a single
    # "block" member to say what the tile is.
    block <- output$block
    stopifnot(!is.null(block))
    chunks <- list()
    for (chunk_idx in seq_along(output)) {
      cn <- names(output)[chunk_idx]
      if (cn != "block") {
        arr <- output[[chunk_idx]]
        attr(arr, "tile") <- block
        chunks[[cn]] <- arr
      }  # else it's the block, blockhead.
    }
    tile_name <- sprintf("%d_%d", block["row"], block["col"])
    save_chunks(chunks_fn, chunks, group_name = tile_name)
  }
}


read_outputs <- function(chunks_fn, only_name = NULL) {
  tile_names <- rhdf5::h5ls(chunks_fn, recursive = FALSE)[, "name"]
  outputs <- list()
  for (group in tile_names) {
    tile_chunks <- load_chunks(chunks_fn, group_name = group, only_name = only_name)
    tile_chunks$block <- attr(tile_chunks[[1]], "tile")
    stopifnot(length(tile_chunks$block) == 2)
    names(tile_chunks$block) <- c("row", "col")
    outputs[[length(outputs) + 1]] <- tile_chunks
  }
  outputs
}


worker <- function(args) {
  # We will want to split this work different ways for development,
  # laptops, and the cluster, so we use an explicit domain decomposition.
  rampdata::initialize_workflow(args$config)
  version_dir <- build_outvars_dir(args$outvars)
  plan_rp <- rampdata::add_path(version_dir, file = "plan.json")
  plan <- load_plan(rampdata::as.path(plan_rp))
  work_df <- data.table::fread(
    rampdata::as.path(rampdata::add_path(version_dir, file = "work.csv")))

  # The part that this particular process will do.
  my_work <- my_tasks(args$task, args$tasks, work_df)
  if (is.null(my_work)) {
    flog.info("Quitting because there is no work to do.")
    quit(save = "no", status = 0)
  }

  flog.debug(paste("tile count", dim(my_work)[2], "size",
                   paste(plan$tiles$blocksize, collapse = ",")))

  # By this point, someone told us what to do, so load all the data.
  process_extent <- raster_extent_from_work(my_work, plan$tiles$blocksize)
  flog.debug(paste(
    "extent", paste(process_extent, collapse = ","), "names",
    paste(names(process_extent), collapse = ","), "\n"))
  # The load_extent is what the process needs, relative to the domain coordinates.
  load_extent <- apply_relative_extent(plan$domain_extent, process_extent)
  data <- load_data(args$config, args$pr2ar, load_extent, args$years)
  # We permute data axes for more efficient running and
  # split the data into chunks for processing.
  timeseries_data <- prepare_timeseries(data, my_work, process_extent, plan$tiles$blocksize)
  data <- NULL  # Help the memory cleanup. This data is now in timeseries chunks.

  params <- timeseries_data$parameters
  # Uses the parallel random generation from the `parallel` package.
  set.seed(params$random_seed, "L'Ecuyer")
  draw_params <- draw_parameters(params, args$draws)
  chunks <- timeseries_data$chunks

  core_cnt <- parallel_core_cnt(args)
  flog.debug("%d chunks to %d cores", length(chunks), core_cnt)
  # Work in parallel on this machine.
  # Thread output usually goes to /dev/null, but we can redirect
  # it to a file. Put it with our data for now. Could be a mess
  # for thread contention, but they don't write a lot.
  thread_out_dir <- build_outvars_dir(args$outvars)
  thread_out <- rampdata::add_path(thread_out_dir, file = "threads.txt")
  thread_fn <- rampdata::as.path(thread_out)

  output <- lapply(
    chunks,
    function(chunk) {
      flog.debug(".")
      over_block_draw(
        chunk,
        draw_params,
        timeseries_data$pr_to_ar_dt,
        params$confidence_percent
      )
    }
  )
  chunks_rp <- rampdata::add_path(version_dir, file = sprintf("chunks%d.h5", args$task))
  chunks_fn <- rampdata::as.path(chunks_rp)
  save_outputs(chunks_fn, output)
}


assemble <- function(args) {
  # We will want to split this work different ways for development,
  # laptops, and the cluster, so we use an explicit domain decomposition.
  rampdata::initialize_workflow(args$config)
  options <- configr::read.config(args$config)[["options"]]
  version_dir <- build_outvars_dir(args$outvars)
  plan_rp <- rampdata::add_path(version_dir, file = "plan.json")
  plan <- load_plan(rampdata::as.path(plan_rp))
  stopifnot(is_tile(plan$tiles$blocksize))

  task_name_fn <- function(task_idx) {
    chunks_rp <- rampdata::add_path(version_dir, file = sprintf("chunks%d.h5", task_idx))
    rampdata::as.path(chunks_rp)
  }

  missing <- list()
  for (miss_idx in 1:args$tasks) {
    find_fn <- task_name_fn(miss_idx)
    if (!file.exists(find_fn)) {
      missing[as.character(miss_idx)] <- find_fn
    }
  }
  if (length(missing) > 0) {
    flog.error(paste(
      sprintf("Cannot find the following %d tasks:", length(missing)),
      paste0(names(missing), collapse = ",")
      ))
    for (miss in missing) {
      flog.error(miss)
    }
    stopifnot(length(missing) == 0)
  }

  ds_names <- read_dataset_names(task_name_fn(1))
  for (ds_name in ds_names) {
    output <- list()
    for (task_idx in 1:args$tasks) {
      more_output <- read_outputs(task_name_fn(task_idx), ds_name)
      output <- c(output, more_output)
    }
    flog.info(sprintf("loaded %d chunks", length(output)))
    # The output chunks need to be reassembled before writing.
    ready_to_write <- combine_output(output, plan$tiles$blocksize, ds_name)
    write_output(ready_to_write, args$years, plan$domain_dimensions, args, options)
  }
}
