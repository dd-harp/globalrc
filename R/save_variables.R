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

  tmap::tmap_options(check.and.fix = TRUE)
  if (kind %in% names(.plot.kinds)) {
    aplot <- tmap::tm_shape(raster_obj) +
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
      tmap::tm_shape(admin0) +
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
     grDevices::dev.off()
  }
}


add_source_scripts <- function(dest_dir) {
  source_files <- list.files(system.file("testdata/R", package="globalrc"),
                             full.names = TRUE)
  if (length(source_files) > 0) {
    dest_r <- file.path(dest_dir, "R")
    dir.create(dest_r, showWarnings = FALSE)
    file.copy(source_files, dest_r, overwrite = TRUE)
  }
}


#' Given outputs for a tile, write them to a folder.
#' @param output A list of results.
#' @param years Years that will save.
#' @param domain_dimensions Resolution of the raster and bounds in lat-long.
#' @param domain_extent Bounds and resolution of the raster.
#' @param args Command-line arguments that have been parsed.
#' @param options Settings from the toml file.
write_output <- function(output, years, domain_dimensions, domain_extent, args, options) {
    flog.debug(paste("writing output", sum(is.na(output[[1]])),
        "na values in first array and",
        sum(output[[1]] > 0 & output[[1]] < 1, na.rm = TRUE),
        "in 0 < x < 1."
    ))

    outline_rp <- rampdata::ramp_path("/inputs/country_outlines/201122")
    # Use capture to quiet the text that sf spews.
    capture.output({admin0 <- sf::st_read(rampdata::as.path(rampdata::add_path(
        outline_rp, file = "ne_10m_admin_0_countries_lakes")))})

    dest_dir <- build_outvars_dir(args$outvars)
    # The data could be larger than the domain because the tiles don't fit evenly.
    # This is where we cut it down.
    row_cnt <- unname(domain_extent["rmax"] - domain_extent["rmin"] + 1)
    col_cnt <- unname(domain_extent["cmax"] - domain_extent["cmin"] + 1)
    for (name in names(output)) {
        if (!name %in% names(output)) {
            msg <- paste("output doesn't have the", name, "array")
            flog.error(msg)
            stop(msg)
        }
        by_year <- aperm(output[[name]], c(2, 3, 1))
        for (year_idx in 1:length(years)) {
            year <- years[year_idx]
            out_rp <- rampdata::add_path(dest_dir, file = paste0(name, "_", year, ".tif"))
            tile_sized_data <- by_year[, , year_idx]
            ready_data <- tile_sized_data[1:row_cnt, 1:col_cnt]
            out_fn <- rampdata::as.path(out_rp)
            exists <- file.exists(out_fn)
            if (exists & args$overwrite) {
                unlink(out_fn)
                exists <- FALSE
            }
            if (!exists) {
                flog.info(paste("writing file", out_fn))
                rampdata::prov.output.file(out_fn, "rc")
                raster_obj <- raster::raster(
                    nrows = row_cnt, ncols = col_cnt,
                    xmn = domain_dimensions$xmin,
                    xmx = domain_dimensions$xmax,  # columns are x for raster.
                    ymn = domain_dimensions$ymin,
                    ymx = domain_dimensions$ymax,
                    crs = domain_dimensions$projection
                )
                raster_obj <- raster::setValues(raster_obj, ready_data)
                raster::writeRaster(raster_obj, filename = out_fn, format = "GTiff")

                png_rp <- rampdata::add_path(out_rp, file = sprintf("%s_%d.png", name, year))
		            png_fn <- rampdata::as.path(png_rp)
                flog.info(paste("writing file", png_fn))
                tryCatch({
                  plot_as_png(raster_obj, png_fn, name, year, options, admin0)
                  rampdata::prov.output.file(png_fn, "image")
                },
                error = function(cond) {
                  # Don't let a failed plot stop us from saving results.
                  message(paste("Could not plot", png_fn))
                  message(cond)
                })
            } else {
                flog.error(paste("cannot overwrite", out_fn))
            }
        }
    }
    file_io_fn <- rampdata::as.path(rampdata::add_path(
        dest_dir, file = "roles.toml"))
    flog.info(paste("Writing provenance toml to", file_io_fn))
    rampdata::write.meta.data(file_io_fn)
    add_source_scripts(rampdata::as.path(dest_dir))
}
