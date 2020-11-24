#' globalrc: Calculate R_c from PfPR and treatment rates.
#'
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
#'
#' @docType package
#' @name globalrc
#' @import configr
#' @import data.table
#' @import docopt
#' @import futile.logger
#' @import parallel
#' @import pracma
#' @import rampdata
#' @import raster
#' @import rhdf5
#' @import rgdal
#' @import sf
#' @import sp
#' @import tmap
NULL
