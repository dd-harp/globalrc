#' Loads a script that has functions to compute pixel values.
#' @param script_name The file with the script.
#' @value An environment containing script functions.
#' @export
load_pixel_script <- function(script_name) {
  contain <- new.env(parent = .BaseNamespaceEnv)
  source(script_name, contain)
  contain
}
