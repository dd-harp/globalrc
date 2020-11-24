#' This file is for saving arrays in HDF5.
#' We need to save chunks of arrays, with lots of metadata,
#' so HDF5 lets us do that and keep everything together,
#' instead of trying to parse filenames to do it.


library(futile.logger)

#' Save a list of arrays to a single file.
#' @param filename The filename to open, save, and close.
#' @param chunks A list, where each member is an array,
#'     and each array has dimnames where each dimension is
#'     named. Any attributes of those arrays will also be
#'     saved, if HDF5 will save them. I don't know the rules
#'     on what it will/won't save, but strings and arrays
#'     should be good.
#' @export
save_chunks <- function(filename, chunks, group_name = "chunks") {
  if (!file.exists(filename)) {
    base_file <- rhdf5::H5Fcreate(filename)
  } else {
    base_file <- rhdf5::H5Fopen(filename)
  }
  tryCatch({
    # parent_group <- create_parent_groups_id(base_file, data_name)
    sandbox <- rhdf5::H5Gcreate(base_file, group_name)
    tryCatch({
      array_names <- names(chunks)
      if (!is.null(array_names)) {
        if (any(vapply(array_names, nchar, integer(1))) == 0) {
          array_names <- as.character(1:length(chunks))
        }  # else use array names as is.
      } else {
        array_names <- as.character(1:length(chunks))
      }
      for(chunk_idx in seq_along(chunks)) {
        chunk_name <- array_names[chunk_idx]
        chunk <- chunks[[chunk_idx]]
        dn <- dimnames(chunk)
        attr(chunk, "dimorder") <- names(dn)
        for (dim_idx in seq_along(dn)) {
          dim_extent <- dn[[dim_idx]]
          attr_name <- names(dn)[dim_idx]
          if (!is.null(dim_extent) && length(dim_extent) > 0) {
            # save as integers if we can, but strings are fine.
            attr(chunk, names(dn)[dim_idx]) <- tryCatch(
              as.integer(dim_extent),
              warning = function(w) dim_extent
            )
          }  # else don't save the extents.
        }
        # The dimnames won't be writeable.
        dimnames(chunk) <- NULL
        rhdf5::h5write(
          chunk, sandbox, chunk_name, write.attributes = TRUE)
      }
    }, finally = {
      rhdf5::H5Gclose(sandbox)
    })
  }, finally = {
    rhdf5::H5Fclose(base_file)
  })
}


load_chunks <- function(filename, group_name = "chunks", only_name = NULL) {
  base_file <- rhdf5::H5Fopen(filename)
  tryCatch({
    # parent_group <- create_parent_groups_id(base_file, data_name)
    sandbox <- rhdf5::H5Gopen(base_file, group_name)
    tryCatch({
      directory <- rhdf5::h5ls(sandbox, recursive = FALSE, all = FALSE)
      ds_names <- directory[, "name"]
      if (!is.null(only_name)) {
        if (only_name %in% ds_names) {
          ds_names <- only_name
        } else {
          flog.error(sprintf(
            "Cannot find dataset %s in file %s. Found datasets %s",
            only_name,
            filename,
            paste0(ds_names, collapse = ", ")
          ))
          stopifnot(only_name %in% ds_names)
        }
      }
      chunks <- list()
      for (ds_name in ds_names) {
        ds <- rhdf5::h5read(sandbox, ds_name, read.attributes = TRUE)
        attrs <- attributes(ds)
        if ("dimorder" %in% names(attrs)) {
          dim_str <- attrs[["dimorder"]]
          dim_list <- vector(mode = "list", length = length(dim_str))
          names(dim_list) <- dim_str
          for (dim_idx in seq_along(dim_str)) {
            one_dim <- dim_str[dim_idx]
            if (one_dim %in% names(attrs)) {
              dim_list[[one_dim]] <- attrs[[one_dim]]
            }
          }
          dimnames(ds) <- dim_list
        }
        chunks[[ds_name]] <- ds
      }
      chunks
    }, finally = {
      rhdf5::H5Gclose(sandbox)
    })
  }, finally = {
    rhdf5::H5Fclose(base_file)
  })
}


read_dataset_names <- function(filename) {
  df <- rhdf5::h5ls(filename, recursive = 2, all = TRUE)
  unique(df[df$rank == 3, "name"])
}
