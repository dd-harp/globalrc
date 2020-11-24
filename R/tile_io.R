
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
    flog.debug(sprintf(
      "save %d chunks to tile %s in file %s",
      length(chunks), tile_name, chunks_fn
    ))
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
