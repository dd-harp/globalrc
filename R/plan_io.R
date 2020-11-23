library(configr)


save_plan <- function(plan, plan_file) {
  plan_list <- list(
    domain_extent = as.list(plan$domain_extent),
    domain_dimensions = plan$domain_dimensions,
    blocksize = as.list(plan$tiles$blocksize)
  )
  nullity <- vapply(plan, function(x) is.null(x), logical(1))
  if(any(nullity)) {
    was_null <- names(plan)[nullity]
    flog.error(sprintf("Could not save plan because entry %s is null",
                       paste0(was_null, collapse = ",")))
    stop()
  }
  configr::write.config(plan_list, file.path = plan_file, write.type = "json")
}


load_plan <- function(plan_file) {
  plan <- configr::read.config(plan_file)
  domain_extent <- unlist(plan$domain_extent)
  names(domain_extent) <- names(plan$domain_extent)
  plan$domain_extent <- domain_extent
  blocksize <- unlist(plan$blocksize)
  names(blocksize) <- names(plan$blocksize)
  plan$tiles <- list(blocksize = blocksize)
  plan$blocksize <- NULL
  plan
}
