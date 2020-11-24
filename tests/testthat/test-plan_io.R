test_that("save_plan saves to yaml", {
  plan <- list(
    domain_extent = c(rmin = 1, rmax = 37, cmin = 1, cmax = 89),
    domain_dimensions = list(
      xmin = -2.0,
      xmax = 0.3,
      xres = 0.003,
      ymax = 32.1,
      ymin = 35.7,
      yres = 0.004,
      projection = "lat-long and stuff"
    ),
    tiles = list(blocksize = c(row = 32, col = 32))
  )
  zplan <- tempfile(fileext = ".json")
  save_plan(plan, zplan)
  plan2 <- load_plan(zplan)
  expect_equal(plan2$domain_extent, plan$domain_extent)
  expect_equal(plan2$domain_dimensions, plan$domain_dimensions)
  expect_equal(plan2$tiles$blocksize, plan$tiles$blocksize)
  unlink(zplan)
})
