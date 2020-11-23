test_that("gadm_country_shapefile finds a country", {
    afile <- gadm_country_shapefile("uga")
    expect_equal(class(afile)[1], "sf")
})
