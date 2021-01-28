test_that("gadm_country_shapefile finds a country", {
    test_requires_data()
    afile <- gadm_country_shapefile("uga")
    expect_equal(class(afile)[1], "sf")
})
