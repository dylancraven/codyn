context("richness")

test_that("richness loads and returns correct result", {
    # Load our example data set
    data(knz_001d)
    expect_that(names(knz_001d)[4], equals("abundance"))

    # Test if richness produces correct values with count data
    # Some abundance values are 0, and some are missing
    x <- c(2,3,4,5,8,9,11,0,0,23,11,2,1,NA,NA,4)
    result <- richness(x)
    expect_that(length(result), equals(1))
    expect_that(class(result), matches("integer"))
    expect_that(result, equals(12))
    
    # Test if richness produces correct values with factor data
    result <- richness(unique(iris$Species))
    expect_that(length(result), equals(1))
    expect_that(class(result), matches("integer"))
    expect_that(result, equals(3))
    
    # Test if richness produces correct values with species lists data
    result <- richness(unique(knz_001d$species))
    expect_that(length(result), equals(1))
    expect_that(class(result), matches("integer"))
    expect_that(result, equals(84))
    
})