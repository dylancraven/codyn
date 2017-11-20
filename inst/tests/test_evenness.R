context("evenness")

test_that("evenness loads and returns correct result", {
    # Load our example data set
    data(knz_001d)
    expect_that(names(knz_001d)[4], equals("abundance"))

    # Test if evenness produces correct values with count data
    # Some abundance values are 0, and some are missing
    x <- c(2,3,4,5,8,9,11,0,0,23,11,2,1,4)
    result <- evenness(x)
    expect_that(length(result), equals(1))
    expect_that(class(result), matches("numeric"))
    expect_that(result, equals(0.2033371, tolerance=0.0001))
    
    # Test if evenness produces NA with only 1 species
    x <- c(99)
    result <- evenness(x)
    expect_that(result, equals(NA))
    
    # Test if evenness produces 1 in a perfectly even community
    x <- c(15,15,15)
    result <- evenness(x)
    expect_that(length(result), equals(1))
    expect_that(class(result), matches("numeric"))
    expect_that(result, equals(1))
    
    # Test if richness produces a error with non-numeric data
    expect_error(evenness(unique(iris$Species)))
    
})