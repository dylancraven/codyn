#' @title Richness
#' @description A count of the number of species in a sample.
#' @param x the vector of abundances of each species in a sample
#' @return an integer
#' @details
#' For numeric data only, richness will count the number of species in vector.
#' @references
#' @examples
#'  richness(c(2,3,4,5,8,9,11,0,0,23,11,2,1,NA,NA,4)) # Numeric count data
#'  richness(subset(knz_001d, year==1983&subplot=="A_1")$abundance) # Numeric data from a datafram
#' @export
richness <- function(x) {
    x1 <- x[x!=0 & !is.na(x)]
    stopifnot(x1==as.numeric(x1))
    length(x1)
}


############################################################################
#
# Private functions: these are internal functions not intended for reuse.
# Future package releases may change these without notice. External callers
# should not use them.
#
############################################################################
