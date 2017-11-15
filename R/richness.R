#' @title Richness
#' @description A measure of the relative change in species rank abundances, which indicates shifts in relative abundances over time (Collins et al. 2008).
#' Mean rank shifts are calculated as the average difference in species' ranks between consecutive time periods, among species that are present across the entire time series.
#' @param x the vector of abundances of each species
#' @return rank_shift returns a data frame with the following columns:
#' \itemize{
#'  \item{time.var_pair: }{A factor column that returns the two time periods being compared, separated by a dash. The name of this column is the same as the time.var column in the input dataframe followed by "_pair".}
#'  \item{MRS: }{A numeric column with the mean rank shift values.}
#'  \item{replicate.var: }{A column that has same name and type as the replicate.var column, if replication is specified.}
#' }
#' @details
#' The input data frame needs to contain columns for time, species and abundance; time.var, species.var and abundance.var are used to indicate which columns contain those variables.
#' If multiple replicates are included in the data frame, that column should be specified with replicate.var. Each replicate should reflect a single experimental unit - there must be a single abundance value per species within each time point and replicate.
#' @references  Collins, Scott L., Katharine N. Suding, Elsa E. Cleland, Michael Batty, Steven C. Pennings, Katherine L. Gross, James B. Grace, Laura Gough, Joe E. Fargione, and Christopher M. Clark.  (2008) "Rank clocks and plant community dynamics." Ecology 89, no. 12: 3534-41.
#' @examples
#'  richness(c(2,3,4,5,8,9,11,0,0,23,11,2,1,NA,NA,4)) # Numeric count data
#'  richness(unique(iris$Species))  # Factor data
#'  richness(as.character(unique(iris$Species)))  # Species list
#' @export
richness <- function(x) {
    x1 <- x[x!=0 & !is.na(x)]
    length(x1)
}


############################################################################
#
# Private functions: these are internal functions not intended for reuse.
# Future package releases may change these without notice. External callers
# should not use them.
#
############################################################################
