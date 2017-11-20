#' @title Evenness measured with EQ 
#' @description A measure of the evenness of a sample using EQ (See Smith and Wilson 1996). 
#' @param x the vector of abundances of each species in a sample
#' @return Evenness returns EQ evenness, a number. 
#' @details
#' For numeric data only, evenness will calculate EQ evenness of the sample. EQ = 2/(pi)*arctan(b'), where b' is the slope of the fitted line of the log of abundance versus the scaled rank of all species in a sample. Evenness is bound between 0-1. It is 1 if all abundances are equal, and NA if there is only 1 species in the sample. [NOTE TO MEGHAN: talk to ian about the equation, if I recall we think thier negative is wrong]
#' @references  Smith, B, and Wilson, J. (1996). A consumer's guide to evenness indices. Oikos 76: 70-82.
#' Avolio et al. IN PREP.
#' @examples
#' evenness(c(2,3,4,5,8,9,11,0,0,23,11,2,1,NA,NA,4)) # Numeric data
#' evenness(subset(knz_001d, year==1983&subplot=="A_1")$abundance) # Numeric data from a dataframe
#' evenness(c(3,3,3)) # will return 1, the community is perfectly even, that is, the abundance of individuals are are equally distrubited among the species
#' evenness(3)# will return NA because evenness as a concept does not apply to a community with only one species
#' @export
evenness <- function(x) {
    x1<-x[x!=0 & !is.na(x)]
    #stopifnot(x1==as.numeric(x1))#is this necessary?
    if (length(x1)==1) {
        return(NA)
    }
    if (abs(max(x1) - min(x1)) < .Machine$double.eps^0.5) { 
        return(1)
    }
    r<-rank(x1, ties.method = "average")
    r_scale<-r/max(r)
    x_log<-log(x1)
    fit<-lm(r_scale~x_log)
    b<-fit$coefficients[[2]]
    2/pi*atan(b)
}



############################################################################
#
# Private functions: these are internal functions not intended for reuse.
# Future package releases may change these without notice. External callers
# should not use them.
#
############################################################################
