#' @title Curve Diff
#' @description A measure of the relative change in species rank abundances, which indicates shifts in relative abundances over time (Collins et al. 2008).
#' Mean rank shifts are calculated as the average difference in species' ranks between consecutive time periods, among species that are present across the entire time series.
#' @param df A data frame containing time, species and abundance columns and an optional column of replicates
#' @param time.var The name of the time column
#' @param species.var The name of the species column
#' @param abundance.var The name of the abundance column
#' @param replicate.var The name of the optional replicate column
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
#'  # Calculate mean rank shifts within replicates
#'  data(knz_001d)
#'
#'  myoutput <- rank_shift(knz_001d,
#'                      time.var = "year",
#'                      species.var = "species",
#'                      abundance.var = "abundance",
#'                      replicate.var = "subplot")
#'
#'  # Calculate mean rank shifts for a data frame with no replication
#'
#'  myoutput_singlerep <- rank_shift(subset(knz_001d, subplot=="A_1"),
#'                            time.var = "year",
#'                            species.var = "species",
#'                            abundance.var = "abundance")
#' @export
CurveDiff <- function(df) {
    ct_rank <- add_ranks_for_non_present_and_present_species_and_average_to_singe_sp_pool_per_treatment(df)
    curve_change_diff <- calculate_curve_diff(ct_rank)
    return(curve_change_diff)
}



############################################################################
#
# Private functions: these are internal functions not intended for reuse.
# Future package releases may change these without notice. External callers
# should not use them.
#
############################################################################

calculate_curve_diff <- function(rank) {
    curvechange_diff=data.frame(treatment=c(), time=c(), CC=c())#expeiment year is year of timestep2
    
    ranks<-rank%>%
        filter(abundance!=0)%>%
        group_by(time, treatment, C_T)%>%
        mutate(maxrank=max(rank),
               relrank=rank/maxrank)%>%
        arrange(-abundance)%>%
        mutate(cumabund=cumsum(abundance))%>%
        ungroup()
    
    control<-ranks%>%
        filter(C_T=="Control")
    
    treat_list<-unique(subset(ranks, C_T=="Treatment")$treatment)
    
    for(i in 1:length(treat_list)) {#minus 1 will keep me in year bounds NOT WORKING
        subset_trt<-ranks%>%
            filter(treatment==treat_list[i])
        
        #dataset of two treatments    
        subset_t12<-rbind(control, subset_trt)
        
        result <- subset_t12 %>%
            group_by(time) %>%
            do({
                y <- unique(.$treatment)###assumption this is a length 2 list
                df1 <- filter(., treatment==y[[1]])
                df2 <- filter(., treatment==y[[2]])
                sf1 <- stepfun(df1$relrank, c(0, df1$cumabund))
                sf2 <- stepfun(df2$relrank, c(0, df2$cumabund))
                r <- sort(unique(c(0, df1$relrank, df2$relrank)))
                h <- abs(sf1(r) - sf2(r))
                w <- c(diff(r), 0)
                data.frame(CC=sum(w*h))#do has to output a dataframe
            })
        
        d_output1=data.frame(treatment=treat_list[i], time=timestep[i+1], CC=result$CC)#expeiment year is year of timestep2
        
        curvechange_diff<-rbind(curvechange_diff, d_output1)
    }
    return(curve_change_diff)
}