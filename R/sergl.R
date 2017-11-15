#' @title SERGL
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
calculate_SERGL <- function(rank) {
    SERGL=data.frame(replicate=c(), time=c(), S=c(), E=c(), R=c(), G=c(), L=c()) #experiment year is year of timestep2
    
    replist<-unique(rank$replicate)
    
    for (i in 1:length(replist)){
        subset<-rank%>%
            filter(replicate==replist[i])
        replicate<-replist[i]
        
        #now get all timestep within an experiment
        timestep<-sort(unique(subset$time))    
        
        for(i in 1:(length(timestep)-1)) {#minus 1 will keep me in year bounds NOT WORKING
            subset_t1<-subset%>%
                filter(time==timestep[i])
            
            subset_t2<-subset%>%
                filter(time==timestep[i+1])
            
            subset_t12<-merge(subset_t1, subset_t2, by=c("species","replicate"), all=T)%>%
                filter(abundance.x!=0|abundance.y!=0)
            #reordering
            MRSc<-mean(abs(subset_t12$rank.x-subset_t12$rank.y))/nrow(subset_t12)
            #ricness and evenness differences
            s_t1 <- S(subset_t12$abundance.x)
            e_t1 <- E_q(as.numeric(subset_t12$abundance.x))
            s_t2 <- S(subset_t12$abundance.y)
            e_t2 <- E_q(as.numeric(subset_t12$abundance.y))
            
            sdiff<-abs(s_t1-s_t2)/nrow(subset_t12)
            ediff<-abs(e_t1-e_t2)/nrow(subset_t12)
            
            #gains and losses
            subset_t12$gain<-ifelse(subset_t12$abundance.x==0, 1, 0)
            subset_t12$loss<-ifelse(subset_t12$abundance.y==0, 1, 0)
            
            gain<-sum(subset_t12$gain)/nrow(subset_t12)
            loss<-sum(subset_t12$loss)/nrow(subset_t12)
            
            metrics<-data.frame(replicate=replicate, time=timestep[i+1], S=sdiff, E=ediff, R=MRSc, G=gain, L=loss)#spc_id
            ##calculate differences for these year comparison and rbind to what I want.
            
            SERGL=rbind(metrics, SERGL)  
        }
    }
    return(SERGL)
}

SERGL_func <- function(df) {
    rank <- add_ranks_for_non_present_and_present_species(df)
    result <- calculate_SERGL(rank)
    return(result)
}



############################################################################
#
# Private functions: these are internal functions not intended for reuse.
# Future package releases may change these without notice. External callers
# should not use them.
#
############################################################################

#' @param df dataset with columns for time, replicate, species, and abundance, (and optionally an id column for grouping and second column for defining the groupins)
add_ranks_for_non_present_and_present_species <- function(df) {
    ##add ranks for present species
    rank_pres<-df%>%
        filter(abundance!=0)%>%
        tbl_df()%>%
        group_by(replicate, time)%>%
        mutate(rank=rank(-abundance, ties.method = "average"))%>%
        tbl_df()
    
    #adding zeros 
    replist<-unique(df$replicate)
    allsp<-data.frame()
    for (i in 1:length(replist)){
        subset <- df %>%
            filter(replicate==replist[i])%>%
            spread(species, abundance, fill=0)
        
        long<-subset%>%
            gather(species, abundance, 5:ncol(subset))###this is depend on whether there are treatment and trt-control definition columns, will either be 3 or 5. I know this is a problem.
        allsp<-rbind(long, allsp)  
    }
    
    ###make zero abundant species have the rank S+1 (the size of the species pool plus 1)
    ##pull out zeros
    zeros<-allsp%>%
        filter(abundance==0)
    ##get species richness for each year
    SpR<-group_by(allsp, replicate, time)%>%
        summarize(S=S(abundance))
    ##merge together make zero abundances rank S+1
    zero_rank<-merge(zeros, SpR, by=c("time","replicate"))%>%
        mutate(rank=S+1)%>%
        select(-S)%>%
        tbl_df()
    ##combine all
    rank<-rbind(rank_pres, zero_rank)
    
    return(rank)
}
