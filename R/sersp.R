#' @title SERSP
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
calculate_SERSp <- function(df) {
    SERSp=data.frame(treatment=c(), time=c(), Sd=c(), Ed=c(), Rd=c(), spd=c())      
    timestep<-sort(unique(rank$time)) 
    
    for(i in 1:(length(timestep))){
        
        time<-rank%>%
            filter(time==timestep[i])
        
        time_id<-timestep[i]
        
        #fitler out control plots
        control<-time%>%
            filter(C_T=="Control")
        
        treat_list<-unique(subset(time, C_T=="Treatment")$treatment)
        
        for (i in 1:length(treat_list)){
            treat<-time%>%
                filter(treatment==treat_list[i])
            
            treat_id<-treat_list[i]
            
            subset_ct<-merge(control, treat, by=c("time","species"), all=T)%>%
                filter(abundance.x!=0|abundance.y!=0)
            
            MRSc_diff<-mean(abs(subset_ct$rank.x-subset_ct$rank.y))/nrow(subset_ct)
            
            spdiff<-subset_ct%>%
                filter(abundance.x==0|abundance.y==0)
            
            spdiffc<-nrow(spdiff)/nrow(subset_ct)
            
            ##eveness richness
            s_c <- S(subset_ct$abundance.x)
            e_c <- E_q(subset_ct$abundance.x)
            s_t <- S(subset_ct$abundance.y)
            e_t <- E_q(subset_ct$abundance.y)
            
            sdiff<-abs(s_c-s_t)/nrow(subset_ct)
            ediff<-abs(e_c-e_t)/nrow(subset_ct)
            
            metrics<-data.frame(treatment=treat_id, time=time_id, Sd=sdiff, Ed=ediff, Rd=MRSc_diff, spd=spdiffc)#spc_id
            ##calculate differences for these year comparison and rbind to what I want.
            
            SERSp=rbind(metrics, SERSp)  
        }
    }
    return(SERSp)
}

SERSp_func <- function(df) {
    ct_rank <- add_ranks_for_non_present_and_present_species_and_average_to_singe_sp_pool_per_treatment(df)
    result <- calculate_SERSp(ct_rank)
    return(result)
}



############################################################################
#
# Private functions: these are internal functions not intended for reuse.
# Future package releases may change these without notice. External callers
# should not use them.
#
############################################################################

#' @param sim foo dataset with columns for time, replicate, species, abundance, a treatment column for grouping replicates, and a column speficfying if the treatment is a control or a treatment
add_ranks_for_non_present_and_present_species_and_average_to_singe_sp_pool_per_treatment <- function(df) {
    
    ###add zeros and average up species pool for control and treatment plots
    wide<-df%>%
        spread(species, abundance, fill=0)
    
    ##make long and get averages of each species by treatment
    long<-wide%>%
        gather(species, abundance, 5:ncol(wide))%>%
        group_by(time, treatment, species, C_T)%>%
        summarize(abundance=mean(abundance))
    
    ##add ranks dropping zeros
    rank_pres<-long%>%
        filter(abundance!=0)%>%
        tbl_df()%>%
        group_by(time, treatment, C_T)%>%
        mutate(rank=rank(-abundance, ties.method = "average"))%>%
        tbl_df()
    
    ###make zero abundant species have the rank S+1 (the size of the species pool plus 1)
    ##pull out zeros
    zeros<-long%>%
        filter(abundance==0)
    ##get species richness for each year
    rich<-group_by(long, time, treatment, C_T)%>%
        summarize(S=S(abundance))
    ##merge together make zero abundances rank S+1
    zero_rank<-merge(zeros, rich, by=c("time", "treatment","C_T"))%>%
        mutate(rank=S+1)%>%
        select(-S)%>%
        tbl_df()
    ##combine all
    ct_rank<-rbind(rank_pres, zero_rank)
    return(ct_rank)
}
