#' @title Mean Diff
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
calculate_Mult_diff <- function(df) {
    year<-unique(df$time)
    
    #makes an empty dataframe
    Mult_Comp_Disp_Diff=data.frame() 
    ##calculating bray-curtis mean change and disperison differecnes
    for(i in 1:length(year)) {
        
        #subsets out each dataset
        subset<-df%>%
            filter(time==year[i])%>%
            select(treatment, time, species, abundance, replicate, C_T)
        
        #need this to keep track of plot mani
        labels=subset%>%
            select(C_T, treatment)%>%
            unique()
        
        #transpose data
        species=subset%>%
            spread(species, abundance, fill=0)
        
        #calculate bray-curtis dissimilarities
        bc=vegdist(species[,5:ncol(species)], method="bray")
        
        #calculate distances of each plot to treatment centroid (i.e., dispersion)
        disp=betadisper(bc, species$treatment, type="centroid")
        
        #getting distances between centroids over years; these centroids are in BC space, so that's why this uses euclidean distances
        cent_dist=as.data.frame(as.matrix(vegdist(disp$centroids, method="euclidean")))
        
        #extracting only the distances we need and adding labels for the comparisons;
        cent_C_T=data.frame(time=year[i],
                            treatment=row.names(cent_dist),
                            mean_change=t(cent_dist[names(cent_dist)==labels$treatment[labels$C_T=="Control"],]))
        
        #renaming column
        colnames(cent_C_T)[3]<-"comp_diff"
        
        #collecting and labeling distances to centroid from betadisper to get a measure of dispersion and then take the mean for a treatment
        disp2=data.frame(time=year[i],
                         treatment=species$treatment,
                         C_T=species$C_T,
                         replicate=species$replicate,
                         dist=disp$distances)%>%
            tbl_df%>%
            group_by(time, treatment, C_T)%>%
            summarize(dispersion=mean(dist))
        
        control<-disp2$dispersion[disp2$C_T=="Control"]
        
        ##subtract control from treatments
        disp_treat=disp2%>%
            mutate(disp_diff=dispersion-control)%>%
            select(-dispersion)
        
        #merge together change in mean and dispersion data
        distances<-merge(cent_C_T, disp_treat, by=c("time","treatment"))
        
        #pasting dispersions into the dataframe made for this analysis
        Mult_Comp_Disp_Diff=rbind(Mult_Comp_Disp_Diff, distances)  
    }
    
    return(Mult_Comp_Disp_Diff)
}



############################################################################
#
# Private functions: these are internal functions not intended for reuse.
# Future package releases may change these without notice. External callers
# should not use them.
#
############################################################################
