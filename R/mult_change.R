#' @title Mult Change
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
calculate_Mult_change <- function(df) {
    #make a new dataframe with just the label;
    replist<-unique(df$replicate)
    
    #makes an empty dataframe
    Mult_Comp_Disp_Change=data.frame(replicate=c(), time=c(), compositional_change=c(), dispersion_change=c()) 
    
    #Calculating bc mean change and dispersion
    for(i in 1:length(replist)) {
        
        #subsets out each dataset
        subset=df%>%
            filter(replicate==replist[i])  
        #get years
        timestep<-sort(unique(subset$time))
        #transpose data
        species=subset%>%
            spread(species, abundance, fill=0)
        
        #calculate bray-curtis dissimilarities
        bc=vegdist(species[,5:ncol(species)], method="bray") #this is 4 here b/c there is a treatment column, but if obs data only would be ,3:ncol()
        
        #calculate distances of each plot to year centroid (i.e., dispersion)
        disp=betadisper(bc, species$time, type="centroid")
        
        #getting distances between centroids over years; these centroids are in BC space, so that's why this uses euclidean distances
        cent_dist=as.matrix(vegdist(disp$centroids, method="euclidean"))
        
        ##extracting only the comparisions we want year x to year x=1.
        ###(experiment_year is year x+1
        cent_dist_yrs=data.frame(replicate=replist[i],
                                 time=timestep[2:length(timestep)],
                                 mean_change=diag(cent_dist[2:nrow(cent_dist),1:(ncol(cent_dist)-1)]))
        
        #collecting and labeling distances to centroid from betadisper to get a measure of dispersion and then take the mean for a year
        disp2=data.frame(replicate=replist[i],
                         time=species$time,
                         dist=disp$distances)%>%
            tbl_df%>%
            group_by(replicate, time)%>%
            summarize(dispersion=mean(dist))
        
        ##subtract consequtive years subtracts year x+1 - x. So if it is positive there was greater dispersion in year x+1 and if negative less dispersion in year x+1
        disp_yrs=data.frame(replicate=replist[i],
                            time=timestep[2:length(timestep)],
                            dispersion_diff=diff(disp2$dispersion))
        
        #merge together change in mean and dispersion data
        distances<-merge(cent_dist_yrs, disp_yrs, by=c("replicate","time"))
        
        #pasting dispersions into the dataframe made for this analysis
        Mult_Comp_Disp_Change=rbind(distances, Mult_Comp_Disp_Change)  
    }
    return(Mult_Comp_Disp_Change)
}


############################################################################
#
# Private functions: these are internal functions not intended for reuse.
# Future package releases may change these without notice. External callers
# should not use them.
#
############################################################################
