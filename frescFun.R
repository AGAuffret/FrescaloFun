### Frescalo extra FUNctions ###
## Alistair Auffret - April 2024 ##

#---------- SUMMARY -----------#

# This file contains some functions that take outputs from the 'sparta' function 'frescalo' to provide some additional metrics

# 1. frescImport(folder) : Takes a folder containing frescalo outputs and makes list that can be used with rest of the functions. The list resembles the output list from the frescalo function, but also brings in some additional files that are also needed.

# 2. frescS_it(frescalo.results): Takes the results list from the previous function and calculates recorder effort (proportion of neighbourhood benchmark species recorded in the target grid cell)

# 3. frescP_ijt(frescalo.results): Takes the results list from the initial function and calculates probability of occurrence for each species in each grid cell in each time period. Currently just gives values where species have a P_ijt - might be useful to also include rows with NAs where species not observed in a grid-square's neighbourhood?

# 4. frescTrends(): Takes the results list and calculates different trend metrics, for example difference in relative occupancy (a la Fox), and uncertainty of trends from linear models (a la Pescott). [not there yet: Also number of grid cells in different time periods and number of 'extirpated' grid cells?]

# 5. frescPescPlot(): Takes the output from frescTrends to plot species a la Pescott.

#---------- frescImport -----------#

frescImport <- function(folder){
	if(is.null(folder)){stop("Please supply a folder name")}
	if(!file.exists(folder)){stop("Folder not found")}
	
fresc.out <- list(
paths=folder, # in sparta, this is locations to all the output files. We just put the whole folder here. 
trend=read.csv(paste0(folder,"/Output/Trend.csv")),
stat=read.csv(paste0(folder,"/Output/Stats.csv")),
freq=read.csv(paste0(folder,"/Output/Freq.csv")),
lm_stats=read.csv(paste0(folder,"/Maps_Results/Frescalo Tfactor lm stats.csv")),
log_file=readLines(paste0(folder,"/Output/Log.txt")),
in_data=read.table(paste0(folder,"/Input/FocDist.txt")),
spe_codes=read.csv(paste0(folder,"/species_names.csv")))

return(fresc.out)
}



#---------- frescS_it -----------#

frescS_it <- function(frescalo.results){

# Simple check to see if the input object has the right dimensions
if(!class(frescalo.results)=="list" | !length(frescalo.results)==8){stop("Please supply a list from the output of frescImport")}

# Find the benchmark limit from the log file!
log_bench <- frescalo.results$log_file[grep("Benchmark limit is", frescalo.results$log_file)] # find part of log file with benchmark info
R_ast <- as.numeric(gsub("[[:alpha:]]", "", log_bench)) # extract the number from the string

# Find the time periods
time.periods <- unique(frescalo.results$trend$Time)

# Make output data frame. Probably a nicer way, but I wanted to make it adjust according to number of time periods.
out.df <- setNames(as.data.frame(matrix(nrow=length(unique(frescalo.results$freq$Location)),ncol=length(time.periods)+1)),c("grid.square",time.periods))
out.df$grid.square <- unique(frescalo.results$freq$Location)

# Now loop for each grid square
for(sq in out.df$grid.square){
freq.sq.bench <- frescalo.results$freq$Species[frescalo.results$freq$Location==sq & frescalo.results$freq$Rank1<R_ast] # identify benchmark species for the neighbourhood
freq.sq.bench.codes <- frescalo.results$spe_codes$SPECIES[frescalo.results$spe_codes$NAME %in% freq.sq.bench] # find what the codes for those species are
in.sq <- frescalo.results$in_data[frescalo.results$in_data$V1==sq,] # subset in data to make the following more readable
out.df[out.df$grid.square==sq,2:ncol(out.df)] <- sapply(time.periods, function(x) sum(freq.sq.bench.codes %in% in.sq$V2[in.sq$V3==x])/length(freq.sq.bench)) # calculate proportion of benchmark species for each time period and add it to the output
}

names(out.df)[names(out.df) %in% time.periods]  <-  paste0("t", names(out.df)[names(out.df) %in% time.periods])
return(out.df)
}


#---------- frescP_ijt -----------#

frescP_ijt <- function(frescalo.results){

# Simple check to see if the input object has the right dimensions
if(!class(frescalo.results)=="list" | !length(frescalo.results)==8){stop("Please supply a list from the output of frescImport")}

# Find the time periods
time.periods <- unique(frescalo.results$trend$Time)

# Make output data frame. Probably a nicer way to do it...	
out.df.names <- c("species","grid.square",time.periods)
out.df <- setNames(as.data.frame(matrix(nrow=nrow(frescalo.results$freq),ncol=length(time.periods)+2)), out.df.names)
out.df[,c("species","grid.square")] <- frescalo.results$freq[order(frescalo.results$freq$Species),c("Species", "Location")]

frescalo.results$freq$Freq1[frescalo.results$freq$Freq1>=0.98] <- 0.98

out.df$fij <- frescalo.results$freq$Freq1[match(paste(out.df$species, out.df$grid.square), paste(frescalo.results$freq$Species,  frescalo.results$freq$Location))]

out.df[,paste0("xjt.",time.periods)] <- sapply(time.periods, function(x) frescalo.results$trend$TFactor[match(paste(out.df$species, x), paste(frescalo.results$trend$Species, frescalo.results$trend$Time))])

out.df[,as.character(time.periods)] <- sapply(time.periods, function(x) 1-exp(-(-log(1-out.df$fij))*out.df[,paste0("xjt.",x)]))

out.df <- out.df[,out.df.names]

names(out.df)[names(out.df) %in% time.periods]  <-  paste0("t", names(out.df)[names(out.df) %in% time.periods])

return(out.df)
		
	}



# ----------- frescTrends -------------- #

frescTrends <- function(frescalo.results, return.all=TRUE){
	
# Simple check to see if the input object has the right dimensions
if(!class(frescalo.results)=="list" | !length(frescalo.results)==8){stop("Please supply a list from the output of frescImport")}

if(!is.logical(return.all)){stop("'return.all' must be logical. TRUE gives output table and list of linear model outputs; FALSE gives the output table only")}

# Find the time periods
time.periods <- unique(frescalo.results$trend$Time)

out.df <- setNames(data.frame(matrix(nrow=length(unique(fresc.res$trend$Species)),ncol=2*(length(time.periods))+1)),c("species", paste0("rel.occ.val.",time.periods),paste0("rel.occ.sd.",time.periods)))
out.df$species <- unique(fresc.res$trend$Species)

lm.draw.list <- list() # make list for the Pescott species level lm outputs

for(spe in out.df$species){

	# Extract the Frescalo trend info for this species just to clean up following code
	spe.trend <- frescalo.results$trend[frescalo.results$trend$Species==spe,]

## Fox stuff ##

	out.df[out.df$species==spe,paste0("rel.occ.val.",time.periods)] <- fresc.res$trend$TFactor[fresc.res$trend$Species==spe & fresc.res$trend$Time==time.periods] # For all time periods, add the time factors
	out.df[out.df$species==spe,paste0("rel.occ.sd.",time.periods)] <- fresc.res$trend$StDev[fresc.res$trend$Species==spe & 	fresc.res$trend$Time==time.periods] # add the sd for time factors

	# Then just for first and last time period, do the differences with the z scores
	rrr1 <- spe.trend$TFactor[spe.trend$Time==time.periods[1]]
	rrr2 <- spe.trend$TFactor[spe.trend$Time==time.periods[length(time.periods)]]
	rrr1.sd <- spe.trend$StDev[spe.trend$Time==time.periods[1]]
	rrr2.sd <- spe.trend$StDev[spe.trend$Time==time.periods[length(time.periods)]]

	out.df[out.df$species==spe,"rel.occ.change"] <- rrr2-rrr1
	out.df[out.df$species==spe,"rel.occ.change.per.year"] <- (rrr2-rrr1)/diff(range(time.periods))
	out.df[out.df$species==spe,"z.val"] <- (rrr2-rrr1)/sqrt(rrr1.sd^2+rrr2.sd^2)
	out.df[out.df$species==spe,"p.val"] <- 2*pnorm(abs(rrr2-rrr1)/sqrt(rrr1.sd^2+rrr2.sd^2),lower.tail=FALSE)


## Pescott stuff ##

	spe.time <- data.frame(species=rep(spe,length(rep(time.periods,100))),time.period=rep(time.periods,100)) # make data frame
	tf.draws <- sapply(spe.trend$Time, function(x) rnorm(100, mean=spe.trend$TFactor[spe.trend$Time==x], sd=spe.trend$StDev[spe.trend$Time==x])) # make table with the different draws
	spe.time$tf.draws <- c(t(tf.draws)) # simplify and rearrange (so that years alternate), and add to data frame
	spe.time.ls <- split(spe.time,rep(1:100, each=length(time.periods))) # split data frame into list of mini data frames
	spe.time.lm <- lapply(spe.time.ls, function(x) lm(tf.draws~time.period, data=x)$coefficients) # perform lm on those
	spe.lm.coeff <- setNames(data.frame(do.call(rbind,spe.time.lm)), c("intercept","estimate")) # collect all the parameter estimates
	out.df[out.df$species==spe, c("Strong decline", "Moderate decline", "Stable", "Moderate increase", "Strong increase")] <- table(cut(spe.lm.coeff$estimate,c(1,0.004,0.001,-0.001,-0.004,-1)))/100 # Then cut these according to Pescott and add them to the output!

lm.draw.list[[spe]] <- spe.lm.coeff # Save these results in a list for plotting later

}

if(!length(unique(frescalo.results$trend$Time))==2){warning("Trend estimates based on z-scores only take into account first and last time period")}

if(return.all==TRUE){return(list(trends=out.df, lm.coeff=lm.draw.list, time.periods=time.periods))}
if(return.all==FALSE){return(out.df)}

}


#---------- frescPescPlot -----------#

frescPescPlot <- function(species, trends, point.col="black", line.col="forestgreen"){
  # Simple check to see if the input object has the right dimensions
  if(!class(trends)=="list" | !length(trends)==3){stop("Please supply a list from the output of frescTrends(return.all=TRUE)")}
  
  # Check that the species is there
  if(length(species)>1){stop("Sorry, this function only plots one species at a time")}
  
    # Check that the species is there
  if(!species %in% trends$trends$species){stop(paste("The species", species, "is not in the trends data frame"))}
  
  # Check colours
  if(!point.col %in% colors()) {stop(paste(point.col, "is not a colour - this function only accepts namned colours"))}
  if(!line.col %in% colors()) {stop(paste(line.col, "is not a colour - this function only accepts namned colours"))}
  
  # Get vectors of the relative occupancy and standard deviations
  rel.occ.spe  <-  sapply(trends$time.periods, function(x) trends$trends[trends$trends$species==species,paste0("rel.occ.val.",x)] )
  rel.occ.sd.spe <- sapply(trends$time.periods, function(x) trends$trends[trends$trends$species==species,paste0("rel.occ.sd.",x)] )
  
  # make axis limits relevant to the species at hand
  ylims <- c(mean(rel.occ.spe)-(3*max(rel.occ.sd.spe)),mean(rel.occ.spe)+(3*max(3*rel.occ.sd.spe)))
  xlims <- c(min(trends$time.periods)-10, max(trends$time.periods)+10)
  
  # plot it  
  plot(0,xlim=xlims,ylim=ylims,type="n", axes=TRUE, frame.plot=FALSE, main=species, xlab="Year", ylab="Relative occupancy", font.axis=2, font.lab=2, cex.lab=1.2, cex.main=1.3) # create empty plot
  apply(trends$lm.coeff[[species]],1,abline,col=adjustcolor(line.col,alpha.f = 0.3) )
  points(trends$time.periods,rel.occ.spe, pch=16, col=point.col,cex=1.75) # create empty plot
  segments(trends$time.periods, rel.occ.spe-rel.occ.sd.spe, trends$time.periods, rel.occ.spe+rel.occ.sd.spe, lwd=2.5, col=point.col) # draw lines for confidence intervals

}
