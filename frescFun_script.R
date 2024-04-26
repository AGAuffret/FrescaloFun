### Frescalo extra FUNctions - script ###
## Alistair Auffret - April 2024 ##

# import functions file
source("...frescFun.R")

# First just choose your folder - usually would have a name like frescalo_210805(1)
out.path<-"my path"

# 1. Make your frescalo results list
fresc.res <- frescImport(out.path)
lapply(fresc.res, head) # have a peek

# 2. Calculate sampling effort (proportion of benchmark species) per sampling period
fresc.samp<-frescS_it(fresc.res)
head(fresc.samp)

# 3. Calculate species-grid-time probability of occurrence
t0<-Sys.time()
fresc_prob_vec <- frescP_ijt(fresc.res) # can take some seconds
difftime(Sys.time(), t0)

#identical(fresc_prob_vec, fresc_prob)

# 4. Produce table with trend values
fresc_trends <- frescTrends(fresc.res, return.all=TRUE) # can take some seconds
head(fresc_trends$trends)
head(fresc_trends$time.periods)
lapply(fresc_trends$lm.coeff, head)


# 5. Make a plot

frescPescPlot("Utricularia minor", fresc_trends, point.col="black", line.col="forestgreen")
