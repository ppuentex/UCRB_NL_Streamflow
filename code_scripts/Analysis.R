rm(list=ls())
source("./code_scripts/function_library.r")
library(fractal) #needed to use timeLag func 
library(zoo)


#output directories
datadir = "./Output_Data/"


#Select the column of data that will be used for the reconstructed 
#this is always the second column for all data from treeflow
recon_col=2 

#window needed for rolling average, rolling variance to match avg LLE 
#roll_window = 39

scale = 1e6 #want to convert everything to MAF


#Threshold for False nearest neighbors for dimension picking
rt = 10 #threshold for declaring false neighbors - if distance between neighbors exceeds rt of that in the previous dimension they are FNN
fnn_th=10 #%FNN allowed

# evolution periods that will be use for Lyapunov exponents
exp_scale=c(1,2,4,8,16,20,32,64)
scale_pick=20
nscale=length(exp_scale)
scalei=which(exp_scale==scale_pick)


names = c("greenriverwy", "greenriverut", "glenwood", "gunnison", "cisco", "leesferry-short", "leesferry-long")
locations = c("greenriverwy","greenriverut","COglenwood","gunnisonriver", "COcisco", "COleesmeko", "COleesmeko")
startY_list = c(1569, 1569, 1569, 1569, 1569, 1569, 762)
i_loc = c(6,7,3,4,5,2,2) #column index of location in observed data 
upper_sig_list = c(0.77, 0.73, 0.61, 0.61, 0.68, 0.92, 0.95)
lower_sig_list = c(0.7, 0.65, 0.55, 0.55, 0.6, 0.85, 0.9)

#create dataframe to store the embedding parameters 
stats_table = data.frame()

#initiate counter to use as index for each list
counter = 0
#location = locations[7]
for(location in locations){
  counter = counter + 1 
  
  file = paste("./Input_Data/",location,".csv", sep="")
  raw_obs = read.csv("./Input_Data/MultiObserved_data.csv", header = TRUE)
  obs_col = i_loc[counter]
  raw_rec=read.csv(file, header = TRUE)
  
  #rename columns so they are consistent for binding
  colnames(raw_rec)[1]='year'
  colnames(raw_rec)[recon_col]='flow'
  colnames(raw_obs)[1]='year'
  colnames(raw_obs)[obs_col]='flow'
  
  startR = startY_list[counter]
  istartR = which(raw_rec[,1]==startR)
  endR = 1905 #end year for recon
  iendR = which(raw_rec[,1]==endR)
  startO = raw_obs[1,1] #start obs year
  istartO = which(raw_obs[,1]==startO)
  endO = raw_obs[length(raw_obs[,1]),1]
  iendO = which(raw_obs[,1]==endO)
  raw_data = rbind(raw_rec[istartR:iendR, c(1, recon_col)], raw_obs[istartO:iendO, c(1,obs_col)])
  
  #simple statistics for each location
  y_range = paste(startR,"-",endO)
  y_num = length(raw_data[,1])
  raw_avg = round(mean(raw_data[,2])/scale,2)
  raw_max = round(max(raw_data[,2])/scale,2)
  raw_min = round(min(raw_data[,2])/scale,2)
  raw_range = paste(raw_min,"-",raw_max)
  
  
  #define wavelet filtering 
  upper_sig = upper_sig_list[counter]
  lower_sig = lower_sig_list[counter]
  
  #Step 1: wavelet transform of streamflow 
  data_wavelet=decompose2(raw_data[,2]/scale, sig=upper_sig, dj=0.025)
  data_signal = data_wavelet$sig
  raw_data$signal <- data_signal
  
  #Step 2: Time delay embedding 
  # Calculate time lag using mutual information
  Tau=timeLag(data_signal,method='mutual',plot.data=F)
  mi=attr(Tau,'data')  # this is the mutual information 
  lags=attr(Tau,'lags')
  lag_pick=which(diff(sign(diff(mi)))==2)[1]+1 #pick the first local minima 
  
  
  #Step 3: Calculate the embedding dimension using false nearest neighbors
  fnn=FNN(data_signal, dimension=10, tlag=lag_pick, rtol=rt, olag=1)[1,]
  dim_pick=as.numeric(which(fnn<fnn_th)[1]) #Pick first dimension where FNN falls below the fnn_threshold
  
  
  #Step 4: Calculate Average Lyaponov exponent 
  
  #figure out how much of the signal can be used based on longest evolution period
  #to be tested
  n.reference=10
  Ne=length(data_signal)-(lag_pick*(dim_pick-1)) #number of evaluation points 
  avgmax.ref=Ne-scale_pick-n.reference-2
  avgLE_spectrum <- lyapunov(data_signal, dimension=dim_pick, tlag=lag_pick,
                             scale=scale_pick, local.dimension=dim_pick, reference=(1:(avgmax.ref)),
                             n.reference=n.reference)
  
  avgLE_val_mat = matrix(unlist(avgLE_spectrum),nrow=avgmax.ref)
  avg.LL=rowMeans(avgLE_val_mat)
  avg_write = matrix(NA, nrow = length(raw_data[,2]), ncol = 1)
  #calculates how far from the starting year in order to be centered
  dist = round(((length(raw_data[,2]) - length(avg.LL))/2), digit = 0)
  
  avg_write[(dist+1):(dist+length(avg.LL))] = avg.LL
  
  raw_data$avgLLE <- avg_write
  
  #find the difference in length of avg LLE and SF
  #avgLLE_dif = length(raw_data[,2]) - length(avg.LL)

  #need to figure this out 
  roll_window = scale_pick + n.reference + (lag_pick*(dim_pick-1)) +2 + 1
  #2 is subtracted in avg.ref and 1 is added to get additional year
  
  #calculate a centered rolling mean and rolling variance
  raw_data$flow_mean <- rollapply((raw_data[,2]/scale), width = roll_window, FUN = mean, align = "center", fill=NA)
  raw_data$flow_var <- rollapply((raw_data[,2]/scale), width=roll_window, FUN = var , align="center", fill = NA)
  
  #mean.diff = length(raw_data[,2]) - length(na.omit(raw_data$flow_mean))
  #var.diff = length(raw_data[,2]) - length(na.omit(raw_data$flow_var))
  
  id = names[counter]
  
  
  stats_output = c(names[counter], y_range, y_num, raw_avg, raw_range, 
                   lag_pick, dim_pick, roll_window)
                   #avgLLE_dif, mean.diff, var.diff, roll_window)
  
  stats_table = rbind(stats_table, stats_output)
  
  write.csv(raw_data, file = paste(datadir,id,"_ts.csv", sep=""), row.names = F)
  
}
#prints table of embedding parameters  for each location 


colnames(stats_table) <- c("flow ts" ,"Year Range", "# of Years", "Mean (MAF)", 
                           "Flow Range (MAF)", "Time Lag", "Dim embedding", 
                           #"avg LLE diff", "mean diff ", "var diff", "roll window")
                           "roll window")
stats_table
write.csv(stats_table, file=paste(datadir,"ts_info.csv", sep=""),row.names=F)
