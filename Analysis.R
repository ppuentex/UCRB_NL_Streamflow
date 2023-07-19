rm(list=ls())
source("./function_library.r")


#output directories
datadir = "./Output_Data/"


#Select the column of data that will be used for the reconstructed 
#this is always the second column for all data from treeflow
recon_col=2 

#window needed for rolling average, rolling variance to match avg LLE 
avg_window = 37


#Threshold for False nearest neighbors for dimension picking
rt = 10 #threshold for declaring false neighbors - if distance between neighbors exceeds rt of that in the previous dimension they are FNN
fnn_th=10 #%FNN allowed

# evolution periods that will be use for Lyapunov exponents
scale=c(1,2,4,8,16,20,32,64)
scale_pick=20
nscale=length(scale)
scalei=which(scale==20)


names = c("greenriverwy", "greenriverut", "glenwood", "gunnison", "cisco", "leesferry-short", "leesferry-long")
locations = c("greenriverwy","greenriverut","COglenwood","gunnisonriver", "COcisco", "COleesmeko", "COleesmeko")
startY_list = c(1569, 1569, 1569, 1569, 1569, 1569, 762)
i_loc = c(6,7,3,4,5,2,2) #column index of location in observed data 
upper_sig_list = c(0.77, 0.73, 0.61, 0.61, 0.68, 0.92, 0.95)
lower_sig_list = c(0.7, 0.65, 0.55, 0.55, 0.6, 0.85, 0.9)

#create dataframe to store the embedding parameters 
embed_table = data.frame()
stats_table = data.frame()

#initiate counter to use as index for each list
counter = 0
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
  
}
print("All data files have been processed")
