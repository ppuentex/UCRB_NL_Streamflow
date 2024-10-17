rm(list=ls())
source("./code_scripts/function_library.r")
library(quantmod) #get peaks and valleys of time series 
library(scatterplot3d) #needed for phase space plots 
library(corrplot) #correlation plots 
library(ggplot2) #needed for scatterplots color spectrum 
library(ggpmisc) #use for stat_poly_eq in  scatterplots
library(patchwork) ### used for plots ie p1| p2 | p3 etc 

#input directories 
datadir = "./Output_Data/"
#output directories
graphdir = "./Figures/"

#specify columns 
yr_col = 1 
flow_col = 2
signal_col = 3
avg_LLE_col = 4
sawp_col = 5
roll_mean_col = 6
roll_var_col = 7

#parameters 
high_pred_col = "#99CCFF"
low_pred_col = "#FF9999"
scale = 1e6 #for MaF conversion 
names = c("greenriverwy", "greenriverut", "glenwood", "gunnison", "cisco", "leesferry-short", "leesferry-long")
names_sub = c("greenriverwy", "greenriverut", "glenwood", "gunnison", "cisco", "leesferry-short")
colrs = c("#ED9121","#FF6347", "#88DEB0", "#377A98", "#1034A6", "#003152", "#003152")
startY_list = c(1569, 1569, 1569, 1569, 1569, 1569, 762)
upper_sig_list = c(0.77, 0.73, 0.61, 0.61, 0.68, 0.92, 0.95)
#"#2573B5", "#4EADAF"
#make labels for area fraction plots 
llabels = c("(c)", "(e)", "(f)", "(b)", "(d)", "(a)", "(g)")

#high and low predictability years [grabbed from high-low_pred_years.csv]
#by getting the years that fall above / below +0.05 / -0.05 and are not close to each other
low_pred = list(c(1933,1998), #greenriverwy
                c(), #greenriverut
                c(1987,1844,1600, 1653), #glenwood
                c(1981,1688,1854,1821), #gunnison
                c(1983,1853,1732,1684), #cisco
                c(1980,1844,1612,1722), #leesferry-short
                c(991,1325,1990,1847)) #leesferry-long 

high_pred = list(c(1820,1872,1692,1650), #greenriverwy
                 c(1955,1634,1751,1777), #greenriverut
                 c(1630,1756,1966,1821), #glenwood 
                 c(1780,1889,1709), #gunnison 
                 c(1769,1628,1895,1939), #cisco
                 c(1787,1882,1941,1656), #leesferry-short
                 c(1161,861,1021,1256)) #leesferry-long

#range_val_list = c()
#roll_window_list = c(42,39,42,39,37,37,39)
#low_list = as.integer(low_pred[[7]])
#high_list = as.integer(high_pred[[7]])

#order c("greenriverwy", "greenriverut", "glenwood", "gunnison", "cisco", "leesferry-short")
lag_pick = c(3,2,3,3,2,2)
dim_pick = 3
#time_step_yr = 20
#Multiple figures - area fraction, and avgLLE, roll mean, roll variance comparison w/ windows
for(name in names_sub){
  i = which(name==names_sub)
  data_df = read.csv(paste(datadir,name,"_ts.csv", sep = ""), header = TRUE)
  
  low_list = as.integer(low_pred[[i]])
  high_list = as.integer(high_pred[[i]])
  
  len_high = length(high_list)
  len_low = length(low_list)
  
  range_val = 9 #40 year epochs based on the surrounding area
  
  #high predictability initial year in 20 year epoch
  LLE_lowerbound1 = high_list[1]-range_val-1 
  LLE_lowerbound2 = high_list[2]-range_val-1 
  LLE_lowerbound3 = high_list[3]-range_val-1
  LLE_lowerbound4 = high_list[4]-range_val-1
  
  #low predictability initial year in 20 year epoch
  LLE_lowerbound5 = low_list[1]-range_val-1
  LLE_lowerbound6 = low_list[2]-range_val-1
  LLE_lowerbound7 = low_list[3]-range_val-1
  LLE_lowerbound8 = low_list[4]-range_val-1
  
  #figure windows - avg LLE, mean, var 
  #initial year of the centering 40 yr epoch
  indexs = c(LLE_lowerbound1,LLE_lowerbound2, LLE_lowerbound3, LLE_lowerbound4,
             LLE_lowerbound5,LLE_lowerbound6, LLE_lowerbound7,LLE_lowerbound8)
  epoch = 1:8
  #lag_pick = 2
  types = c(1,1,1,1,2,2,2,2) 
  type_names = c("high", "low")
  low_pred_col = "#993300" #this matches peak colors 
  high_pred_col = "#336699" #this matches valley colors 
  high_list = which(types==1) #list to create high pred. phase space 
  low_list = which(types==2) #list to create low pred. phase space
  ind = which(data_df[,yr_col]%in%indexs, arr.ind = T)[rank(indexs)] #rank(preseves the order indexs list)
  #create index matrix to make the window 
  index.mat = matrix(0,8,2)
  index.mat[,1] = ind 
  index.mat[,2] = ind+19
  data_signal = data_df$signal
  #make embedding signal 
  embb = embedd(data_signal,3,lag_pick[i])
  #ticks <- round(seq(round(min(embb[,1]),1), round(max(embb[,1]),1), length.out=4),1)
  
  #individual plots 
  xlb <- "Y[t]"
  ylb <- paste("Y[t+",lag_pick[i],"]",sep="")
  zlb <- paste("Y[t+",lag_pick[i]+lag_pick[i],"]",sep="")
  #possibly need to change these
  ltys <- 1
  lwds <- 1
  
  #labels = c('', '','', '(d)')
  #labelss = c('', '','', '(h)')
  
  
  png(paste(graphdir,names_sub[i],"_high_pred_trajectories.png", sep=""), width = 8, height = 8, units='in', res = 600)
  par(mfrow=c(2,2), oma = c(0.5,0.5,0.5,0.5), mar= c(1,1,1,1))
  for(j in 1:len_high){
    k=high_list[j]
    rnge=index.mat[k,1]:index.mat[k,2]
    my.plot<-scatterplot3d(embb[,1],embb[,2],embb[,3],type="l",xlab=xlb ,ylab=ylb,zlab=zlb,color="gray", 
                           cex.lab = 0.8,cex.axis = 0.8,lwd=2, 
                           axis = TRUE ,tick.marks = TRUE,grid = TRUE, box = FALSE, 
                           main = paste(indexs[k],'-', indexs[k]+19, sep=" "))
    
    # Customize axis labels with fewer ticks
    
    #y_ticks <- c(min(embb[,1]), median(embb[,1]), max(embb[,1]))
    #z_ticks <- c(min(embb[,1]), median(embb[,1]), max(embb[,1]))
    
    #axis(1, at = x_ticks, labels = format(x_ticks, nsmall = 1, scientific = FALSE))
    #axis(2, at = y_ticks, labels = format(y_ticks, nsmall = 1, scientific = FALSE))
    #axis(3, at = z_ticks, labels = format(z_ticks, nsmall = 1, scientific = FALSE))
    
    my.plot$points3d(embb[rnge,1],embb[rnge,2],embb[rnge,3],col=high_pred_col,type="l",lwd=2,lty = ltys)
    my.plot$points3d(embb[(index.mat[k,1]),1],embb[(index.mat[k,1]),2],embb[(index.mat[k,1]),3],col=high_pred_col,type="p",lwd=2,lty = ltys)
    
    #mtext(side=3,paste(labels[j]),adj=0,line = -1.5,cex=1.3)
    #axis(1, )
  }
  #mtext("High Predictability trajectories", side = 1, outer = T, line=0)
  dev.off()
  
  
  if (len_low == 0) {
    print("no low pred figure")
  } else if (len_low>0) {
    png(paste(graphdir,names_sub[i],"_low_pred_trajectories.png", sep=""), width = 8, height = 8, units='in', res = 600)
    par(mfrow=c(2,2), oma = c(0.5,0.5,0.5,0.5), mar= c(1,1,1,1))
    for(j in 1:len_low){
      k=low_list[j]
      rnge=index.mat[k,1]:index.mat[k,2]
      my.plot<-scatterplot3d(embb[,1],embb[,2],embb[,3],type="l",xlab=xlb ,ylab=ylb,zlab=zlb,color="gray",cex.lab = 0.8,cex.axis = 0.8, lwd=2,axis = TRUE, box = FALSE, main = paste(indexs[k],'-', indexs[k]+19, sep=" "))
      my.plot$points3d(embb[rnge,1],embb[rnge,2],embb[rnge,3],col=low_pred_col,type="l",lwd=2,lty = ltys)
      my.plot$points3d(embb[(index.mat[k,1]),1],embb[(index.mat[k,1]),2],embb[(index.mat[k,1]),3],col=low_pred_col,type="p",lwd=2,lty = ltys)
      #mtext(side=3,paste(labelss[j]),adj=0,line = -1.5,cex=1.3)
    }
    dev.off()
    print ("continue")
  }
  
}

