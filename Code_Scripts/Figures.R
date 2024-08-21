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

#####
#Columns in all data sets are 
#df[,1] = year
#df[,2] = flow
#df[,3] = signal 
#df[,4] = avg_LLE
#df[,5] = sawp_avg 
#df[,6] = rolling mean
#df[,7] = rolling variance 

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
colrs = c("#88DEB0", "#ED9121", "#FF6347", "#377A98", "#1034A6", "#003152", "#003152")
startY_list = c(1569, 1569, 1569, 1569, 1569, 1569, 762)
upper_sig_list = c(0.77, 0.73, 0.61, 0.61, 0.68, 0.92, 0.95)

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

#single figures data load 
leesferry_long_df = read.csv("./Output_Data/leesferry-long_ts.csv", header = TRUE)
leesferry_short_df = read.csv("./Output_Data/leesferry-short_ts.csv", header = TRUE)

sf_mean = mean(leesferry_long_df[,flow_col]/scale)
sf_sd = sd(leesferry_long_df[,flow_col]/scale)
#this data signal is expressed in MaF 
data_signal = scale(leesferry_long_df[,signal_col]) * sf_sd + sf_mean

#indexes for non-na values
val_i=which(!is.na(leesferry_long_df$avgLLE))
#grab index of pre, post, and observed years [used for box plots and scatterplots]
#762 - 1489
pre_i <- c(val_i[1]:which(leesferry_long_df$year==1489))
#1490 - 1905
post_i <- c((which(leesferry_long_df$year==1490)):which(leesferry_long_df$year==1905))
#1906 - 2020
obs_i <- c(which(leesferry_long_df$year==1906):tail(val_i, n=1))

#get the years of high and low predictability for leesferry-long only 
low_list = as.integer(low_pred[[7]])
high_list = as.integer(high_pred[[7]])

#######Figure 2 (a)#######
#leesferry-long time series for showing flow only 
png(paste(graphdir,"leesferry-long_flowts_1.png",sep = ""), width = 400, height = 150, units='mm', res = 300)
par(mfrow=c(1,1),mar=c(4,5,3,1), mgp = c(2.5, 1, 0), cex.axis = 2, cex.lab =2)
plot(leesferry_long_df[,yr_col], leesferry_long_df[,flow_col]/scale, ylim = c(0, range(leesferry_long_df[,flow_col]/scale)[2]), 
     type = 'l', lwd = 2.5, 
     col= colrs[6], 
     panel.first = grid(NULL,NULL, lty = 3, col = "grey"),
     xlab = "Years", ylab = "Streamflow (MaF)", cex.lab = 2, 
     cex.axis = 2)
abline(v = leesferry_long_df[which(leesferry_long_df[,yr_col]==1569),yr_col], col = 'darkgrey', lwd = 6, lty = 2)
grid (NULL,NULL, lty = 3, col = "grey")
mtext(side=3,"(a)",adj=0,cex=2)
legend("bottomright", legend="leesferry-long", lwd = 3, lty = 1, col = colrs[6], bty = 'n', cex =2)
dev.off()

#######Figure 2 (b)#######
#all gauge time series for showing flow only
png(paste(graphdir,"allgauge_flowts.png",sep = ""), width = 400, height = 150, units='mm', res = 300)
par(mfrow=c(1,1),mar=c(4,5,3,1), mgp = c(2.5, 1, 0), cex.axis = 2, cex.lab =2)
plot(leesferry_short_df[,yr_col], leesferry_short_df[,flow_col]/scale, ylim = c(0, range(leesferry_short_df[,flow_col]/scale)[2]), 
     type = 'l', lwd = 2.5, col= colrs[6],panel.first = grid(NULL,NULL, lty = 3, col = "grey"),
     xlab = "Years", ylab = "Streamflow (MaF)")
for(i in 1:5){
  data_df = read.csv(paste(datadir,names[i],"_ts.csv", sep = ""), header = TRUE)
  lines(data_df[,yr_col], data_df[,flow_col]/scale, col=colrs[i], lwd = 2.5)
}
mtext(side=3,"(b)",adj=0,cex=2)
par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0), mar=c(0, 0, 0, 0), new=TRUE)
legend("topright",horiz = T, legend=c(names[1],names[2],names[3],names[4], names[5],names[6]), lwd = 3, lty = 1, col = c(colrs[1],colrs[2],colrs[3],colrs[4],colrs[5],colrs[6]), bty = 'n', cex =1.3)
dev.off()

#######Figure 4 (a)#######
png(paste(graphdir,"leesferry-long_flowts.png",sep = ""), width = 400, height = 150, units='mm', res = 300)
par(mfrow=c(1,1),mar=c(4,5,3,1), mgp = c(2.5, 1, 0), cex.axis = 2, cex.lab =2)
plot(leesferry_long_df[,yr_col], leesferry_long_df[,flow_col]/scale, ylim = c(0, range(leesferry_long_df[,flow_col]/scale)[2]), 
     type = 'l', lwd = 2.5, 
     col= 'black', 
     panel.first = grid(NULL,NULL, lty = 3, col = "grey"),
     xlab = "Years", ylab = "Streamflow (MaF)", cex.lab = 2, 
     cex.axis = 2)
grid (NULL,NULL, lty = 3, col = "grey")
mtext(side=3,"(a)",adj=0,cex=2)
dev.off()

#######Figure 4 (b) & (c)#######
#creates wavelet plot 
wplot = wplot(leesferry_long_df[,flow_col]/scale, 762, name = "leesferry-long_PS", graphdir=graphdir)

#######Figure 4 (d)#######
#leesferry-long signal time series for wavelet figure 
png(paste(graphdir,"leesferry-long_signalts.png",sep = ""), width = 400, height = 150, units='mm', res = 300)
par(mfrow=c(1,1),mar=c(4,5,3,1), mgp = c(2.5, 1, 0), cex.axis = 2, cex.lab =2)
plot(leesferry_long_df[,yr_col], data_signal, ylim = c(0, range(data_signal)[2]), 
     type = 'l', lwd = 2.5, col= 'black', panel.first = grid(NULL,NULL, lty = 3, col = "grey"),
     xlab = "Years", ylab = "Signal (MaF)", cex.lab = 2, 
     cex.axis = 2)
mtext(side=3,"(d)",adj=0,cex=2)
dev.off()

#######Figure 4 (e)#######
png(paste(graphdir,"leesferry-long_sawp.png",sep = ""), width = 400, height = 150, units='mm', res = 300)
par(mfrow=c(1,1),mar=c(4,5,3,1), mgp = c(2.5, 1, 0), cex.axis = 2, cex.lab =2)
plot(leesferry_long_df[,yr_col], leesferry_long_df[,sawp_col], 
     type = 'l', lwd = 2.5, col= 'black', panel.first = grid(NULL,NULL, lty = 3, col = "grey"),
     xlab = "Years", ylab = "SAWP", cex.lab = 2, 
     cex.axis = 2)
mtext(side=3,"(e)",adj=0,cex=2)
dev.off()

#######Figure 7 (a),(b),(c)#######
#temp: leesferry-long streaflow, signal, avg LLE 
png("./Figures/leesferry-long_ts.png", width = 9, height = 9, units='in', res = 300)
par(mfrow=c(3,1),mar=c(4,4,2,1), mgp = c(2, .6, 0))
plot(leesferry_long_df[,yr_col], leesferry_long_df[,flow_col]/scale, ylim = c(0, range(leesferry_long_df[,flow_col]/scale)[2]), 
     type = 'l', lwd = 1.5, 
     xlab = "Years", ylab = "Streamflow (MaF)", cex.lab = 1.8, panel.first = grid(NULL,NULL, lty = 3, col = "grey"), 
     cex.axis = 1.8)
abline(v=leesferry_long_df[which(leesferry_long_df[,yr_col]==1906),yr_col], lwd = 4, col= "darkorange3")
mtext(side=3,"(a)",adj=0,cex=1.3)

plot(leesferry_long_df[,yr_col], data_signal, ylim = c(0, range(data_signal)[2]), 
     type = 'l', lwd = 1.5,
     xlab = "Years", ylab = "Signal (MaF)", cex.lab = 1.8, panel.first = grid(NULL,NULL, lty = 3, col = "grey"),
     cex.axis = 1.8)
mtext(side=3,"(b)",adj=0,cex=1.3)

plot(leesferry_long_df[,yr_col], leesferry_long_df[,avg_LLE_col], ylim = c(-0.4,0.3), 
     type = 'l', lwd = 1.5,
     xlab = "Years", ylab = "Average LLE", cex.lab = 1.8, panel.first = grid(NULL,NULL, lty = 3, col = "grey"),
     cex.axis = 1.8)
abline(h = 0.05, lty =2)
abline(h = -0.05)


#plot the peaks and valleys
for(i in 1:4){
  points(leesferry_long_df[(which(leesferry_long_df[,yr_col]==low_list[i])),yr_col], leesferry_long_df[(which(leesferry_long_df[,yr_col]==low_list[i])),avg_LLE_col], col="#993300", pch=15, cex= 2)
  points(leesferry_long_df[(which(leesferry_long_df[,yr_col]==high_list[i])),yr_col], leesferry_long_df[(which(leesferry_long_df[,yr_col]==high_list[i])),avg_LLE_col], col="#336699", pch=15, cex= 2)
}

legend("bottomright", legend=c("<-0.05 (high predictability)", ">0.05 (low predictability)"), 
       col=c("black"), lwd=c( 1, 1), lty = c(1,2), cex=1.5, bty="n", border=F)
legend("bottom", legend=c("high predictability", "low predictability"),
       col=c("#336699","#993300"),cex=1.5, pch=15, pt.cex = 2, bty="n", border=F)
mtext(side=3,"(c)",adj=0,cex=1.3)

dev.off()

#######Figure 12 (a)-(h)#######
#leesferry - long scatterplots 
all_avg<-ggplot(data=NULL, aes(leesferry_long_df$avgLLE[val_i], leesferry_long_df$flow_mean[val_i], col = leesferry_long_df$year[val_i])) + 
  stat_poly_line(se = FALSE, color = 'black') + stat_correlation() + geom_point(size=1) +
  scale_color_gradient(high = "#000033", low = "#CCCCFF") + 
  labs(title = '(a)' , x ="", y = "Flow Average (MAF)") + 
  theme(legend.title =element_text(size = 13), axis.title=element_text(size=15),plot.title = element_text(size=15), legend.text =element_text(size = 13),axis.text.x = element_text(size = 13), axis.text.y = element_text(size = 13)) +
  guides(colour="none")

pre_avg<- ggplot(data=NULL, aes(leesferry_long_df$avgLLE[pre_i], leesferry_long_df$flow_mean[pre_i], col = leesferry_long_df$year[pre_i])) + 
  xlim(range(na.omit(leesferry_long_df$avgLLE))) + ylim(range(na.omit(leesferry_long_df$flow_mean))) + 
  stat_poly_line(se = FALSE, color = 'black') + stat_correlation() + geom_point(size=1) +
  scale_color_gradient(high = "#003333", low = "#CCFFFF") +
  labs(title = '(c)' , x ="", y = "Flow Average (MAF)") + 
  theme(legend.title =element_text(size = 13), axis.title=element_text(size=15),plot.title = element_text(size=15), legend.text =element_text(size = 13),axis.text.x = element_text(size = 13), axis.text.y = element_text(size = 13)) +
  guides(colour="none") 

post_avg<- ggplot(data=NULL, aes(leesferry_long_df$avgLLE[post_i], leesferry_long_df$flow_mean[post_i], col = leesferry_long_df$year[post_i])) + 
  xlim(range(na.omit(leesferry_long_df$avgLLE))) + ylim(range(na.omit(leesferry_long_df$flow_mean))) + 
  stat_poly_line(se = FALSE, color = 'black') + stat_correlation() + geom_point(size=1) +
  scale_color_gradient(high = "#FF3300", low = "#FFCC99") + 
  labs(title = '(e)' , x =" ", y = "Flow Average (MAF)") + 
  theme(legend.title =element_text(size = 13), axis.title=element_text(size=15),plot.title = element_text(size=15), legend.text =element_text(size = 13),axis.text.x = element_text(size = 13), axis.text.y = element_text(size = 13)) +
  guides(colour="none") #turns off spectrum on side 

obs_avg <- ggplot(data=NULL, aes(leesferry_long_df$avgLLE[obs_i], leesferry_long_df$flow_mean[obs_i], col = leesferry_long_df$year[obs_i])) + 
  xlim(range(na.omit(leesferry_long_df$avgLLE))) + ylim(range(na.omit(leesferry_long_df$flow_mean))) + 
  stat_poly_line(se = FALSE, color = 'black') + stat_correlation() + geom_point(size=1) +
  scale_color_gradient(high = "#CC33CC", low = "#FFCCFF") + 
  labs(title = '(g)' , x ="Average LLE", y = "Flow Average (MAF)") +
  theme(legend.title =element_text(size = 13), axis.title=element_text(size=15),plot.title = element_text(size=15), legend.text =element_text(size = 13),axis.text.x = element_text(size = 13), axis.text.y = element_text(size = 13)) +
  guides(colour="none") #turns off spectrum on side 

#var SF and avg LLE 
all_var<-ggplot(data=NULL, aes(leesferry_long_df$avgLLE[val_i], leesferry_long_df$flow_var[val_i], col = leesferry_long_df$year[val_i])) + 
  xlim(range(leesferry_long_df$avgLLE[val_i])) + ylim(c(0,30)) + 
  stat_poly_line(se = FALSE, color = 'black') + stat_correlation() + geom_point(size=1) + 
  scale_color_gradient(high = "#000033", low = "#CCCCFF") + 
  labs(title = '(b)' , x ="", y = "Flow Variance (MAF)", color = "All Years") + 
  theme(legend.title =element_text(size = 13), axis.title=element_text(size=15),plot.title = element_text(size=15), legend.text =element_text(size = 13),axis.text.x = element_text(size = 13), axis.text.y = element_text(size = 13))

pre_var<-ggplot(data=NULL, aes(leesferry_long_df$avgLLE[pre_i], leesferry_long_df$flow_var[pre_i], col = leesferry_long_df$year[pre_i])) + 
  xlim(range(leesferry_long_df$avgLLE[val_i])) + ylim(c(0,30)) + 
  stat_poly_line(se = FALSE, color = 'black') + stat_correlation() + geom_point(size=1) + 
  scale_color_gradient(high = "#003333", low = "#CCFFFF") + 
  labs(title = '(d)' , x ="", y = "Flow Variance (MAF)", color = "762 - 1489") + 
  theme(legend.title =element_text(size = 13), axis.title=element_text(size=15),plot.title = element_text(size=15), legend.text =element_text(size = 13),axis.text.x = element_text(size = 13), axis.text.y = element_text(size = 13))

post_var<-ggplot(data=NULL, aes(leesferry_long_df$avgLLE[post_i], leesferry_long_df$flow_var[post_i], col = leesferry_long_df$year[post_i])) + 
  xlim(range(leesferry_long_df$avgLLE[val_i])) + ylim(c(0,30)) + 
  stat_poly_line(se = FALSE, color = 'black') + stat_correlation() + geom_point(size=1) + 
  scale_color_gradient(high = "#FF3300", low = "#FFCC99") + 
  labs(title = '(f)' , x ="", y = "Flow Variance (MAF)", color = "1490 - 1905") + 
  theme(legend.title =element_text(size = 13), axis.title=element_text(size=15),plot.title = element_text(size=15), legend.text =element_text(size = 13),axis.text.x = element_text(size = 13), axis.text.y = element_text(size = 13))

obs_var<-ggplot(data=NULL, aes(leesferry_long_df$avgLLE[obs_i],leesferry_long_df$flow_var[obs_i], col = leesferry_long_df$year[obs_i])) + 
  xlim(range(leesferry_long_df$avgLLE[val_i])) + ylim(c(0,30)) + 
  stat_poly_line(se = FALSE, color = 'black') + stat_correlation() + geom_point(size=1) + 
  scale_color_gradient(high = "#CC33CC", low = "#FFCCFF") + 
  labs(title = '(h)' , x ="Average LLE", y = "Flow Variance (MAF)", color = "1906 - 2019") + 
  theme(legend.title =element_text(size = 13), axis.title=element_text(size=15),plot.title = element_text(size=15), legend.text =element_text(size = 13),axis.text.x = element_text(size = 13), axis.text.y = element_text(size = 13))

png(paste(graphdir,"leesferry-long_scatterplots.png", sep=""), width = 220, height = 320, units='mm', res = 300)
(all_avg | all_var) / (pre_avg | pre_var) / (post_avg | post_var) / (obs_avg | obs_var)
dev.off()


#######Figure 14 (a)-(h)#######
range_val = 19 #half of 39 from Table 2

#high predictability initial year in 40 year epoch
LLE_lowerbound1 = high_list[1]-range_val-1 
LLE_lowerbound2 = high_list[2]-range_val-1 
LLE_lowerbound3 = high_list[3]-range_val-1
LLE_lowerbound4 = high_list[4]-range_val-1

#low predictability initial year in 40 year epoch
LLE_lowerbound5 = low_list[1]-range_val-1
LLE_lowerbound6 = low_list[2]-range_val-1
LLE_lowerbound7 = low_list[3]-range_val-1
LLE_lowerbound8 = low_list[4]-range_val-1

#figure windows - avg LLE, mean, var 
#initial year of the centering 40 yr epoch
indexs = c(LLE_lowerbound1,LLE_lowerbound2, LLE_lowerbound3, LLE_lowerbound4,
           LLE_lowerbound5,LLE_lowerbound6, LLE_lowerbound7,LLE_lowerbound8)
epoch = 1:8
lag_pick = 2
types = c(1,1,1,1,2,2,2,2) 
type_names = c("high", "low")
low_pred_col = "#993300" #this matches peak colors 
high_pred_col = "#336699" #this matches valley colors 
high_list = which(types==1) #list to create high pred. phase space 
low_list = which(types==2) #list to create low pred. phase space
ind = which(leesferry_long_df[,yr_col]%in%indexs, arr.ind = T)[rank(indexs)] #rank(preseves the order indexs list)
#create index matrix to make the window 
index.mat = matrix(0,8,2)
index.mat[,1] = ind 
index.mat[,2] = ind+39
#make embedding signal 
embb = embedd(data_signal,3,lag_pick)

#individual plots 
xlb <- "Y[t]"
ylb <- paste("Y[t+",lag_pick,"]",sep="")
zlb <- paste("Y[t+",lag_pick+lag_pick,"]",sep="")
#possibly need to change these
ltys <- 1
lwds <- 1

labels = c('(a)', '(b)','(c)', '(d)')
labelss = c('(e)', '(f)','(g)', '(h)')

png(paste(graphdir,"leesferry-long_high_pred_trajectories.png", sep=""), width = 8, height = 8, units='in', res = 300)
par(mfrow=c(2,2), oma = c(0.5,0.5,0.5,0.5))
for(j in 1:4){
  i=high_list[j]
  rnge=index.mat[i,1]:index.mat[i,2]
  my.plot<-scatterplot3d(embb[,1],embb[,2],embb[,3],type="l",xlab=xlb ,ylab=ylb,zlab=zlb,color="gray", cex.lab = 1.1,cex.axis = 1.1,lwd=2,axis = TRUE, box = FALSE, main = paste(indexs[i],'-', indexs[i]+39, sep=" "))
  my.plot$points3d(embb[rnge,1],embb[rnge,2],embb[rnge,3],col=high_pred_col,type="l",lwd=2,lty = ltys)
  my.plot$points3d(embb[(index.mat[i,1]),1],embb[(index.mat[i,1]),2],embb[(index.mat[i,1]),3],col=high_pred_col,type="p",lwd=2,lty = ltys)
  mtext(side=3,paste(labels[j]),adj=0,line = -1.5,cex=1.3)
}
#mtext("High Predictability trajectories", side = 1, outer = T, line=0)
dev.off()

png(paste(graphdir,"leesferry-long_low_pred_trajectories.png", sep=""), width = 8, height = 8, units='in', res = 300)
par(mfrow=c(2,2), oma = c(0.5,0.5,0.5,0.5))
for(j in 1:4){
  i=low_list[j]
  rnge=index.mat[i,1]:index.mat[i,2]
  my.plot<-scatterplot3d(embb[,1],embb[,2],embb[,3],type="l",xlab=xlb ,ylab=ylb,zlab=zlb,color="gray",cex.lab = 1.1,cex.axis = 1.1, lwd=2,axis = TRUE, box = FALSE, main = paste(indexs[i],'-', indexs[i]+39, sep=" "))
  my.plot$points3d(embb[rnge,1],embb[rnge,2],embb[rnge,3],col=low_pred_col,type="l",lwd=2,lty = ltys)
  my.plot$points3d(embb[(index.mat[i,1]),1],embb[(index.mat[i,1]),2],embb[(index.mat[i,1]),3],col=low_pred_col,type="p",lwd=2,lty = ltys)
  mtext(side=3,paste(labelss[j]),adj=0,line = -1.5,cex=1.3)
}
dev.off()

#######Figure 10 (a)&(b)#######
#predictability for all locations 
epoch = 1:12
#order c("greenriverwy", "greenriverut", "glenwood", "gunnison", "cisco", "leesferry-short")
lag_pick = c(3,2,3,3,2,2)
dim_pick = 3
time_step_yr = 20

png(paste(graphdir,"high_pred_trajectories.png", sep=""), width = 300, height = 200, units='mm', res = 300)
par(mfrow=c(2,3), oma = c(6, 6, 4, 4))
for(name in names_sub){
  data_df = read.csv(paste(datadir,name,"_ts.csv", sep = ""), header = TRUE)
  i = which(name==names_sub)
  signal = data_df$signal
  #sf_mean = mean(data_df$flow/scale)
  #sf_sd = sd(data_df$flow/scale)
  #this data signal is expressed in MaF 
  #data_signal = signal * sf_sd + sf_mean
  embb = embedd(signal, dim_pick, lag_pick[i])
  start_yr = which(data_df$year==1889)
  rnge=start_yr:(start_yr+19)
  xlim = range(data_df$signal)
  ylim = range(data_df$signal)
  zlim = range(data_df$signal)
  xlb <- "Y[t]"
  ylb <- paste("Y[t+",lag_pick[i],"]",sep="")
  zlb <- paste("Y[t+",lag_pick[i]+lag_pick[i],"]",sep="")
  my.plot<-scatterplot3d(embb[,1],embb[,2],embb[,3],type="l",xlab=xlb ,ylab=ylb,zlab=zlb,xlim = xlim,ylim = ylim, zlim = zlim, box=F, color="gray",lwd=2, cex.main = 2, main = paste(names_sub[i]))
  my.plot$points3d(embb[rnge,1],embb[rnge,2],embb[rnge,3],col=high_pred_col,type="l",lwd=2,lty = 1)
  my.plot$points3d(embb[(rnge[1]),1],embb[(rnge[1]),2],embb[(rnge[1]),3],col=high_pred_col,type="p",lwd=2,lty = 1)
}
mtext("High Predictability trajectory [1889-1905]", side = 1, outer = T, line=0, cex=1.5)
mtext("(a)", side = 3, outer = T, line = 0, adj = 0, cex=2)
dev.off()

png(paste(graphdir,"low_pred_trajectories.png", sep=""), width = 300, height = 200, units='mm', res = 300)
par(mfrow=c(2,3), oma = c(6, 6, 4, 4))
for(name in names_sub){
  data_df = read.csv(paste(datadir,name,"_ts.csv", sep = ""), header = TRUE)
  i = which(name==names_sub)
  signal = data_df$signal
  #sf_mean = mean(data_df$flow/scale)
  #sf_sd = sd(data_df$flow/scale)
  #this data signal is expressed in MaF 
  #data_signal = signal * sf_sd + sf_mean
  embb = embedd(signal, dim_pick, lag_pick[i])
  start_yr = which(data_df$year==1980)
  rnge=start_yr:(start_yr+19)
  xlim = range(data_df$signal)
  ylim = range(data_df$signal)
  zlim = range(data_df$signal)
  xlb <- "Y[t]"
  ylb <- paste("Y[t+",lag_pick[i],"]",sep="")
  zlb <- paste("Y[t+",lag_pick[i]+lag_pick[i],"]",sep="")
  my.plot<-scatterplot3d(embb[,1],embb[,2],embb[,3],type="l",xlab=xlb ,ylab=ylb,zlab=zlb,color="gray",box=F,lwd=2, cex.main = 2, main = paste(names_sub[i]))
  my.plot$points3d(embb[rnge,1],embb[rnge,2],embb[rnge,3],col=low_pred_col,type="l",lwd=2,lty = 1)
  my.plot$points3d(embb[(rnge[1]),1],embb[(rnge[1]),2],embb[(rnge[1]),3],col=low_pred_col,type="p",lwd=2,lty = 1)
}
mtext("Low Predictability trajectory [1980-1999]", side = 1, outer = T, line=0, cex=1.5)
mtext("(b)", side = 3, outer = T, line = 0, adj = 0, cex=2)
dev.off()


#######Figure 13 (a)&(b)#######
#boxplots

#grab indices of where avg LLE falls below -0.05 and above 0.05 
high_pred_i <- which(leesferry_long_df[,avg_LLE_col] < -0.05)
low_pred_i <- which(leesferry_long_df[,avg_LLE_col] > 0.05)

#break each section of high predictability into pre, post, obs years 
#these are indexes only 
high_pre_i <- high_pred_i[high_pred_i %in% pre_i]
high_post_i <- high_pred_i[high_pred_i %in% post_i]
high_obs_i <- high_pred_i[high_pred_i %in% obs_i]

#break each section of low predictability into pre, post, obs years 
#these are indexes only 
low_pre_i <- low_pred_i[low_pred_i %in% pre_i]
low_post_i <- low_pred_i[low_pred_i %in% post_i]
low_obs_i <- low_pred_i[low_pred_i %in% obs_i]

#mean 
high_pred_avg <- c(leesferry_long_df$flow_mean[high_pred_i],matrix(NA, nrow=(length(high_pred_i)), ncol=1))
high_pred_avg_pre <- c(leesferry_long_df$flow_mean[high_pre_i],matrix(NA, nrow=(length(high_pre_i)), ncol=1))
high_pred_avg_post <- c(leesferry_long_df$flow_mean[high_post_i],matrix(NA, nrow=(length(high_post_i)), ncol=1))
high_pred_avg_obs <- c(leesferry_long_df$flow_mean[high_obs_i],matrix(NA, nrow=(length(high_obs_i)), ncol=1))

low_pred_avg <- c(leesferry_long_df$flow_mean[low_pred_i],matrix(NA, nrow=(length(low_pred_i)), ncol=1))
low_pred_avg_pre <- c(leesferry_long_df$flow_mean[low_pre_i],matrix(NA, nrow=(length(low_pre_i)), ncol=1))
low_pred_avg_post <- c(leesferry_long_df$flow_mean[low_post_i],matrix(NA, nrow=(length(low_post_i)), ncol=1))
low_pred_avg_obs <- c(leesferry_long_df$flow_mean[low_obs_i],matrix(NA, nrow=(length(low_obs_i)), ncol=1))

#variance
high_pred_var <- c(leesferry_long_df$flow_var[high_pred_i],matrix(NA, nrow=(length(high_pred_i)), ncol=1))
high_pred_var_pre <- c(leesferry_long_df$flow_var[high_pre_i],matrix(NA, nrow=(length(high_pre_i)), ncol=1))
high_pred_var_post <- c(leesferry_long_df$flow_var[high_post_i],matrix(NA, nrow=(length(high_post_i)), ncol=1))
high_pred_var_obs <- c(leesferry_long_df$flow_var[high_obs_i],matrix(NA, nrow=(length(high_obs_i)), ncol=1))

low_pred_var <- c(leesferry_long_df$flow_var[low_pred_i],matrix(NA, nrow=(length(low_pred_i)), ncol=1))
low_pred_var_pre <- c(leesferry_long_df$flow_var[low_pre_i],matrix(NA, nrow=(length(low_pre_i)), ncol=1))
low_pred_var_post <- c(leesferry_long_df$flow_var[low_post_i],matrix(NA, nrow=(length(low_post_i)), ncol=1))
low_pred_var_obs <- c(leesferry_long_df$flow_var[low_obs_i],matrix(NA, nrow=(length(low_obs_i)), ncol=1))

png(paste(graphdir,"pred_var_boxplot.png", sep=""), width = 6, height = 6, units='in', res = 300)
par(mfcol=c(1,4), mar = numeric(4), oma = c(6, 6, 4, 4),  mgp = c(2, .6, 0))
boxplot(high_pred_var,low_pred_var, axes = FALSE,
        at = c(1,2), ylim=c(0,25), border = "#000033", col = "#CCCCFF") 
axis(2, cex.axis = 1.8)
mtext("Streamflow Variance (MAF)", side = 2, outer = T, line = 2, adj = 0.5, cex = 1.2)
mtext("(b)       All", side = 3, outer = T, line = 0.5, adj = 0, cex=1.2)
mtext("High", side = 1, outer = T, line = 0.5, adj = 0.02, cex=1.2)
mtext("Low", side = 1, outer = T, line = 0.5, adj = 0.16, cex=1.2)
box()

boxplot(high_pred_var_pre,low_pred_var_pre, axes = FALSE,
        at = c(1,2), ylim=c(0,25), border = "#003333", col = "#CCFFFF")
mtext("762-1489", side = 3, outer = T, line = 0.5, adj = 0.35, cex = 1.2)
mtext("High", side = 1, outer = T, line = 0.5, adj = 0.3, cex=1.2)
mtext("Low", side = 1, outer = T, line = 0.5, adj = 0.43, cex=1.2)
box()

boxplot(high_pred_var_post,low_pred_var_post, axes = FALSE,
        at = c(1,2), ylim=c(0,25), border = "#FF3300", col = "#FFCC99")
mtext("1490-1905", side = 3, outer = T, line = 0.5, adj = 0.65, cex=1.2)
mtext("High", side = 1, outer = T, line = 0.5, adj = 0.57, cex=1.2)
mtext("Low", side = 1, outer = T, line = 0.5, adj = 0.71, cex=1.2)
box()

boxplot(high_pred_var_obs,low_pred_var_obs, axes = FALSE,
        at = c(1,2), ylim=c(0,25), border = "#CC33CC", col = "#FFCCFF")
mtext("1906-2019", side = 3, outer = T, line = 0.5, adj = 0.98, cex=1.2)
mtext("High", side = 1, outer = T, line = 0.5, adj = 0.85, cex=1.2)
mtext("Low", side = 1, outer = T, line = 0.5, adj = 0.98, cex=1.2)
mtext("Predictability", side = 1, outer = T, line = 2.5, adj = 0.5, cex=1.2)
box()
dev.off()

png(paste(graphdir,"pred_avg_boxplot.png", sep=""), width = 6, height = 6, units='in', res = 300)
par(mfcol=c(1,4), mar = numeric(4), oma = c(6, 6, 4, 4),  mgp = c(2, .6, 0))
boxplot(high_pred_avg,low_pred_avg, axes = FALSE,
        at = c(1,2), ylim=c(12,18), border = "#000033", col = "#CCCCFF")
axis(2, cex.axis=1.8)
mtext("Streamflow Average (MAF)", side = 2, outer = T, line = 2, adj = 0.5, cex = 1.2)
mtext("(a)       All", side = 3, outer = T, line = 0.5, adj = 0, cex=1.2)
mtext("High", side = 1, outer = T, line = 0.5, adj = 0.02, cex=1.2)
mtext("Low", side = 1, outer = T, line = 0.5, adj = 0.16, cex=1.2)
box()

boxplot(high_pred_avg_pre,low_pred_avg_pre, axes = FALSE,
        at = c(1,2), ylim=c(12,18), border = "#003333", col = "#CCFFFF")
mtext("762-1489", side = 3, outer = T, line = 0.5, adj = 0.35, cex = 1.2)
mtext("High", side = 1, outer = T, line = 0.5, adj = 0.3, cex=1.2)
mtext("Low", side = 1, outer = T, line = 0.5, adj = 0.43, cex=1.2)
box()

boxplot(high_pred_avg_post,low_pred_avg_post, axes = FALSE,
        at = c(1,2), ylim=c(12,18), border = "#FF3300", col = "#FFCC99")
mtext("1490-1905", side = 3, outer = T, line = 0.5, adj = 0.65, cex=1.2)
mtext("High", side = 1, outer = T, line = 0.5, adj = 0.57, cex=1.2)
mtext("Low", side = 1, outer = T, line = 0.5, adj = 0.71, cex=1.2)
box()

boxplot(high_pred_avg_obs,low_pred_avg_obs, axes = FALSE,
        at = c(1,2), ylim=c(12,18), border = "#CC33CC", col = "#FFCCFF")
mtext("1906-2019", side = 3, outer = T, line = 0.5, adj = 0.98, cex=1.2)
mtext("High", side = 1, outer = T, line = 0.5, adj = 0.85, cex=1.2)
mtext("Low", side = 1, outer = T, line = 0.5, adj = 0.98, cex=1.2)
mtext("Predictability", side = 1, outer = T, line = 2.5, adj = 0.5, cex=1.2)
box()
dev.off()

#######Figure 3, Figure 6, Figure 9#######
#correlation plots of streamflow, signal, and average LLE
#Loading all gauge data and storing it in a dataframe without opening all
cols = c(flow_col, signal_col, avg_LLE_col)
corr_name = c('flow', 'signal', 'avgLLE')
plot_limit = c(0,0,-.1)
for(i in 1:3){
  col_num = cols[i]
  Multi_df = data.frame(matrix(nrow = 451,ncol = 0))
  for(name in names_sub){
    #i = which(name==names)
    data_df = read.csv(paste(datadir,name,"_ts.csv", sep = ""), header = TRUE)[col_num]
    Multi_df[,ncol(Multi_df)+1] <- c(data_df)
  }
  colnames(Multi_df) = names_sub
  
  png(paste(graphdir,corr_name[i],"_corr_plot.png", sep=""), width = 150, height = 150, units='mm', res = 300)
  Multp_plot = corrplot(cor(na.omit(Multi_df)), method = 'number',tl.col = 'black',is.corr = T, col.lim = c(plot_limit[i],1), tl.srt = 45, type = 'lower', col= COL2("BrBG"))
  dev.off()
}
#######Figure S2#######
#scatterplots for all locations 
for(name in names_sub){
  data_df = read.csv(paste(datadir,name,"_ts.csv", sep = ""), header = TRUE)
  val_i=which(!is.na(data_df$avgLLE))
  
  all_avg<-ggplot(data=NULL, aes(data_df$avgLLE[val_i], data_df$flow_mean[val_i], col = data_df$year[val_i])) + 
    xlim(range(data_df$avgLLE[val_i])) + ylim(range(data_df$flow_mean[val_i])) + 
    stat_poly_line(se = FALSE, color = 'black') + stat_correlation() + geom_point(size=1) +
    scale_color_gradient(high = "#000033", low = "#CCCCFF") + 
    labs(title = paste('(a)',name, sep = " ") , x ="Average LLE", y = "Flow Average (MAF)") + 
    theme(legend.title =element_text(size = 13), axis.title=element_text(size=15),plot.title = element_text(size=15), legend.text =element_text(size = 13),axis.text.x = element_text(size = 13), axis.text.y = element_text(size = 13)) +
    guides(colour="none")
  
  all_var<-ggplot(data=NULL, aes(data_df$avgLLE[val_i], data_df$flow_var[val_i], col = data_df$year[val_i])) + 
    xlim(range(data_df$avgLLE[val_i])) + ylim(range(data_df$flow_var[val_i])) + 
    stat_poly_line(se = FALSE, color = 'black') + stat_correlation() + geom_point(size=1) + 
    scale_color_gradient(high = "#000033", low = "#CCCCFF") + 
    labs(title = '(b)' , x ="Average LLE", y = "Flow Variance (MAF)", color = "1569 - 2019") + 
    theme(legend.title =element_text(size = 13), axis.title=element_text(size=15),plot.title = element_text(size=15), legend.text =element_text(size = 13),axis.text.x = element_text(size = 13), axis.text.y = element_text(size = 13))
  
  png(paste(graphdir, name, "_scatterplots.png", sep=""), width = 250, height = 100, units='mm', res = 300)
    show(all_avg | all_var)
  dev.off()
}

#######Figure 8#######
#predictability windows only for each gauge
png(paste(graphdir,"pred_windows.png", sep=""), width = 200, height = 200, units='mm', res = 300)
par(mfcol=c(6,1), mar=numeric(4), oma = c(6, 6, 1, 4),  mgp = c(2, .6, 0)) #oma(bottom,left,top,right)
for(name in names){
  i = which(name==names)
  data_df = read.csv(paste(datadir,name,"_ts.csv", sep = ""), header = TRUE)
  #graph y-axis ranges 
  avglle_l = -0.3
  avglee_u = 0.7
  range_val = 9 
  if(i<=6){
    low_list = as.integer(low_pred[[i]])
    high_list = as.integer(high_pred[[i]])
    
    #high predictability 20 year epochs
    LLE_lowerbound1 = high_list[1]-range_val-1; LLE_upperbound1 = high_list[1]+range_val
    LLE_lowerbound2 = high_list[2]-range_val-1; LLE_upperbound2 = high_list[2]+range_val
    LLE_lowerbound3 = high_list[3]-range_val-1; LLE_upperbound3 = high_list[3]+range_val
    LLE_lowerbound4 = high_list[4]-range_val-1; LLE_upperbound4 = high_list[4]+range_val
    
    #low predictability 20 year epochs
    LLE_lowerbound5 = low_list[1]-range_val-1; LLE_upperbound5 = low_list[1]+range_val
    LLE_lowerbound6 = low_list[2]-range_val-1; LLE_upperbound6 = low_list[2]+range_val
    LLE_lowerbound7 = low_list[3]-range_val-1; LLE_upperbound7 = low_list[3]+range_val
    LLE_lowerbound8 = low_list[4]-range_val-1; LLE_upperbound8 = low_list[4]+range_val
    
    plot(data_df[,yr_col], data_df[,avg_LLE_col], type = 'l', 
         panel.first = rect(c(LLE_lowerbound5, LLE_lowerbound6, LLE_lowerbound7, LLE_lowerbound8), -2, c(LLE_upperbound5, LLE_upperbound6, LLE_upperbound7, LLE_upperbound8), 32, col=low_pred_col, border=NA),
         ylim = c(avglle_l,avglee_u),
         axes = FALSE, lwd = 2)
    lines(data_df[,yr_col], data_df[,avg_LLE_col],col='black', lwd=2, 
          panel.first = rect(c(LLE_lowerbound1, LLE_lowerbound2, LLE_lowerbound3, LLE_lowerbound4), -2, c(LLE_upperbound1, LLE_upperbound2, LLE_upperbound3, LLE_upperbound4), 32, col=high_pred_col, border=NA))
    grid (NULL,NULL, lty = 6, col = "grey") 
    abline(h=-0.05,lty=2)
    abline(h=0.05,lty=2)
    axis(2, tck = -0.01, cex.axis = 1.5, las = 2)
    mtext(paste(name), side = 4, line=1)
    box()
    
  } else{
    axis(1,cex.axis = 1.8)
    mtext(side=1, "Years", line = 3, cex = 1.3)
    mtext(side=2, "Average LLE", line = 3, cex = 1.3, outer=T)
    par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
    
    legend("bottomleft", inset = c(0,0), legend=c("low predictability window", "high preditability window"), cex=1.7,
           col=c(low_pred_col, high_pred_col), pch=15, pt.cex = 3, bty="n", border=F)
    dev.off()
    print("This figure is complete")
  }
}

#######Figure 5(a)-(f), Figure 11, and Figure S1 #######
#Multiple figures - area fraction, and avgLLE, roll mean, roll variance comparison w/ windows
for(name in names){
  i = which(name==names)
  data_df = read.csv(paste(datadir,name,"_ts.csv", sep = ""), header = TRUE)
  
  #area fraction plots 
  areafrac = wplot_area(data_df[,flow_col],startY_list[i],dj=0.025,u.sig=upper_sig_list[i],name=name, graphdir, letter_label = llabels[i])
  
  #graph y-axis ranges 
  avglle_l = round((range(na.omit(data_df[,avg_LLE_col]))[1]-0.08), digit = 1)
  avglee_u = round((range(na.omit(data_df[,avg_LLE_col]))[2]+0.08), digit = 1)
  mean_l = round((range(na.omit(data_df[,roll_mean_col]))[1] - 0.2), digit = 1)
  mean_u = round((range(na.omit(data_df[,roll_mean_col]))[2] + 0.2), digit = 1)
  var_l = round((range(na.omit(data_df[,roll_var_col]))[1] - 0.2), digit = 1)
  var_u = round((range(na.omit(data_df[,roll_var_col]))[2] + 0.2), digit = 1)
  
  if(i <= 6){
    #add a window centered from year +- 9 years for the first gauge locations
    range_val = 9 
    print(range_val)
    
  } else{
    #add a window centered from year +- 19 years for leesferry-long
    range_val = 19 
    print(range_val)
  }
  
  low_list = as.integer(low_pred[[i]])
  high_list = as.integer(high_pred[[i]])
  
  #high predictability 20(if i<=6) or 40(if i=7) year epochs
  LLE_lowerbound1 = high_list[1]-range_val-1; LLE_upperbound1 = high_list[1]+range_val
  LLE_lowerbound2 = high_list[2]-range_val-1; LLE_upperbound2 = high_list[2]+range_val
  LLE_lowerbound3 = high_list[3]-range_val-1; LLE_upperbound3 = high_list[3]+range_val
  LLE_lowerbound4 = high_list[4]-range_val-1; LLE_upperbound4 = high_list[4]+range_val
  
  #low predictability 20(if i<=6) or 40(if i=7) year epochs
  LLE_lowerbound5 = low_list[1]-range_val-1; LLE_upperbound5 = low_list[1]+range_val
  LLE_lowerbound6 = low_list[2]-range_val-1; LLE_upperbound6 = low_list[2]+range_val
  LLE_lowerbound7 = low_list[3]-range_val-1; LLE_upperbound7 = low_list[3]+range_val
  LLE_lowerbound8 = low_list[4]-range_val-1; LLE_upperbound8 = low_list[4]+range_val
  
  
  #figure with avgLLE, rolling average and rolling variance with windows
  png(paste(graphdir, name, "_ts_comp.png", sep=""), width = 9, height = 9, units='in', res = 300)
  par(mfcol=c(3,1), mar=numeric(4), oma = c(6, 6, 4, 4),  mgp = c(2, .6, 0)) #oma(bottom,left,top,right)
  plot(data_df[,yr_col], data_df[,avg_LLE_col], type = 'l', 
       panel.first = rect(c(LLE_lowerbound5, LLE_lowerbound6, LLE_lowerbound7, LLE_lowerbound8), -2, c(LLE_upperbound5, LLE_upperbound6, LLE_upperbound7, LLE_upperbound8), 32, col=low_pred_col, border=NA),
       ylim = c(avglle_l,avglee_u),
       axes = FALSE, lwd = 2)
  lines(data_df[,yr_col], data_df[,avg_LLE_col],col='black', lwd=2, 
        panel.first = rect(c(LLE_lowerbound1, LLE_lowerbound2, LLE_lowerbound3, LLE_lowerbound4), -2, c(LLE_upperbound1, LLE_upperbound2, LLE_upperbound3, LLE_upperbound4), 32, col=high_pred_col, border=NA))
  grid (NULL,NULL, lty = 6, col = "grey") 
  abline(h=-0.05,lty=2)
  abline(h=0.05,lty=2)
  axis(2, tck = -0.01, cex.axis = 1.5, las = 2)
  mtext(side=3, paste(name), line = 2, cex = 1.3)
  mtext(side=2, "Average LLE", line = 3, cex = 1.3)
  mtext(side=4, "high pred.",adj = 0,line =1, cex = 1.1)
  mtext(side=4, "low pred.",adj = 1,line =1, cex = 1.1)
  mtext(side=3,"(a)",adj=0,line = -1.5,cex=1.3)
  box()
  
  plot(data_df[,yr_col], data_df[,roll_mean_col], type = 'l', 
       panel.first = rect(c(LLE_lowerbound5, LLE_lowerbound6, LLE_lowerbound7, LLE_lowerbound8), -2, c(LLE_upperbound5, LLE_upperbound6, LLE_upperbound7, LLE_upperbound8), 32, col=low_pred_col, border=NA),
       ylim = c(mean_l,mean_u),
       axes = FALSE, lwd = 2)
  lines(data_df[,yr_col], data_df[,roll_mean_col],col='black', lwd=2, 
        panel.first = rect(c(LLE_lowerbound1, LLE_lowerbound2, LLE_lowerbound3, LLE_lowerbound4), -2, c(LLE_upperbound1, LLE_upperbound2, LLE_upperbound3, LLE_upperbound4), 32, col=high_pred_col, border=NA))
  grid (NULL,NULL, lty = 6, col = "grey") 
  axis(2, cex.axis = 1.5, tck = -0.01, las =2)
  mtext(side=2, "Flow Average (MaF)", line = 3, cex = 1.3)
  mtext(side=3,"(b)",adj=0,line = -1.5,cex=1.3)
  box()
  
  plot(data_df[,yr_col], data_df[,roll_var_col], type = 'l', 
       ylim = c(var_l,var_u),
       panel.first = rect(c(LLE_lowerbound5, LLE_lowerbound6, LLE_lowerbound7, LLE_lowerbound8), -2, c(LLE_upperbound5, LLE_upperbound6, LLE_upperbound7, LLE_upperbound8), 32, col=low_pred_col, border=NA),
       axes = FALSE, lwd = 2)
  lines(data_df[,yr_col], data_df[,roll_var_col],col='black', lwd=2, 
        panel.first = rect(c(LLE_lowerbound1, LLE_lowerbound2, LLE_lowerbound3, LLE_lowerbound4), -2, c(LLE_upperbound1, LLE_upperbound2, LLE_upperbound3, LLE_upperbound4), 32, col=high_pred_col, border=NA))
  grid (NULL,NULL, lty = 6, col = "grey") 
  axis(2, cex.axis = 1.5, tck = -0.01, las =2)
  axis(1, cex.axis = 1.5, tck = -0.01, las =1)
  mtext(side=2, "Flow Variance (MaF)", line = 3, cex = 1.3)
  mtext(side=1, "Years", line = 3, cex = 1.3)
  mtext(side=3,"(c)",adj=0, line = -1.5, cex=1.3)
  box()
  
  par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
  legend("bottomleft", inset = c(0,0), legend=c("low predictability window", "high preditability window"), cex=1.7,
         col=c(low_pred_col, high_pred_col), pch=15, pt.cex = 3, bty="n", border=F)
  
  
  dev.off()
}
