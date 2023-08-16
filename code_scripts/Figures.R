rm(list=ls())
source("./code_scripts/function_library.r")
library(quantmod) #get peaks and valleys of time series 
library(scatterplot3d) #needed for phase space plots 

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
#df[,5] = rolling mean
#df[,6] = rolling variance 

#specify columns 
yr_col = 1 
flow_col = 2
signal_col = 3
avg_LLE_col = 4
roll_mean_col = 5
roll_var_col = 6

#parameters 
high_pred_col = "#99CCFF"
low_pred_col = "#FF9999"
scale = 1e6 #for MaF conversion 
names = c("greenriverwy", "greenriverut", "glenwood", "gunnison", "cisco", "leesferry-short", "leesferry-long")
colrs = c("#88DEB0", "#2573B5", "#4EADAF", "#377A98", "#1034A6", "#003152", "#003152")
startY_list = c(1569, 1569, 1569, 1569, 1569, 1569, 762)
upper_sig_list = c(0.77, 0.73, 0.61, 0.61, 0.68, 0.92, 0.95)

#####
#single figures 
leesferry_long_df = read.csv("./Output_Data/leesferry-long_ts.csv", header = TRUE)
leesferry_short_df = read.csv("./Output_Data/leesferry-short_ts.csv", header = TRUE)

#Figure 1 (a) [once the paper is ready for submission, this is how it should be labeled]
#temp: leesferry-long time series for showing flow only 
#uncomment comments to label the first figure 
png(paste(graphdir,"leesferry-long_flowts.png",sep = ""), width = 400, height = 150, units='mm', res = 300)
par(mfrow=c(1,1),mar=c(4,5,3,1), mgp = c(2.5, 1, 0), cex.axis = 2, cex.lab =2)
plot(leesferry_long_df[,yr_col], leesferry_long_df[,flow_col]/scale, ylim = c(0, range(leesferry_long_df[,flow_col]/scale)[2]), 
     type = 'l', lwd = 2.5, 
     #col= colrs[6], 
     col = 'black',
     panel.first = grid(NULL,NULL, lty = 3, col = "grey"),
     xlab = "Years", ylab = "Streamflow (MaF)", cex.lab = 2, 
     cex.axis = 2)
#abline(v = leesferry_long_df[which(leesferry_long_df[,yr_col]==1569),yr_col], col = 'darkgrey', lwd = 6, lty = 2)
#grid (NULL,NULL, lty = 3, col = "grey")
mtext(side=3,"(a)",adj=0,cex=2)
#legend("bottomright", legend="leesferry-long", lwd = 3, lty = 1, col = colrs[6], bty = 'n', cex =2)
dev.off()

#temp: all time series for showing flow only
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

#creates wavelet plot 
wplot = wplot(leesferry_long_df[,flow_col]/scale, 762, name = "leesferry-long_PS", graphdir=graphdir)

#temp: leesferry-long signal time series for wavelet figure 
sf_mean = mean(leesferry_long_df[,flow_col]/scale)
sf_sd = sd(leesferry_long_df[,flow_col]/scale)
#this data signal is expressed in MaF 
data_signal = scale(leesferry_long_df[,signal_col]) * sf_sd + sf_mean

png(paste(graphdir,"leesferry-long_signalts.png",sep = ""), width = 400, height = 150, units='mm', res = 300)
par(mfrow=c(1,1),mar=c(4,5,3,1), mgp = c(2.5, 1, 0), cex.axis = 2, cex.lab =2)
plot(leesferry_long_df[,yr_col], data_signal, ylim = c(0, range(data_signal)[2]), 
     type = 'l', lwd = 2.5, col= 'black', panel.first = grid(NULL,NULL, lty = 3, col = "grey"),
     xlab = "Years", ylab = "Signal (MaF)", cex.lab = 2, 
     cex.axis = 2)
mtext(side=3,"(d)",adj=0,cex=2)
dev.off()



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

#finds peak indexes 
p_i=findPeaks(leesferry_long_df[,avg_LLE_col]) 
#get the values of LLE for peak indexes 
peak_values=leesferry_long_df[p_i,avg_LLE_col] 
#sort these values in decreasing order (+ lle to - lle)
sort_pval = sort(peak_values, index.return = TRUE, decreasing = TRUE)
#create peak matrix to store year and avg LLE value
peak_mat = matrix(0,length(sort_pval$x),2)
peak_mat[,1] = leesferry_long_df[p_i[sort_pval$ix], yr_col]
peak_mat[,2] = sort_pval$x

#finds valley indexes 
v_i = findValleys(leesferry_long_df[,avg_LLE_col])
#get the values of LLE for valley indexes 
valley_values = leesferry_long_df[v_i,avg_LLE_col] 
#sort these values in decreasing order (+ lle to - lle)
sort_vval = sort(valley_values, index.return = TRUE, decreasing = FALSE)
#create peak matrix to store year and avg LLE value
valley_mat = matrix(0,length(sort_vval$x),2)
valley_mat[,1] = leesferry_long_df[v_i[sort_vval$ix], yr_col]
valley_mat[,2] = sort_vval$x

#I think we should do this ourselves and change it for every case
#looking at the matrix created above
peaks = list()
valleys = list()
#this works only for leesferry-long 
#check that the peaks and valleys are not next to each other (alteast 19 years apart)
for(i in 1:12){
  if(i==1 & abs(v_i[sort_vval$ix[i]] - v_i[sort_vval$ix[i+1]]) < 19){
    valleys <- append(valleys, i)
  }
  else if(abs(v_i[sort_vval$ix[i]] - v_i[sort_vval$ix[i+1]]) > 19){
    valleys <- append(valleys, i)
  }
}
for(i in 1:12){
  if(abs(p_i[sort_pval$ix[i]] - p_i[sort_pval$ix[i+1]]) > 19){
    peaks <- append(peaks, i)
  }
}
valleys <- as.numeric(valleys)
peaks <- as.numeric(peaks)
#plot the peaks and valleys
for(i in 1:4){
  points(leesferry_long_df[p_i[sort_pval$ix[peaks[i]]],yr_col], leesferry_long_df[p_i[sort_pval$ix[peaks[i]]],avg_LLE_col], col="#993300", pch=15, cex= 2)
  points(leesferry_long_df[v_i[sort_vval$ix[valleys[i]]],yr_col], leesferry_long_df[v_i[sort_vval$ix[valleys[i]]],avg_LLE_col], col="#336699", pch=15, cex= 2)
}

legend("bottomright", legend=c("<-0.05 (high predictability)", ">0.05 (low predictability)"), 
       col=c("black"), lwd=c( 1, 1), lty = c(1,2), cex=1.5, bty="n", border=F)
#    col=c("black","#999966", "#999999"), lwd=c(2, 1, 1), lty = c(1,1,2), cex=1.0, bty="n", border=F)
legend("bottom", legend=c("high predictability", "low predictability"),
       col=c("#336699","#993300"),cex=1.5, pch=15, pt.cex = 2, bty="n", border=F)
mtext(side=3,"(c)",adj=0,cex=1.3)

dev.off()


#temp: leesferry-long avg LLE, roll average, roll variance
#ranges for y axis
avglle_l = round((range(na.omit(leesferry_long_df[,avg_LLE_col]))[1]-0.08), digit = 1)
avglee_u = round((range(na.omit(leesferry_long_df[,avg_LLE_col]))[2]+0.08), digit = 1)
mean_l = round((range(na.omit(leesferry_long_df[,roll_mean_col]))[1] - 0.2), digit = 1)
mean_u = round((range(na.omit(leesferry_long_df[,roll_mean_col]))[2] + 0.2), digit = 1)
var_l = round((range(na.omit(leesferry_long_df[,roll_var_col]))[1] - 0.2), digit = 1)
var_u = round((range(na.omit(leesferry_long_df[,roll_var_col]))[2] + 0.2), digit = 1)
high_pred_col = "#99CCFF"
low_pred_col = "#FF9999"
range_val = 19
#range_val = 19 #+ - for window of center avgLLE, var, avg 

#these are in years, not indexes
#high predictability 40 year epochs
LLE_lowerbound1 = leesferry_long_df[v_i[sort_vval$ix[valleys[2]]],yr_col]-range_val-1; LLE_upperbound1 = leesferry_long_df[v_i[sort_vval$ix[valleys[2]]],yr_col]+range_val
LLE_lowerbound2 = leesferry_long_df[v_i[sort_vval$ix[valleys[1]]],yr_col]-range_val-1; LLE_upperbound2 = leesferry_long_df[v_i[sort_vval$ix[valleys[1]]],yr_col]+range_val
LLE_lowerbound3 = leesferry_long_df[v_i[sort_vval$ix[valleys[3]]],yr_col]-range_val-1; LLE_upperbound3 = leesferry_long_df[v_i[sort_vval$ix[valleys[3]]],yr_col]+range_val
LLE_lowerbound4 = leesferry_long_df[v_i[sort_vval$ix[valleys[4]]],yr_col]-range_val-1; LLE_upperbound4 = leesferry_long_df[v_i[sort_vval$ix[valleys[4]]],yr_col]+range_val

#low predictability 40 year epochs
LLE_lowerbound5 = leesferry_long_df[p_i[sort_pval$ix[peaks[1]]],yr_col]-range_val-1; LLE_upperbound5 = leesferry_long_df[p_i[sort_pval$ix[peaks[1]]],yr_col]+range_val
LLE_lowerbound6 = leesferry_long_df[p_i[sort_pval$ix[peaks[3]]],yr_col]-range_val-1; LLE_upperbound6 = leesferry_long_df[p_i[sort_pval$ix[peaks[3]]],yr_col]+range_val
LLE_lowerbound7 = leesferry_long_df[p_i[sort_pval$ix[peaks[2]]],yr_col]-range_val-1; LLE_upperbound7 = leesferry_long_df[p_i[sort_pval$ix[peaks[2]]],yr_col]+range_val
LLE_lowerbound8 = leesferry_long_df[p_i[sort_pval$ix[peaks[4]]],yr_col]-range_val-1; LLE_upperbound8 = leesferry_long_df[p_i[sort_pval$ix[peaks[4]]],yr_col]+range_val

#figure windows - avg LLE, mean, var 
#initial year of the centering 40 yr epoch
indexs = c(LLE_lowerbound1,LLE_lowerbound2, LLE_lowerbound3, LLE_lowerbound4,
           LLE_lowerbound5,LLE_lowerbound6, LLE_lowerbound7,LLE_lowerbound8)
epoch = 1:8
lag_pick = 2
types = c(1,1,1,1,2,2,2,2) 
type_names = c("high", "low")
high_pred_col = "#457B9D"#presentation color  "#99CCFF"#figure color 
low_pred_col = "#FF9999" #"#E63946"presentation color "#FF9999"#figure color 
high_list = which(types==1) #list to create high pred. phase space 
low_list = which(types==2) #list to create low pred. phase space
ind = which(leesferry_long_df[,yr_col]%in%indexs, arr.ind = T)[rank(indexs)] #rank(preseves the order indexs list)
#create index matrix to make the window 
index.mat = matrix(0,8,2)
index.mat[,1] = ind 
index.mat[,2] = ind+38
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

png(paste(graphdir,"high_pred_trajectories.png", sep=""), width = 8, height = 8, units='in', res = 300)
par(mfrow=c(2,2), oma = c(0.5,0.5,0.5,0.5))
for(j in 1:4){
  i=high_list[j]
  rnge=index.mat[i,1]:index.mat[i,2]
  my.plot<-scatterplot3d(embb[,1],embb[,2],embb[,3],type="l",xlab=xlb ,ylab=ylb,zlab=zlb,color="gray",lwd=2,axis = TRUE, box = FALSE, main = paste(indexs[i],'-', indexs[i]+38, sep=" "))
  my.plot$points3d(embb[rnge,1],embb[rnge,2],embb[rnge,3],col=high_pred_col,type="l",lwd=2,lty = ltys)
  my.plot$points3d(embb[(index.mat[i,1]),1],embb[(index.mat[i,1]),2],embb[(index.mat[i,1]),3],col=high_pred_col,type="p",lwd=2,lty = ltys)
  mtext(side=3,paste(labels[j]),adj=0,line = -1.5,cex=1.3)
}
#mtext("High Predictability trajectories", side = 1, outer = T, line=0)
dev.off()

png(paste(graphdir,"low_pred_trajectories.png", sep=""), width = 8, height = 8, units='in', res = 300)
par(mfrow=c(2,2), oma = c(0.5,0.5,0.5,0.5))
for(j in 1:4){
  i=low_list[j]
  rnge=index.mat[i,1]:index.mat[i,2]
  my.plot<-scatterplot3d(embb[,1],embb[,2],embb[,3],type="l",xlab=xlb ,ylab=ylb,zlab=zlb,color="gray",lwd=2,axis = TRUE, box = FALSE, main = paste(indexs[i],'-', indexs[i]+38, sep=" "))
  my.plot$points3d(embb[rnge,1],embb[rnge,2],embb[rnge,3],col=low_pred_col,type="l",lwd=2,lty = ltys)
  my.plot$points3d(embb[(index.mat[i,1]),1],embb[(index.mat[i,1]),2],embb[(index.mat[i,1]),3],col=low_pred_col,type="p",lwd=2,lty = ltys)
  mtext(side=3,paste(labels[j]),adj=0,line = -1.5,cex=1.3)
}
dev.off()

#figure windows - avg LLE, mean, var 
png(paste(graphdir,"leesferry-long_ts_comp.png",sep = ""), width = 9, height = 9, units='in', res = 300)
par(mfcol=c(3,1), mar=numeric(4), oma = c(6, 6, 1, 4),  mgp = c(2, .6, 0)) #oma(bottom,left,top,right)
#par(mfcol=c(3,1), mar=numeric(4), oma = c(6, 6, 4, 4),  mgp = c(2, .6, 0)) #oma(bottom,left,top,right)
plot(leesferry_long_df[,yr_col], leesferry_long_df[,avg_LLE_col], type = 'l', 
     #panel.first = rect(c(LLE_lowerbound1, LLE_lowerbound2, LLE_lowerbound3, LLE_lowerbound4), -2, c(LLE_upperbound1, LLE_upperbound2, LLE_upperbound3, LLE_upperbound4), 32, col=high_pred_col, border=NA),
     panel.first = rect(c(LLE_lowerbound5, LLE_lowerbound6, LLE_lowerbound7, LLE_lowerbound8), -2, c(LLE_upperbound5, LLE_upperbound6, LLE_upperbound7, LLE_upperbound8), 32, col=low_pred_col, border=NA),
     ylim = c(avglle_l,avglee_u),
     axes = FALSE, lwd = 2)
lines(leesferry_long_df[,yr_col], leesferry_long_df[,avg_LLE_col],col='black', lwd=2, 
      panel.first = rect(c(LLE_lowerbound1, LLE_lowerbound2, LLE_lowerbound3, LLE_lowerbound4), -2, c(LLE_upperbound1, LLE_upperbound2, LLE_upperbound3, LLE_upperbound4), 32, col=high_pred_col, border=NA))
abline(h=-0.05,lty=2)
abline(h=0.05,lty=2)
axis(2, tck = -0.01, cex.axis = 1.5, las = 2)
#mtext(side=3, "greenriverut", line = 2, cex = 1.3)
mtext(side=2, "Average LLE", line = 3, cex = 1.3)
mtext(side=4, "high pred.",adj = 0,line =1, cex = 1.1)
mtext(side=4, "low pred.",adj = 1,line =1, cex = 1.1)
mtext(side=3,"(a)",adj=0,line = -1.5,cex=1.3)
box()

plot(leesferry_long_df[,yr_col], leesferry_long_df[,roll_mean_col], type = 'l', 
     panel.first = rect(c(LLE_lowerbound5, LLE_lowerbound6, LLE_lowerbound7, LLE_lowerbound8), -2, c(LLE_upperbound5, LLE_upperbound6, LLE_upperbound7, LLE_upperbound8), 32, col=low_pred_col, border=NA),
     #panel.first = rect(c(LLE_lowerbound1, LLE_lowerbound2, LLE_lowerbound3, LLE_lowerbound4), -2, c(LLE_upperbound1, LLE_upperbound2, LLE_upperbound3, LLE_upperbound4), 32, col=high_pred_col, border=NA),
     ylim = c(mean_l,mean_u),
     axes = FALSE, lwd = 2)
lines(leesferry_long_df[,yr_col], leesferry_long_df[,roll_mean_col],col='black', lwd=2, 
      panel.first = rect(c(LLE_lowerbound1, LLE_lowerbound2, LLE_lowerbound3, LLE_lowerbound4), -2, c(LLE_upperbound1, LLE_upperbound2, LLE_upperbound3, LLE_upperbound4), 32, col=high_pred_col, border=NA))
axis(2, cex.axis = 1.5, tck = -0.01, las =2)
mtext(side=2, "Flow Average (MaF)", line = 3, cex = 1.3)
mtext(side=3,"(b)",adj=0,line = -1.5,cex=1.3)
box()

plot(leesferry_long_df[,yr_col], leesferry_long_df[,roll_var_col], type = 'l', 
     ylim = c(var_l,var_u),
     #panel.first = rect(c(LLE_lowerbound1, LLE_lowerbound2, LLE_lowerbound3, LLE_lowerbound4), -2, c(LLE_upperbound1, LLE_upperbound2, LLE_upperbound3, LLE_upperbound4), 32, col=high_pred_col, border=NA),
     panel.first = rect(c(LLE_lowerbound5, LLE_lowerbound6, LLE_lowerbound7, LLE_lowerbound8), -2, c(LLE_upperbound5, LLE_upperbound6, LLE_upperbound7, LLE_upperbound8), 32, col=low_pred_col, border=NA),
     axes = FALSE, lwd = 2)
lines(leesferry_long_df[,yr_col], leesferry_long_df[,roll_var_col],col='black', lwd=2, 
      panel.first = rect(c(LLE_lowerbound1, LLE_lowerbound2, LLE_lowerbound3, LLE_lowerbound4), -2, c(LLE_upperbound1, LLE_upperbound2, LLE_upperbound3, LLE_upperbound4), 32, col=high_pred_col, border=NA))
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


#####
#make labels for area fraction plots 
llabels = c("(c)", "(e)", "(f)", "(b)", "(d)", "(a)", "(g)")
#Multiple figures
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
  
  #figure with avgLLE, rolling average and rolling variance
  #png(paste(graphdir, name, "_ts_comp.png", sep=""), width = 9, height = 9, units='in', res = 300)
  par(mfcol=c(3,1), mar=numeric(4), oma = c(6, 6, 4, 4),  mgp = c(2, .6, 0)) #oma(bottom,left,top,right)
  plot(data_df[,yr_col], data_df[,avg_LLE_col], type = 'l', 
       ylim = c(avglle_l,avglee_u),
       axes = FALSE, lwd = 2)
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
       ylim = c(mean_l,mean_u),
       axes = FALSE, lwd = 2)
  axis(2, cex.axis = 1.5, tck = -0.01, las =2)
  mtext(side=2, "Flow Average (MaF)", line = 3, cex = 1.3)
  mtext(side=3,"(b)",adj=0,line = -1.5,cex=1.3)
  box()
  
  plot(data_df[,yr_col], data_df[,roll_var_col], type = 'l', 
       ylim = c(var_l,var_u),
       axes = FALSE, lwd = 2)
  axis(2, cex.axis = 1.5, tck = -0.01, las =2)
  axis(1, cex.axis = 1.5, tck = -0.01, las =1)
  mtext(side=2, "Flow Variance (MaF)", line = 3, cex = 1.3)
  mtext(side=1, "Years", line = 3, cex = 1.3)
  mtext(side=3,"(c)",adj=0, line = -1.5, cex=1.3)
  box()
  
  #dev.off()
}
