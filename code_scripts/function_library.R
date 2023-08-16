remove_outliers <- function(x, na.rm = TRUE, ...) {
  qnt <- quantile(x, probs=c(.25, .75), na.rm = na.rm, ...)
  H <- 1.5 * IQR(x, na.rm = na.rm)
  y <- x
  y[x < (qnt[1] - H)] <- NA
  y[x > (qnt[2] + H)] <- NA
  y
}
########################################
proj.NLts <-function(X,ref,prj,n.sim=200,sig=0.9,dta=2){
  tser=X; tser[,2]=scale(tser[,2])
  mn=mean(X[,2])
  sdd=sd(X[,2])
  # get the signal and noise compponents of the time series
  decomp=decompose2(tser[,2],sig)
  signal=decomp$sig; sc.signal=rowSums(decomp$sig.sc.comp);sig.sawp=rowSums(decomp$sig.sawp);signal.components=decomp$sig.comp
  noise=decomp$noise; sc.noise=rowSums(decomp$noise.sc.comp);noise.sawp=rowSums(decomp$noise.sawp);noise.components=decomp$noise.comp
  Year=X[,1]
  ################# Simulate!!! ###################
  ref.year=ref
  sim.index<- 1:which(Year==ref.year)
  prj.index<- which(Year==ref.year+1):(which(Year==ref.year)+prj)
  sim.yr<-Year[sim.index]
  prj.yr<-Year[prj.index]
  
  ## get the Tau and d of the simulation time series and simulate
  
  Tau=timeLag(signal[sim.index],method='mutual',plot.data=T)[1]
  fnn=FNN(signal[sim.index], dimension=10, tlag=Tau, olag=1)[1,]
  plot(fnn,t='b');d=which(fnn==0)[1];d
  
  if(dta==1){
    xx=signal[sim.index];yy=noise[sim.index]
  }
  if(dta==2){
    xx=signal[-prj.index];yy=noise[-prj.index]
  }
  if(dta==3){
    xx=signal;yy=noise
  }
  
  emb=embedd(signal[sim.index],d,Tau)
  sim=sim.nldynamics(X=xx,X.mat=emb,Tau=Tau,prj,n.sim)
  sig.sim=sim$sim
  sig.sim.median=apply(sig.sim,1,median)
  
  nemb=embedd(noise[sim.index],d,Tau)
  nsim=sim.nldynamics(X=yy,X.mat=nemb,Tau=Tau,prj,n.sim)
  noise.sim=nsim$sim
  noise.sim.median=apply(noise.sim,1,median)
  
  sim.ts=sig.sim+noise.sim
  proj.ts=apply(sim.ts,1,median)
  
  simm=(sim.ts*sdd)+mn
  
  list(sim=simm,ssim=sig.sim,nsim=noise.sim)
}

####################################################################
sim.nldynamics=function(X,X.mat,Tau,prj,n.sim){
  n_sim=2*n.sim
  d=length(X.mat[1,])
  
  ## *********Sampling function*******************
  
  sample.index=function(index,n.sim){
    nn=length(unique(index))
    xx=1
    count=0
    while(xx < nn){
      sam=sample(index,n.sim,replace=T)
      xx=length(unique(sam))
      count=count+1
      if(count>10000000) stop("too many iterations")
    }
    return(sam)
  }
  
  vec=length(X)-(d-1)*Tau-1
  max.index=vec-prj
  sim.index=NLkindex(X=X.mat,Y=X,Tau=Tau,nsim=n_sim)
  # remove index values greater than max. index
  X.index=sim.index[sim.index<=max.index]
  #length(X.index)
  
  
  res.index=sample.index(X.index,n.sim)
  pr.ind=matrix(NaN,prj,n.sim)
  for(i in 1:n.sim){
    pr.ind[,i]=seq(res.index[i],(res.index[i]+prj-1))
  }
  
  sim.vec= X[((d-1)*Tau+2):length(X)]
  simX=sim.vec[pr.ind]
  X.sim=matrix(simX,prj,n.sim)
  #X.sim
  Result=list(sim=X.sim,simindex=pr.ind)
  return(Result)
  #   
}
##################################
NLkindex = function(X,Y,Tau,nsim)
{
  # X is the feature matrix: here it is the multivariate version of the time series
  # Y is the time series from which resampling is taking place
  # d is embedding dimension
  # Tau is the delay time
  # nsim = number of simulations
  # sim.year is the starting year of projection
  
  # This function require updateFM function
  
  # 1. Determine the weight matrix
  d=length(X[1,])
  #   if(d!=dim){
  #     print("Please insert the correct dimension")
  #     stop
  #   }
  N=length(X[,1])
  N1 = N-1
  K=round(sqrt(N1))
  W=1:K
  W=1/W
  W=W/sum(W)
  W=cumsum(W)
  
  n = 1 #initial point
  
  index=c()
  for ( i in 1:nsim){
    # 2. calculate the euclidean distance
    edist=as.matrix(dist(scale(X)))
    #dist = edist[1:length(edist[1,])]
    dist=edist[N,]
    # 3. select the first K nearest neighbors
    xdist = order(dist)[1:K] 
    
    # 4. Resample based on the closest distance
    # 4.1  Select the  euclidean distance for the current time
    #dist = edist[N,] 
    
    # 4.2 pick a uniformly distributed random number(udrn) between 0 and 1
    udrn=runif(1,0,1)
    
    # 4.3 rank the random number with the weighting metrics
    xy=c(udrn,W)
    xx=rank(xy)
    
    # 4.4 Get the resampling index
    indx=xdist[xx[1]]
    
    # 5. create the y vector to resample from
    begin=2+(d-1)*Tau
    Y.vec = Y[begin:length(Y)]
    # 6. resample one point and update the feature vector 
    #         yhat= Y.vec[indx]  
    #         
    #         end = begin  
    #         initial = begin-(d-1)*Tau
    #         data=Y[initial:end]
    #       
    #         #update the last row of feature vector
    #         update.row=updateFM(data,d,Tau)
    #         X[N,]=update.row
    index=c(index,indx) 
    #     # 7. resample from record
    
    yres=Y[index]
    
  }
  index
}
###################################
get.NL.params <-function(X,is.sig="yes",sig=0.9,dj=0.25){
  #X=scale(X)
  
  library('tseriesChaos')
  library('scatterplot3d')
  library('fractal')
  library(sm)
  
  if(is.sig == 'yes'){
    
    comps=decompose2(X,sig,dj)
    signal=comps$sig
    
    Tau=as.numeric(timeLag(signal,method='mutual',plot.data=T))#[1]
    
    fnn=FNN(signal, dimension=10, tlag=Tau, olag=1)[1,]
    plot(fnn,t='b', xlab='Embedding dimensions',ylab='% FNN')
    d=which(fnn==0)[1]
  }else{
    Tau=as.numeric(timeLag(X,method='mutual',plot.data=T))#[1]
    
    fnn=FNN(X, dimension=10, tlag=Tau, olag=1)[1,]
    plot(fnn,t='b', xlab='Embedding dimensions',ylab='False Nearest Neighbors (%)')
    d=which(fnn==0)[1]
    
  }
  
  list(t=Tau,d=d)
}


###########################################################

local.lyapunov <- function(X,is.sig='yes',sig=0.9,dj=0.28,plot="T",ttle=""){
  Year=X[,1]
  library('tseriesChaos')
  params=get.NL.params(X[,2],is.sig,sig,dj)
  Tau=params$t
  d=params$d;d
  
  if(is.sig == 'yes'){
    comps=decompose2(X[,2],sig,dj)
    signal=comps$sig
    saw=comps$ssawp
  }else{
    signal=X[,2]
  }
  
  n.reference=10
  N=length(signal)
  strt=(Tau*d)
  Ne=length(signal)-(d-1)*Tau
  end=Ne-n.reference
  indx=seq((strt),(end-1))
  scale.max = Ne - 2 - n.reference
  sc=end - strt
  lyp=lyap_k(signal, m=d, d=Tau, t=0, k=1, ref=n.reference, s=sc, eps=4)
  len=length(indx)
  if(plot=="T"){
    par(mar=c(5,4,4,4))
    plot(Year[indx],lyp,t='l',ylab='',xlab='Year',col=1,lwd=2,las=1,main=ttle,xlim=range(Year));abline(h=0,lty=2,lwd=2)
    
    mtext("Lyapunov Exponent",side=2,line=3,col=1)
    par(new=T)
    plot(Year,saw,axes=F,xlab="",ylab="",t='l',col=2,lwd=2)
    axis(side=4)
    mtext("Signal SAWP",side=4,line=2,col=2)
    abline(h=0,lty=2)
  }
  
  list(lyap=lyp, sawpp=saw,inx=indx)
}

######################################################################
ts.simplot <- function(obs,sim,indx,conf=5,ttle,cexx=1.5,adjj=0,ylm=c(0,50),unit="Flow (MAF)",statt=median){
  
  # Get the quantiles
  n=length(sim[,1])
  qts=matrix(0,length(obs),2)
  prob=c(conf/100,(1-conf/100))
  for(i in 1:n){
    qts[i,]=quantile(sim[i,],prob)
  }
  
  # Get the simulation median
  sim.median=apply(sim,1,statt)
  
  # Plot
  plot(Year[indx],sim.median,t='l',lwd=2,col=4,ylim=ylm,ylab=" ",xlab=" ")#, xlab="Year", )
  axis(side=2)
  mtext(unit,side=2,line=2)
  mtext("Year",side=1,line=2)
  mtext(ttle, side=3, adj=adjj, cex=cexx)
  y=c(qts[,1],rev(qts[,2]))
  x=c(Year[indx],rev(Year[indx]))
  
  polygon(x,y, col=8, border=NA)
  lines(Year[indx],sim.median,lwd=2,col=4,ylim=c(-3,3))
  
  lines(Year[indx],obs,col=2,lwd=2)
  
  
  # Add Legend
  
  lgd=c("Simulation median","Observed","Simulations")  
  l.wid=c(2,2,2)
  colr=c(4,2,8)  
  legend("topleft",col=colr,lwd=l.wid,legend=lgd,bty="n",cex=1)  
  sim.median
}
###########################################################################

boxplots<-function(obs,sim,yr.id){
  
  boxplot(t(simm[,,yr.id]));lines(obss[,yr.id],col=2) 
}
#######################################################

pdf=function(obs,sim,conf=0,unit='flow (MAF)',title="",cexx=1.5,adjj=0,ylm=c(0,0.2), is.emsemble='yes'){
  require(sm)
  require(Rlab)
  nsim=length(sim[1,])
  n=length(obs)
  
  std=sd(obs)
  xeval=seq(min(obs)-0.25*std,max(obs)+0.25*std,length=n)
  nevals=length(xeval)
  
  simpdf=matrix(0,n,nsim)
  for(k in 1:nsim){
    simpdf[,k]=sm.density(sim[,k],eval.points=xeval,display="none")$estimate
  }
  
  den=sm.density(obs,eval.points=xeval,display="none")$estimate
  simline=apply(simpdf,1,median)
  
  
  
  # Get the quantiles
  n=length(simpdf[,1])
  qts=matrix(0,length(obs),2)
  prob=c(conf/100,(1-conf/100))
  for(i in 1:n){
    qts[i,]=quantile(simpdf[i,],prob)
  }
  
  
  
  
  if(is.emsemble=='yes'){
    # plot the pdf
    nsim=length(simpdf[1,])
    par(font.main=4,font.axis=2,font.lab=2)
    plot(simpdf[,1],col=8,ylab="",xlab=unit,xaxt="n",ylim=ylm,t='l')
    #     for(i in 1:nsim){
    #       lines(simpdf[,i], col=8,lwd=2)
    #     }
    
    y=c(qts[,1],rev(qts[,2]))
    x=c(round(1:n),rev(1:n))
    
    polygon(x,y, col=8, border=NA)
    
  }else{
    boxplot(t(simpdf), col="gray",ylab="",xlab=unit,lab.font=2,xaxt="n",ylim=ylm)#,main=title)
  }
  
  #mtext(paste(unit),side=1,line=2)
  mtext("Density",side=2,line=2)
  mtext(title, side=3, adj=adjj, cex=cexx)
  axis(1,at=1:n,labels=round(xeval,2))
  lines(den,col="red",lwd=2)
  lines(simline,col="blue",lwd=2)
  
  ## add legend
  
  ltype=c(1,1,lwd=2)
  lgd=c("simulation median","Observed")
  colr=c("blue","red")
  legend("topleft",lgd,lty=ltype,col=colr,lwd=2,cex=1,bty="n")
}  

################################################

stats=function(obs,sim,unit="MAF"){
  require(sm)
  #   require(Rlab)
  #   library(Rlab)
  #   library(moments)
  nsim=length(sim[1,])
  n=length(obs)
  
  
  simmean=simvar=simskew=simlag1=1:nsim
  
  
  for(isim in 1:nsim){
    
    zsim=sim[,isim]
    simmean[isim]=median(zsim,na.rm=T)
    simvar[isim]= var(zsim,na.rm=T)
    simskew[isim]=skew(zsim)
    
    simlag1[isim]=cor(zsim[1:(n-1)], zsim[2:n])
    
  }
  
  ###  stats
  
  par(mfrow=c(1,4),mar=c(2,4,2,1),font.main=4,font.lab=2,font.axis=2)
  
  boxplot(simmean, main="Mean",ylab=unit,ylim=c(min(simmean,mean(obs)),max(simmean,mean(obs))))
  points(median(obs,na.rm=T),lwd=6,col="red")
  
  boxplot(simvar, main="Variance",ylab=paste(unit,'^2'),ylim=c(min(simvar,var(obs)),max(simvar,var(obs))))
  points(var(obs,na.rm=T),lwd=6,col="red")
  
  boxplot(simskew, main="Skew")
  points(skew(obs),lwd=6,col="red",ylim=c(min(simskew,skew(obs)),max(simskew,skew(obs))))
  
  boxplot(simlag1, main="lag1",ylim=c(min(simlag1,cor(obs[1:(n-1)], obs[2:n])),max(simlag1,cor(obs[1:(n-1)], obs[2:n]))))
  points(cor(obs[1:(n-1)], obs[2:n]),lwd=6,col="red")
}


#############################################
drandsurp=function(obs,sim,units="flow (MAF)",plot="yes"){
  # obs= observed data
  # sim= simulated data
  treshold = median(obs)   # treshold 
  ## simulated data  
  surp= drght = sim#; sim[sim<0]=0
  surp.index = which(sim > treshold)   # surplus index
  dr.index   = which(sim< treshold)  # drought index
  
  surp[dr.index]=NaN   # surplus
  drght[surp.index]=NaN   # drought
  
  surplus=surp-treshold  # surplus value
  drought=treshold - drght # drought value
  
  ## stats for simulated data
  # A. Total
  spls.total=colSums(surplus,na.rm=T) # surplus total
  drght.total=colSums(drought,na.rm=T) # drought total
  
  # B. Maximum
  surp.max=apply(surplus,2,max,na.rm=T)  # surplus maximum
  drght.max=apply(drought,2,max,na.rm=T)  # drought maximum
  
  # C. Mimimum
  surp.min=apply(surplus,2,min,na.rm=T)  # surplus minimum
  drght.min=apply(drought,2,min,na.rm=T)  # drought minimum
  
  
  # Observed data
  
  osurp=obs[obs>treshold] # surplus 
  odrght =obs[obs<treshold]  # drought
  osurplus=osurp-treshold # surplus value
  odrought=treshold-odrght 
  
  # stats for observed data
  # A. Total
  osurp.total = sum(osurplus)
  odrght.total = sum(odrought)
  # B. Maximum
  osurp.max=max(osurplus)
  odrght.max = max(odrought)
  # C. minimum
  osurp.min=min(osurplus)
  odrght.min = min(odrought)
  
  
  ### boxplots
  
  if(plot=="yes"){
    par(mfrow=c(2,3),font.main=4,font.lab=2,font.axis=2)
    
    boxplot(spls.total,main="Surplus Total",ylab=units)
    points(osurp.total,pch=19,col="red",lwd=4)
    
    boxplot(drght.total,main="Drought Total",ylab=units)
    points(odrght.total,pch=19,col="red",lwd=4)
    
    boxplot(surp.max,main="Surplus Maximum",ylab=units)
    points(osurp.max,pch=19,col="red",lwd=4)
    
    boxplot(drght.max,main="Drought Maximum",ylab=units)
    points(odrght.max,pch=19,col="red",lwd=4)
    
    boxplot(surp.min,main="Surplus Minimum",ylab=units)
    points(osurp.min,pch=19,col="red",lwd=4)
    
    boxplot(drght.min,main="Drought Minimum",ylab=units)
    points(odrght.min,pch=19,col="red",lwd=4)
    
    box()
    
  }
  
  
  list(stot=spls.total,smax=surp.max,smin=surp.min,dtot=drght.total,dmax=drght.max,dmin=drght.min,ostot=osurp.total,osmax=osurp.max,osmin=osurp.min,odtot=odrght.total,odmax=odrght.max,odmin=odrght.min)
  
}
######################################################

### sequent peak analysis 
sequent.peak=function(obs,sim,stat.par){
  n=length(obs)
  nsim=length(sim[1,])
  R=seq(13.5,20,by=1.5)    # demand
  l=length(R)
  
  sequent=function(x){  # sequent peak function
    kold=0
    K=matrix(NaN,l,n)
    for (r in 1:l){
      
      for(i in 1:n){
        knew=kold+R[r]-x[i]
        knew[knew<0]=0
        kold=knew
        K[r,i]=knew
      }
    }
    capacity=apply(K,1,stat.par)
    return(capacity)
    
  }
  
  sim.capacity= array(0,dim=c(nsim,l))
  
  for(i in 1:nsim){
    
    sim.capacity[i,]=sequent(sim[,i])
  }
  
  # observed data
  Obs.capacity=sequent(obs)
  
  par(mfrow=c(1,1))
  boxplot(sim.capacity,xaxt="n",xlab="Demand (MAF)", ylab=" Res. Capacity (MAF)")
  axis(1,at=1:length(R),labels= R)
  points(Obs.capacity,pch=19,col="red",lwd=6)
  
}
##########################################################
# Dominant Periodicities:

pdcty2=function(X,sig=0.95,dj=0.25){
  
  # requires wavelet function
  #requires significant level function
  
  wlt=wavelet(X,dj)
  local.power=wlt$power
  n=length(wlt$period)
  XL   = as.vector(X) # Format data
  signific = CI(sig,XL, "w",dj)
  
  
  
  divf=mean(wlt$power)/5
  DIVF=signif(divf, digits=2)
  
  pow=wlt$power/divf
  g.power=apply(pow,1,mean)
  #sig=signif$sig/divf
  
  sig.level = signific$sig/divf
  period = wlt$period
  sig.power = g.power; sig.power[sig.power<sig.level]=0
  xn=seq(1,length(sig.power))
  index = cbind(xn,sig.power)
  pds=index[(index[,2]==0),1]
  period[pds]=1000
  max.period=max(wlt$period)
  min.period = min(wlt$period)
  
  
  
  ind1=which((period)==1000);ind2=which((period)!=1000)
  # significant power
  sig.pwr.index=c(ind2[1],ind2[length(ind2)],ind2[which(diff(ind2)>1)],ind2[which(diff(ind2)>1)+1])
  sig.prd=wlt$period[sig.pwr.index[order(sig.pwr.index)]]
  yy=diff(sig.prd)
  indx.sig=which(yy==0)
  if(length(indx.sig)==0) sig.pd=matrix(sig.prd,ncol=2,byrow=T)
  if(length(indx.sig)>0) sig.pd=matrix(sig.prd[-c(indx.sig,(indx.sig+1))],ncol=2,byrow=T)
  
  
  
  # noise component
  if(length(indx.sig)>0) sigg=sig.prd[-c(indx.sig,(indx.sig+1))] else sigg=sig.prd
  
  noi.prd=c(min.period,sigg,max.period)
  xx=diff(noi.prd)
  indx=which(xx==0)
  if(length(indx)==0) noi.pd=matrix(noi.prd,ncol=2,byrow=T)
  if(length(indx)>0) noi.pd=matrix(noi.prd[-c(indx,(indx+1))],ncol=2,byrow=T)
  #noi.prd=noi.prd[-c(indx,(indx+1))]
  #noi.pd=matrix(noi.prd,ncol=2,byrow=T)
  
  list(pd=sig.pd,npd=noi.pd)
  
}


##########################################
#WAVELET  1D Wavelet transform with optional singificance testing
#
#   [WAVE,PERIOD,SCALE,COI] = wavelet(Y,DT,PAD,DJ,S0,J1,MOTHER,PARAM)
#
#   Computes the wavelet transform of the vector Y (length N),
#   with sampling rate DT.
#
#   By default, the Morlet wavelet (k0=6) is used.
#   The wavelet basis is normalized to have total energy=1 at all scales.
#
#
# INPUTS:
#
#    Y = the time series of length N.
#    DT = amount of time between each Y value, i.e. the sampling time.
#
# OUTPUTS:
#
#    WAVE is the WAVELET transform of Y. This is a complex array
#    of dimensions (N,J1+1). FLOAT(WAVE) gives the WAVELET amplitude,
#    ATAN(IMAGINARY(WAVE),FLOAT(WAVE) gives the WAVELET phase.
#    The WAVELET power spectrum is ABS(WAVE)^2.
#    Its units are sigma^2 (the time series variance).
#
#
# OPTIONAL INPUTS:
# 
# *** Note *** setting any of the following to -1 will cause the default
#               value to be used.
#
#    PAD = if set to 1 (default is 0), pad time series with enough zeroes to get
#         N up to the next higher power of 2. This prevents wraparound
#         from the end of the time series to the beginning, and also
#         speeds up the FFT's used to do the wavelet transform.
#         This will not eliminate all edge effects (see COI below).
#
#    DJ = the spacing between discrete scales. Default is 0.25.
#         A smaller # will give better scale resolution, but be slower to plot.
#
#    S0 = the smallest scale of the wavelet.  Default is 2*DT.
#
#    J1 = the # of scales minus one. Scales range from S0 up to S0*2^(J1*DJ),
#        to give a total of (J1+1) scales. Default is J1 = (LOG2(N DT/S0))/DJ.
#
#    MOTHER = the mother wavelet function.
#             The choices are 'MORLET', 'PAUL', or 'DOG'
#
#    PARAM = the mother wavelet parameter.
#            For 'MORLET' this is k0 (wavenumber), default is 6.
#            For 'PAUL' this is m (order), default is 4.
#            For 'DOG' this is m (m-th derivative), default is 2.
#
#
# OPTIONAL OUTPUTS:
#
#    PERIOD = the vector of "Fourier" periods (in time units) that corresponds
#           to the SCALEs.
#
#    SCALE = the vector of scale indices, given by S0*2^(j*DJ), j=0...J1
#            where J1+1 is the total # of scales.
#
#    COI = if specified, then return the Cone-of-Influence, which is a vector
#        of N points that contains the maximum period of useful information
#        at that particular time.
#        Periods greater than this are subject to edge effects.
#        This can be used to plot COI lines on a contour plot by doing:
#
#              contour(time,log(period),log(power))
#              plot(time,log(coi),'k')
#
#----------------------------------------------------------------------------
#   Copyright (C) 1995-2004, Christopher Torrence and Gilbert P. Compo
#
#   This software may be used, copied, or redistributed as long as it is not
#   sold and this copyright notice is reproduced on each copy made. This
#   routine is provided as is without any express or implied warranties
#   whatsoever.
#
# Notice: Please acknowledge the use of the above software in any publications:
#    ``Wavelet software was provided by C. Torrence and G. Compo,
#      and is available at URL: http://paos.colorado.edu/research/wavelets/''.
#
# Reference: Torrence, C. and G. P. Compo, 1998: A Practical Guide to
#            Wavelet Analysis. <I>Bull. Amer. Meteor. Soc.</I>, 79, 61-78.
#
# Please send a copy of such publications to either C. Torrence or G. Compo:
#  Dr. Christopher Torrence               Dr. Gilbert P. Compo
#  Research Systems, Inc.                 Climate Diagnostics Center
#  4990 Pearl East Circle                 325 Broadway R/CDC1
#  Boulder, CO 80301, USA                 Boulder, CO 80305-3328, USA
#  E-mail: chris[AT]rsinc[DOT]com         E-mail: compo[AT]colorado[DOT]edu
#----------------------------------------------------------------------------
wavelet=function(Y,dj=0.025){
  
  #Y is time series to be analyzed
  DT=1# is timestep for annual data, 1
  pad=1
  #dj=0.025
  param=6
  #pad data ? 0=F, 1=T
  #dj= spacing between discrete scales (.025)
  #param = wavenumber (6)
  
  s0=2*DT
  
  n1 = length(Y)
  J1=floor((log2(n1*DT/s0))/dj)
  
  
  #....construct time series to analyze, pad if necessary
  x = Y - mean(Y)
  
  
  if (pad == 1){
    base2 = trunc(log(n1)/log(2) + 0.4999)   # power of 2 nearest to N
    x = c(x, rep(0, 2^(base2 + 1) - n1))
  }
  n = length(x)
  
  #....construct wavenumber array used in transform [Eqn(5)]
  k = (1:trunc(n/2))
  k = k*((2*pi)/(n*DT))
  k = c(0, k, -rev(k[1:floor((n-1)/2)]))
  
  #....compute FFT of the (padded) time series
  f = fft(x)    # [Eqn(3)]
  
  #....construct SCALE array & empty PERIOD & WAVE arrays
  scale = s0*2^((0:J1)*dj)
  period = scale;
  wave = matrix(data=0, ncol=n, nrow=J1+1)  # define the wavelet array
  wave = as.complex(wave)  # make it complex
  wave=matrix(data=wave, ncol=n, nrow=J1+1)
  
  # loop through all scales and compute transform
  for(a1 in 1:(J1+1)){
    scl=scale[a1]  	
    
    nn = length(k);
    k0 = param
    expnt = -(scl*k - k0)^2/(2*(k > 0))
    norm = sqrt(scl*k[2])*(pi^(-0.25))*sqrt(nn)    # total energy=N   [Eqn(7)]
    daughter = norm*exp(expnt)
    daughter = daughter*(k > 0)    # Heaviside step function
    fourier_factor = (4*pi)/(k0 + sqrt(2 + k0^2)) # Scale-->Fourier [Sec.3h]
    coi = fourier_factor/sqrt(2)                  # Cone-of-influence [Sec.3g]
    dofmin = 2                                   # Degrees of freedom
    
    out <- list(daughter=daughter, fourier_factor=fourier_factor,coi=coi,dofmin=dofmin)
    
    daughter=out$daughter
    fourier_factor=out$fourier_factor
    coi=out$coi
    dofmin=out$dofmin	
    wave[a1,] = fft((f*daughter), inverse = TRUE)/(length(f*daughter))  # wavelet transform[Eqn(4)]
  }
  
  period = fourier_factor*scale
  
  coi = coi*c(1:(floor(n1 + 1)/2), rev(1:floor(n1/2))) * DT
  
  wave = wave[,1:n1]  # get rid of padding before returning
  power=abs(wave)^2
  ncol=length(power[1,])
  nrow=length(scale)
  avg.power=apply(power,1,mean)
  result=list(wave=wave, period=period, scale=scale, power=power, coi=coi,nc=ncol,nr=nrow,p.avg=avg.power)
  return(result)
}
###################################################################################
filterts=function(X, lp, up,dj=0.025){
  
  ## X is the time series to be decomposed
  ## lp is the lower period
  ## up is the upper Period
  
  ## Ref: 
  # A practical guide to Wavelet Analysis; Christopher. Torrence and Gilbert P. Compo
  # Nowak
  XD=X#-mean(X)
  n=length(X)
  
  #dj=0.025
  dt=1
  cd=0.776
  psy0 = pi^(-1/4)  ## page 67 table 2
  
  # source and run the wavelet function
  
  wlt=wavelet(X,dj)
  wave=wlt$wave
  scale=wlt$period
  
  L=length(scale)
  scale.mat=cbind(1:L,scale)
  scale.bounds=(scale.mat[,2]<=up & scale.mat[,2]>=lp)  # filter the lower and upper bounds
  filter.bound = (scale.mat[scale.bounds,])             # filtered bound
  index=filter.bound[,1]  # index
  l=length(index)
  
  xn=matrix(ncol=n, nrow=l)                             # individual time series component
  
  for( i in 1: l){
    row.no=index[i]
    #xn[i,] = (dj*sqrt(dt)/(cd*psy0))*(as.real(wave[row.no,])/sqrt(scale[row.no]))  #eq. 11, page 68
    xn[i,] = (dj*sqrt(dt)/(cd*psy0))*(Re(wave[row.no,])/sqrt(scale[row.no]))  #eq. 11, page 68
  }
  
  filtered.ts=colSums(xn)
  #result=list(filter=t(filtered.ts),comp=t(xn))
  return( filtered.ts)
  
}
###################################

sawp=function(X,lp,up,dj=0.025){
  
  ## Ref: A practical guide to Wavelet Analysis; Christopher. Torrence and Gilbert P. Compo
  XD=X-mean(X)
  n=length(X)
  #dj=0.25
  dt=1
  cd = 0.776;                                          #for morlet ; page 67, table 2
  constant=dj*dt/cd
  s0=2*dt
  J1=floor((log2(n*dt/s0))/dj)
  s=s0*2^(J1*dj)
  
  wlt=wavelet(XD,dj)
  wave=wlt$wave
  scale=wlt$period
  
  L=length(scale)
  scale.mat=cbind(1:L,scale)
  scale.bounds=(scale.mat[,2]<=up & scale.mat[,2]>=lp)  # filter the lower and upper bands
  filter.bound = (scale.mat[scale.bounds,])             # filtered band
  index=filter.bound[,1]                                # index no. to search from the output of wavelet function 
  l=length(index)
  
  wn=matrix(ncol=n, nrow=l)                             # individual time series component
  
  for( i in 1: l){
    row.no=index[i]
    wn[i,] = (Mod(wave[row.no,])^2/(scale[row.no]))     #eq. 24, page 73 " Averaging in Scale"
  }
  
  band=sqrt(colSums(wn))
  sawp=sqrt(constant)*band
  
  return(sawp)
  
  
}
############################
bknn <- function(X,lagg,nsim)
{
  #source block_sim_3
  ###########################################################################################################  
  lagknn = function(X,lagg,nsim)
  {
    rem=length(X)%%lagg
    if(rem==0){
      elements=seq((lagg+1),length(X),by=lagg)
      rem.sim=c()
      
    }else{
      elements=seq((lagg+1),(length(X)-rem),by=lagg)
      X=X[1:(length(X)-rem)]
      rem.sim=matrix(0,rem,nsim)
    }
    
    simmat=mav(X,2*lagg)$datas
    xprev=as.matrix(simmat[,1:lagg])
    xt=as.matrix(simmat[,-c(1:lagg)])
    
    N=length(xprev[,1])
    N1 = N-lagg
    
    ## number of simulations desired is nsim
    ## each simulation will be of length N
    simflow = matrix(0,ncol=nsim, nrow=length(X))
    
    ## the weight metric to do the K-NN resampling
    K=round(sqrt(N1))
    W=1:K
    W=1/W
    W=W/sum(W)
    W=cumsum(W)
    
    for(i in 1:nsim){
      ind=1
      #pick a flow at random to start
      fst.dist=order(as.matrix(dist(xprev))[ind,])
      xx=runif(1,0,1)
      xy=c(xx,W)
      xx=rank(xy)
      indx=fst.dist[xx[1]]
      simflow[1:lagg,i] = xprev[indx,]
      simm=c()
      
      xsam=xprev[1,]
      
      for(j in elements) {
        xsim=xprev
        # find the dist of xprev to all the N-1 years in the historical record
        #
        xsim[ind,]=xsam
        xdist = order(as.matrix(dist(xsim))[ind,])
        
        #select one of the nearest neighbor with the weight function
        # above.
        
        xx=runif(1,0,1)
        xy=c(xx,W)
        xx=rank(xy)
        index=xdist[xx[1]]
        
        ind=j
        xsam=xt[index,]
        simm=c(simm,xsam)
        
      }
      simflow[(lagg+1):length(X),i]=simm[1:(length(X)-lagg)]
      
    }
    
    simflow
  }
  
  simm=matrix(0,length(X),nsim) 
  rem=length(X)%%lagg
  if(rem==0){
    my.sim=lagknn(X,lagg,nsim)
    rem.sim=c()
    simm=my.sim
  }else{
    xprime= X[1:(length(X)-rem)]
    my.sim1=lagknn(xprime,lagg,nsim)
    n=length(X)-rem+1
    xprime2=c(X[n:length(X)],xprime[-c((n-rem):n)])
    rem.sim=lagknn(xprime2,rem,nsim)[1:rem,]
    my.sim=abind(my.sim1,rem.sim,along=1)
  }
  my.sim
}

#########################################

# moving average

mav=function(X,win,param='mean'){
  
  N=length(X)
  n=N-win+1
  datas=matrix(NaN,n,win)
  for(i in 1:n){ datas[i,]=X[i:(i+win-1)] }
  av=apply(datas,1,param)
  list(avg=av,datas=datas)
  
}

###########################
decompose2=function(X,sig=0.9,dj=0.25){
  library(abind)
  periodocities=pdcty2(X,sig,dj)
  pd=periodocities$pd
  noi.pd=periodocities$npd
  nr=length(X)
  
  ## Get the significant components
  nc=length(pd[,1])
  X.comps=X.sawp=X.sc.comp=matrix(NaN,nr,nc)
  for( i in 1:nc){
    X.comps[,i]=filterts(X,pd[i,1],pd[i,2],dj)
    X.sawp[,i]=sawp(X,pd[i,1],pd[i,2],dj)
    X.sc.comp[,i]=X.comps[,i]/X.sawp[,i]
  }
  
  # Get Noise components
  n.nc=length(noi.pd[,1])
  n.comps=n.sawp=n.sc.comp=matrix(NaN,nr,n.nc)
  for( j in 1:n.nc){
    n.comps[,j]=filterts(X,noi.pd[j,1],noi.pd[j,2],dj)
    n.sawp[,j]=sawp(X,noi.pd[j,1],noi.pd[j,2],dj)
    n.sc.comp[,j]=n.comps[,j]/n.sawp[,j]
  }
  if(nc>1){signal=rowSums(X.comps);ssawp=rowSums(X.sawp) }else {signal=X.comps; ssawp=X.sawp}
  if(n.nc>1){noise=rowSums(n.comps);nsawp=rowSums(n.sawp)} else {noise=n.comps; nsawp=n.sawp}
  mean.periodicities=rowMeans(abind(pd,noi.pd,along=1))
  list(sig=signal,noise=noise, sig.comp=X.comps,noise.comp=n.comps,sig.sawp=X.sawp,noise.sawp=n.sawp,sig.sc.comp=X.sc.comp,noise.sc.comp=n.sc.comp, ssawp=ssawp,nsawp=nsawp,prd=mean.periodicities)
}

############################
#decompose3=function(X,per1start=4,per1end=16,per2start=32,per2end=120,dj=0.25){
decompose3=function(X,Period = c(4,16,32,120),dj=0.25){
  library(abind)
  #pd=cbind(c(per1start,per2start),c(per1end,per2end))
  pd=matrix(Period, ncol=2, byrow=T)
  noise=c(1,Period,500)
  noi.pd=matrix(noise,ncol=2,byrow=T)
  #noi.pd=cbind(c(1,per1end,per2end),c(per1start,per2start,500))
  nr=length(X)
  
  ## Get the significant components
  nc=length(pd[,1])
  X.comps=X.sawp=X.sc.comp=matrix(NaN,nr,nc)
  for( i in 1:nc){
    X.comps[,i]=filterts(X,pd[i,1],pd[i,2],dj)
    X.sawp[,i]=sawp(X,pd[i,1],pd[i,2],dj)
    X.sc.comp[,i]=X.comps[,i]/X.sawp[,i]
  }
  
  # Get Noise components
  n.nc=length(noi.pd[,1])
  n.comps=n.sawp=n.sc.comp=matrix(NaN,nr,n.nc)
  for( j in 1:n.nc){
    n.comps[,j]=filterts(X,noi.pd[j,1],noi.pd[j,2],dj)
    n.sawp[,j]=sawp(X,noi.pd[j,1],noi.pd[j,2],dj)
    n.sc.comp[,j]=n.comps[,j]/n.sawp[,j]
  }
  if(nc>1){signal=rowSums(X.comps);ssawp=rowSums(X.sawp) }else {signal=X.comps; ssawp=X.sawp}
  if(n.nc>1){noise=rowSums(n.comps);nsawp=rowSums(n.sawp)} else {noise=n.comps; nsawp=n.sawp}
  mean.periodicities=rowMeans(abind(pd,noi.pd,along=1))
  list(sig=signal,noise=noise, sig.comp=X.comps,noise.comp=n.comps,sig.sawp=X.sawp,noise.sawp=n.sawp,sig.sc.comp=X.sc.comp,noise.sc.comp=n.sc.comp, ssawp=ssawp,nsawp=nsawp,prd=mean.periodicities)
}

############################

comps_dec=function(X,sig=0.9,dj=0.025){
  comps=decompose2(X,sig)
  signals=comps$sig.comp
  noises=comps$noise.comp
  prds=ceiling(comps$prd)
  my.compos=cbind(signals,noises)
  
  list(comp=signals,pr=prds)
}
########################################################################
# a function for component block projection
wknn_prj=function(X,bl,prjj,nsim,is.sig='yes'){
  if(is.sig=='yes'){
    compos=comps_dec(X,sig=0.9,dj=0.025);my.compos=as.matrix(compos$comp) 
    NN=length(my.compos[1,])
    prds=compos$pr[1:NN]
    my.simm=array(0,dim=c(prjj,nsim,NN))
    sim=matrix(0,prjj,nsim)
    for(i in 1:NN){
      prd=max(2,min(bl,as.integer(prds[i]/2)))
      bsim=prj(my.compos[,i],prd,prjj,nsim)
      #pdf(my.compos[,i],bsim)
      my.simm[,,i]=bsim$sim
      sim=sim+bsim$sim
      indx=bsim$indx
    }
  }else{
    simm=prj(X,bl,prjj,nsim)
    sim=simm$sim
    indx=simm$indx
  }
  list(sim=sim,inx=indx)
}


#############################
# Block projection function

prj <- function(X,bl,prjj,nsim){
  
  K.nearest<-function(n){
    
    K=floor(sqrt(n))
    w=1:K
    w=1/w
    w=w/sum(w)
    w=cumsum(w)
    w
  }
  #############################
  embedd=function(X,d,t){
    
    L=length(X)-(d-1)*t
    new.ts=matrix(NaN,L,d)
    for(i in 1:d){
      indx=1+(i-1)*t
      index=seq(indx,(indx+(L-1)))
      new.ts[,i]=X[index]
    }
    new.ts
  }
  
  #######################
  
  #if(bl < prjj){bl=prjj}
  N=length(X)
  NN=N-prjj
  
  # Embedd
  emb.data=embedd(X,2,bl)
  pre.dta=emb.data[,1]  # predictor block 
  predee.dta=emb.data[,2] # predictee block 
  
  
  
  xx=mav(pre.dta,prjj)$datas
  yy=mav(predee.dta,prjj)$datas
  
  #indxx=max(bl,prjj)
  N1=length(xx[,1])
  wt.metric=K.nearest(N)
  ## create matric
  mat <- rbind(yy[N1,],xx)
  
  eucl.dist=order(as.matrix(dist(scale(yy)))[1,][-1])
  
  # Generate ensemble
  sim=matrix(0,prjj,nsim)
  #inxs=matrix(0,prjj,nsim)
  for(i in 1:nsim){
    # Generate random number
    rand=runif(1,0,1)
    indx=eucl.dist[rank(c(rand,wt.metric))[1]]
    
    inxx <- rank(c(rand,wt.metric))[1]
    
    sim[,i]=yy[inxx,]
    
  }
  sim
  
}


######################
CI=function(conf, dat,type,dj){
  
  #### enter confidence as decimal 0-1
  #### two types of tests available 1) red noise enter: "r" , white noise enter: "w"
  # requires the wavelet function
  
  na=length(dat)
  wlt=wavelet(dat,dj)
  
  if(type=="r"){
    
    zz=arima(dat/10^10, order = c(1, 0, 0))#arima(dat/(10^10), order = c(1, 0, 0))
    alpha=zz$coef[1]
    print(alpha)
  } else{
    alpha=0
  }
  
  ps=wlt$period
  LP= length(ps)
  
  freq = 1/ps
  
  CI=1:LP    ## confidence interval
  
  for(i in 1:LP){
    
    P=(1-(alpha^2))/(1+(alpha^2)-(2*alpha*cos(2*pi*freq[i])))    # descrete fourier power spectrum page 69 [qn 16] ( torrence and compo)... null hypothesis test
    df=2*sqrt(1+((na/(2.32*ps[i]))^2))
    CI[i] =P*(qchisq(conf, df)/df)*var(dat)          #divide by 2 removes power of 2.....for mean no chi dist.[ eqn 17]
  }
  
  #plot(ps,(CI), log="x", xlim=rev(range(ps)))
  
  list(sig=CI)
  
}
#####################################
fcontour <- function (x = seq(0, 1, length.out = nrow(z)), y = seq(0, 1, 
                                                                   length.out = ncol(z)), z, left=TRUE, xlim = range(x, finite=TRUE), 
                      ylim = range(y, finite = TRUE), zlim = range(z, finite = TRUE), 
                      levels = pretty(zlim, nlevels), nlevels = 20, color.palette = rainbow, 
                      col = color.palette(length(levels) - 1), plot.title, plot.axes, 
                      xaxs = "i", yaxs = "i", las = 1, axes = TRUE, frame.plot = axes, ...) 
{
  if (missing(z)) {
    if (!missing(x)) {
      if (is.list(x)) {
        z <- x$z
        y <- x$y
        x <- x$x
      }
      else {
        z <- x
        x <- seq.int(0, 1, length.out = nrow(z))
      }
    }
    else stop("no 'z' matrix specified")
  }
  else if (is.list(x)) {
    y <- x$y
    x <- x$x
  }
  
  if (any(diff(x) <= 0) || any(diff(y) <= 0)) 
    stop("increasing 'x' and 'y' values expected")
  
  mar.orig <- (par.orig <- par(c("mar", "las")))$mar
  par(las = las)
  mar <- mar.orig
  on.exit(par(par.orig))
  if (left){
    mar[4]<- 1	
  }
  else{
    mar[2]<- 1
    mar[4]<- 4.1
  }
  par(mar=mar)
  plot.new()
  
  
  plot.window(xlim, ylim, "", xaxs = xaxs, yaxs = yaxs)
  if (!is.matrix(z) || nrow(z) <= 1 || ncol(z) <= 1) 
    stop("no proper 'z' matrix specified")
  if (!is.double(z)) 
    storage.mode(z) <- "double"
  # 	.Internal(filledcontour(as.double(x), as.double(y), z, as.double(levels), 
  # 		col = col))
  .filled.contour(as.double(x), as.double(y), z, as.double(levels), 
                  col = col)
  if (missing(plot.axes)) {
    if (axes && left) {
      title(main = "", xlab = "", ylab = "")
      Axis(x, side = 1)
      Axis(y, side = 2)
    }
    else if(axes){
      title(main = "", xlab = "", ylab = "")
      Axis(x, side = 1)
      Axis(y, side = 4)	
    }
  }
  else plot.axes
  if (frame.plot) 
    box()
  if (missing(plot.title)) 
    title(...)
  else plot.title
  invisible()
}



cbar <- function ( nlevels, zlim, rotate=TRUE, levels = pretty(zlim, nlevels), 
                   color.palette =
                     rainbow, col = color.palette(length(levels) - 1), 
                   key.title, key.axes, xaxs = "i", yaxs = "i", las = 1, axes = TRUE) 
{
  mar.orig <- (par.orig <- par(c("mar", "las")))$mar
  par(las = las)
  mar <- mar.orig
  on.exit(par(par.orig))
  plot.new()
  if(rotate){
    mar[4] <- 2
    mar[2] <- 4
    par(mar = mar)
    plot.window(xlim = c(0, 1), ylim = range(levels), xaxs = "i", 
                yaxs = "i")
    rect(0, levels[-length(levels)], 1, levels[-1], col = col)
    if (missing(key.axes)) {
      if (axes) 
        axis(2)
    }
    else key.axes
  }
  else {
    mar[1] <- mar[3]
    mar[3] <- 2
    par(mar = mar)
    plot.window(xlim = range(levels), ylim = c(0,1), xaxs = "i", 
                yaxs = "i")
    rect( levels[-length(levels)], 0, levels[-1], 1, col = col)
    if (missing(key.axes)) {
      if (axes) 
        axis(1)
    }
    else  key.axes
  }
  box()
  if (!missing(key.title)) 
    key.title
}



sizer <- function(rotate=TRUE){
  
  mar <- par("mar")
  
  if(rotate){		
    w <- (3 + mar[2]) * par("csi") * 2.54
  }
  else{
    w <- (3+mar[3]) * par("csi") * 2.54
  }
  return(w)
}
######################################################

wplot=function(X,start,dj=0.025,u.sig=0.95,l.sig=0.9,name='Power_Spectrum', graphdir = graphdir){
  # X= time series to be analyzed
  # start = start year of the time series
  
  # requires function: fcontour
  #X=X/10^6
  library(fields)
  library(DescTools) #to find the area below the curve
  yr = seq(start,start+length(X)-1)  # No. of years for the starting year of "start"
  
  data=matrix(0,length(X),2)
  
  data[,1]=yr
  
  data[,2]=X
  
  wlt=wavelet(X,dj)
  
  #lvl  = sig #significant level
  
  XL   = as.vector(X) # Format data
  
  xx=CI(u.sig,XL, "w",dj)
  xxx=CI(l.sig,XL,"w",dj)
  #xx=noise(lvl,XL, "w",wlt,length(X))
  
  #########first row from wave/power is smallest period
  
  units="MaF"
  #units="acre-ft"
  #units="index"
  
  coi=wlt$coi
  
  divf=mean(wlt$power)/5
  
  DIVF=signif(divf, digits=2)
  
  pow=wlt$power/(DIVF)
  
  
  my.levels <- quantile(pow, probs = seq(0, 1, len=64))
  
  my.cols = tim.colors(64)
  #my.cols = rainbow(64)
  
  p=t(pow)
  
  for( i in 1:length(p[,1])){
    p[i,]=rev(p[i,])
  }  
  
  xper=subset(wlt$period, wlt$period<(.5*length(data[,1])))
  
  y=1:length(xper)
  
  del=length(wlt$period)-length(xper)
  
  P=p[,(1+del):(length(wlt$period))]
  
  label=c(2,4,8,16,32,64,128,256,512)
  
  l=length(wlt$power[1,])
  
  val=1:length(wlt$period)
  
  for( i in 1:length(wlt$period)){
    
    val[i]=mean(pow[i,])
    
    
  }
  
  ts=1:l
  
  VAL=val
  
  ### plot the local and global power spectrums with key bar
  
  #device = png
  
  #ext = "png"
  
  png(paste(graphdir,name,".png", sep = ""), width=1400, height=600)
  
  par(mfrow=c(1,3), cex=2,cex.axis = 2.5, cex.lab =2.5)
  
  layout(matrix(c(1,2,3),ncol=3, nrow=1),widths = c(.69,.07,.25))#c(.575,.1,.325)
  
  par(mar=c(6,6,6,1))
  
  fcontour(data[,1], log(wlt$period), t(pow), levels= my.levels,col=my.cols, ylim=c(log(2), max(log(wlt$period))), nlevels=5,ylab="", xlab="",
           
           plot.axes = {
             
             axis(1,cex.axis=2.5)
             axis(2 , at = log(label), labels = label, cex.axis=2.5)
             
           }
  )
  
  contour(data[,1], log(wlt$period), t(pow), nlevels=5,levels=(xx$sig[1]/DIVF), add=T, drawlabels=F,lwd=1)
  
  ## plot Cone of Influence
  lines(data[,1], log(wlt$coi), lty=2,lwd=2,col="black")
  mtext("Period (Years)", side = 2, line = 3.7, cex=2)
  
  mtext("Year", side = 1, line = 3, cex=2)
  mtext(side=3,"(b)",adj=0, line =1, cex=2.5)
  ### plot key bar
  n=seq(1,length(my.levels),len=5)
  
  lbl=my.levels[n]
  
  zz=seq(0,length(my.levels)-1)
  
  cbar(nlevels=64,col=my.cols, levels=zz, key.title = "",key.axes = {axis(2,at =n , labels = round(lbl, digits=1), cex.axis=2.5)})
  
  mtext(parse(text=paste("(",units,"^2)")), side = 2, adj=.47,line = 5)
  
  mtext("Power /", side = 2, adj=.35,line = 5) #.27
  
  mtext(paste(DIVF), side = 2, adj=.41,line = 5)
  
  
  ### plot Global power
  
  par(usr=c(min(VAL), max(VAL), min(log(wlt$period)), max(log(wlt$period))))
  #sig.level = xxx$sig/DIVF #original 
  sig.level = xx$sig/DIVF #modified for getting the upper sig 
  sig.power = VAL; sig.power[sig.power>sig.level]=NaN #original
  
  #modified to get what we want #################################3
  #this is the time series above the upper sig level 
  sig.power_upper = VAL; sig.power_upper[sig.power_upper<=sig.level]=NaN
  #this is all the time series 
  sig.power_all = VAL #################################3
  
  #original 
  #plot(sig.power,log(wlt$period), ylab="", xlab="", type="l", ylim=range(log(wlt$period)), xlim=range(VAL,xx$sig/DIVF),lwd=2,main=, axes=F, yaxs="i")#,col="blue")
  
  #modified to show the area fraction 
  plot(sig.power_all,log(wlt$period), ylab="", xlab="", type="l", ylim=range(log(wlt$period)), xlim=range(VAL,xx$sig/DIVF),lwd=2,main=, axes=F, yaxs="i")#,col="blue")
  
  #lines(sig.power,log(wlt$period),col="red",lwd=2)
  axis(1, cex.axis=2.5)
  
  axis(2 , at = log(label), labels = label,las=2, cex.axis=2.5)
  
  box()
  mtext(side=3,"(c)",adj=0,cex=2.5, line =1)
  mtext("Period (Years)", side = 2, line = 4, cex=2)
  
  mtext("Global Power", side = 1, adj=.2,line = 3, cex=2)
  
  mtext(parse(text=paste("(",units,"^2)")), side = 1, adj=.85,line =3.2, cex=2) #.8
  
  
  ### plot confidence level
  lines(xx$sig/DIVF, log(wlt$period), col="red",lty=2,lwd=2)
  
  ## add significant power with different color
  sig.period = log(wlt$period)
  
  wvlt_per = wlt$period #needed for area calculation 
  
  #modifying to get the area 
  area_list = which(is.na(sig.power_upper)==FALSE)
  split_list = split(area_list, cumsum(c(1, diff(area_list) != 1)))
  
  area = AUC(log(wvlt_per), sig.power_all, method = "trapezoid")
  area_sections = c()
  
  for(split in split_list){
    i_start = split[1]; i_end = split[length(split)]
    polygon(c(0,0,sig.power_all[i_start], sig.power_all[i_end], sig.power_all[i_start:i_end]), 
            c(sig.period[i_end], sig.period[i_start],sig.period[i_start], sig.period[i_end],sig.period[i_start:i_end]), col = "gray", border = NA)
    area_section <- AUC(log(wvlt_per[i_start:i_end]), sig.power_upper[i_start:i_end],  method = "trapezoid")
    area_sections[length(area_sections) + 1] <- area_section
  }
  
  area_fraction = sum(area_section)/area
  
  #return this value as part of the wPlot()
  
  
  
  sig.level1 = xx$sig/DIVF
  sig.level2 = xxx$sig/DIVF
  sig.power1 = VAL; sig.power1[sig.power1<sig.level1]=NaN
  sig.power2 = VAL; sig.power2[sig.power2<sig.level2]=NaN; sig.power2[sig.power2>sig.level1]=NaN
  
  
  # add legend
  ltype=c(2) 
  
  displvlu=u.sig*100
  #displvll=l.sig*100
  
  lgd=c(paste(Re(displvlu), "% confidence level"))#, "% confidence level") 
  
  colr=c("red")#"darkgreen")
  
  legend("bottomright",lgd,lty=ltype,col=colr,lwd=c(2,2),cex=2.5,bty="n")
  
  dev.off()
  
  #Get the period ranges for significant signals
  sig.power2.save = VAL; sig.power2.save[sig.power2.save<sig.level2]=NaN;
  periods=cbind(sig.power2.save, wlt$period)
  
  #Result=list(per=periods, sig.power2=sig.power2) #original 
  #modified 
  Result=list(sig.power_upper = sig.power_upper, sig.power_all = sig.power_all, wlt = wlt, sig.level1=sig.level1, sig.period = sig.period, area_fraction = area_fraction)
  return(Result)
  
}
#############################
wplot_area=function(X,start,dj=0.025,u.sig=0.95, name='area-fraction', graphdir = graphdir, letter_label = '(a)'){
  # requires fcontour and sig_label functions
  # X= time series to be analyzed
  # start = start year of the time series
  
  # requires function: fcontour
  #X=X/10^6
  library(fields)
  library(DescTools) #to find the area below the curve
  yr = seq(start,start+length(X)-1)  # No. of years for the starting year of "start"
  
  data=matrix(0,length(X),2)
  
  data[,1]=yr
  
  data[,2]=X
  
  wlt=wavelet(X,dj)
  
  #lvl  = sig #significant level
  
  XL   = as.vector(X) # Format data
  
  xx=CI(u.sig,XL, "w",dj)
  #xxx=CI(l.sig,XL,"w",dj)
  #xx=noise(lvl,XL, "w",wlt,length(X))
  
  #########first row from wave/power is smallest period
  
  units="MaF"
  #units="acre-ft"
  #units="index"
  
  #coi=wlt$coi
  
  divf=mean(wlt$power)/5
  
  DIVF=signif(divf, digits=2)
  
  pow=wlt$power/(DIVF)
  
  
  #my.levels <- quantile(pow, probs = seq(0, 1, len=64))
  
  #my.cols = tim.colors(64)
  #my.cols = rainbow(64)
  
  p=t(pow)
  
  for( i in 1:length(p[,1])){
    p[i,]=rev(p[i,])
  }  
  
  xper=subset(wlt$period, wlt$period<(.5*length(data[,1])))
  
  y=1:length(xper)
  
  del=length(wlt$period)-length(xper)
  
  P=p[,(1+del):(length(wlt$period))]
  
  label=c(2,4,8,16,32,64,128,256,512)
  
  l=length(wlt$power[1,])
  
  val=1:length(wlt$period)
  
  for( i in 1:length(wlt$period)){
    
    val[i]=mean(pow[i,])
    
    
  }
  
  
  VAL=val
  
  ### plot the local and global power spectrums with key bar
  
  #device = png
  
  #ext = "png"
  
  png(paste(graphdir,name,"_AF.png",sep=""), width=4, height=6,units='in', res = 300)
  
  
  ### plot Global power
  
  par(usr=c(min(VAL), max(VAL), min(log(wlt$period)), max(log(wlt$period))))
  sig.level = xx$sig/DIVF #modified for getting the upper sig 
  sig.power = VAL; sig.power[sig.power>sig.level]=NaN #original
  
  #modified to get what we want #################################3
  #this is the time series above the upper sig level 
  sig.power_upper = VAL; sig.power_upper[sig.power_upper<=sig.level]=NaN
  #this is all the time series 
  sig.power_all = VAL #################################3
  
  #modified to show the area fraction 
  plot(sig.power_all,log(wlt$period), ylab="", xlab="", type="l", ylim=range(log(wlt$period)), xlim=range(VAL,xx$sig/DIVF),lwd=2,main=, axes=F, yaxs="i")#,col="blue")
  
  #lines(sig.power,log(wlt$period),col="red",lwd=2)
  axis(1, cex.axis=1.5)
  
  axis(2 , at = log(label), labels = label,las=2, cex.axis=1.5)
  
  box()
  mtext("Period (Years)", side = 2, line = 3, cex=1.5)
  mtext(text=paste("Global Power"), side = 1, adj=.3,line = 2.75, cex=1.5)
  mtext(parse(text=paste("(",units,"^2)")), side = 1, adj=1,line =3.2, cex=1.5) #.8

  #mtext(parse(text=paste("(",units,"^2)")), side = 1, adj=.85,line =3.2, cex=1) #.8
  
  
  ### plot confidence level
  lines(xx$sig/DIVF, log(wlt$period), col="red",lty=2,lwd=2)
  
  ## add significant power with different color
  sig.period = log(wlt$period)
  
  wvlt_per = wlt$period #needed for area calculation 
  
  #modifying to get the area 
  area_list = which(is.na(sig.power_upper)==FALSE)
  split_list = split(area_list, cumsum(c(1, diff(area_list) != 1)))
  
  area = AUC(log(wvlt_per), sig.power_all, method = "trapezoid")
  area_sections = c()
  
  for(split in split_list){
    i_start = split[1]; i_end = split[length(split)]
    polygon(c(0,0,sig.power_all[i_start], sig.power_all[i_end], sig.power_all[i_start:i_end]), 
            c(sig.period[i_end], sig.period[i_start],sig.period[i_start], sig.period[i_end],sig.period[i_start:i_end]), col = "gray", border = NA)
    area_section <- AUC(log(wvlt_per[i_start:i_end]), sig.power_upper[i_start:i_end],  method = "trapezoid")
    area_sections[length(area_sections) + 1] <- area_section
  }
  
  area_fraction = sum(area_sections)/area
  
  #return this value as part of the wPlot()
  mtext(paste(letter_label), side=3,adj=-0.15, line = 1, cex=1.5)
  mtext(paste(name,"AF =",round(area_fraction,4)), side = 3, line = 0, cex=1.5)
  
  sig.level1 = xx$sig/DIVF
  sig.power1 = VAL; sig.power1[sig.power1<sig.level1]=NaN
  
  # add legend
  ltype=c(2) 
  
  displvlu=u.sig*100
  
  
  lgd=c(paste(Re(displvlu),"% confidence level"))#, "% confidence level") 
  
  colr=c("red")#"darkgreen")
  
  legend("topright",lgd,lty=ltype,col=colr,lwd=c(2),cex=1.1,bty="n")
  
  dev.off()
  
  return(area_fraction)
  
}
#############################


## NL-sim from different starting points and conditions
proj.NLts <- function(X,ref,proj,n.sim=200,sig=0.95,is.signal="yes"){
  Year=X[,1]
  tser=X#; tser[,2]=scale(tser[,2])
  mn=mean(X[,2])
  sdd=sd(X[,2])
  
  # get the signal and noise compponents of the time series
  decomp=decompose2(tser[,2],sig)
  signal=decomp$sig; sc.signal=rowSums(decomp$sig.sc.comp);sig.sawp=rowSums(decomp$sig.sawp);signal.components=decomp$sig.comp
  noise=decomp$noise; sc.noise=rowSums(decomp$noise.sc.comp);noise.sawp=rowSums(decomp$noise.sawp);noise.components=decomp$noise.comp
  
  ##### Preresuisite Functions ####################
  
  my.knn=function(X,nsim){
    N=length(X[,1])
    K=round(sqrt(N))
    W=1:K
    W=1/W
    W=W/sum(W)
    W=cumsum(W)
    
    # calculate euclidean distance
    kdist=as.matrix(dist(X,method = "euclidean"))[N,]#[-N]
    
    indexx=(order(kdist))[1:K]
    index=c()
    for(i in 1:nsim){
      udrn=runif(1,0,1)
      #rank the random number among the weighting metrics
      xy=c(udrn,W)
      xx=rank(xy)[1]
      index=c(index,indexx[xx])
    }
    
    index
  }
  ################################
  ### NLKNN simulation
  NLKNN <- function(X,Tau,d,proj,nsim){
    N=length(X)
    #1. Embedd the timeseries
    x=embedd(X,Tau,d)
    n=N-(d-1)*Tau
    
    #2. Generate the resamplee vector
    start.index=(d-1)*Tau+1
    end.index=N
    xt=X[start.index:end.index]
    
    #3. reference time index
    ref.index=n-(proj+2)
    xx=rbind(x[-c(ref.index:n),],x[n,])
    
    #   ref.index=(d-1)*Tau+proj
    #   xx=x[-c(1:ref.index),]
    
    # 4. perform K-NN
    knn.index=my.knn(xx,nsim)
    # 5. Simulate !!!
    #create resampling matrix
    sim=matrix(0,proj,nsim)
    for(i in 1:nsim){
      sim[,i]=xt[seq(knn.index[i],(knn.index[i]+proj-1))]
    }
    sim
  }
  
  
  ################# Simulate!!! ###################
  ref.year=ref
  sim.index<- 1:which(Year==ref.year)
  prj.index<- which(Year==ref.year+1):(which(Year==ref.year)+proj)
  sim.yr<-Year[sim.index]
  prj.yr<-Year[prj.index]
  
  #   ## get the Tau and d of the simulation time series and simulate
  Tau=timeLag(signal[sim.index],method='mutual',plot.data=T)[1]
  fnn=FNN(signal[sim.index], dimension=10, tlag=Tau, olag=1)[1,]
  plot(fnn,t='b');d=which(fnn==0)[1];d
  
  if(is.signal=="yes"){
    
    # simulate the signal
    x=signal[sim.index]
    ssim=NLKNN(x,Tau,d,proj,n.sim)
    
    # simulate the noise ( randomly resample segment of noise with length = proj from the past data)
    # create matrix of noise with ncol=proj
    #noise.mat=mav(noise[sim.index],proj)$datas
    
    nsim= noise.sampler(noise[sim.index],proj)
    #     nn=length(noise.mat[,1])
    #     ## sample rows randomly
    #     #row.indx.sample=sample(nn,n.sim,replace=T)
    #     
    #     ## generate uniformly distributed random number
    #     ref.vec=seq(0,1,length=nn)
    #     zz=runif(n.sim,0,1)
    #     row.indx.sample=c()
    #     for( j in 1:n.sim){
    #       xy=c(zz[j],ref.vec)
    #       
    #       row.indx.sample=c(row.indx.sample,rank(xy)[1])
    #       
    #     }
    # #     nse=rnorm(proj*nsim,mean=mean(noise[sim.index]),sd=sd(noise[sim.index]))
    # #     nsim=matrix(nse,nrow=proj)
    #     nsim=t(noise.mat[row.indx.sample,])
    
    # add the signal and the noise and back standardize the simulation
    ts.sim=(ssim+nsim)+mn#*sdd +mn
  }
  
  if(is.signal=="no"){
    x=X[sim.index,2]
    ts.sim=NLKNN(x,Tau,d,proj,n.sim)
  }
  
  ts.sim
}


embedd=function(X,d,T){
  
  L=length(X)-(d-1)*T
  new.ts=matrix(NaN,L,d)
  for(i in 1:d){
    indx=1+(i-1)*T
    index=seq(indx,(indx+(L-1)))
    new.ts[,i]=X[index]
  }
  new.ts
}