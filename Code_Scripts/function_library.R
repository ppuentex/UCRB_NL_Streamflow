#these functions were written by 
# Rajagopalan, B., Erkyihun, S. T., Lall, U., Zagona, E., & Nowak, K. (2019). 
# A nonlinear dynamical systems‐based modeling approach for stochastic simulation of streamﬂow and understanding predictability.
# Water Resources Research, 55, 6268–6284. https://doi.org/10.1029/2018WR023650

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
