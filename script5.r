############################################################
##run the within-host model to get D and Vtot/S
#does the loop over alpha to determine overall-tradeoff
#the computes fitness for every alpha-gamma combination
############################################################
rm(list=ls()) #this clears the workspace to make sure no leftover variables are floating around. Not strictly needed
graphics.off(); #close all graphics windows
require(deSolve)  #loads ODE solver package

###################################################################
#function that specificies the ode model 
#by setting q=0 it's the simple model without immune response
###################################################################
odeequations=function(t,y,parms) #lsoda requires parms to be there, even though it's empty
{ 
  Utc=y[1]; Itc=y[2]; Vir=y[3]; B=y[4];  #uninfected cells, infected cells, virus
  lV=y[7]; #just trying to track log of virus, not really needed/used
  
  dUtcdt=-b*Vir*Utc;
  dItcdt=b*Vir*Utc-d*Itc;
  dVirdt=p*Itc-cw*Vir-q*Vir*B;
  dB=r*B;
  dVc=cw*Vir; #tracking virus removal due to clearance and IR killing
  dVq=q*Vir*B;
  dlV=0; #s*(p*Itc-cw*exp(lV))/(lV*exp(lV))
  
  return(list(c(dUtcdt,dItcdt,dVirdt,dB,dVc,dVq,dlV))); 
  
  
} #end function specifying the ODEs

###################################################################  
#function specifying events that happen during ODE integration, 
#namely stopping integration when virus load is <1 
###################################################################
rootfun <- function(t, y, x)
{
  (y[3]-Vlow); 
}    


###################################################################
#main program
###################################################################
  tmax=100;              
  timevec=seq(0,tmax,0.01); #vector of times for which integration is evaluated, units of days 
	filename="fitness-tradeoff.Rdata";

  filename.load="decay-fitting-results.Rdata"
  xx=load(file = filename.load); #load data from fitting script

  #values for within-host model parameters, units are assumed to be 1/days
  #data comes from within-host best fit
	
  #set initial conditions 
	Utc0=2.5e7; #initial number of uninfected cells  
	Itc0=0; #initial number for free virus V 
  #values for within-host model parameters, units are assumed to be 1/days
  #data comes from within-host best fit
  Vir0=1.000000e+00
  b=2.220917e-06
  p=1.178816e+00
  d=1.906049e+00

	Vir0=1; #initial number for free virus V 
  #immune response components
  r=1;
  B0=1; #initial B-cell/AB immune response
  Y0=c(Utc0, Itc0, Vir0,B0,0,0,log10(Vir0));  #combine initial conditions into a vector 


  #level of virus at which integration stops
	Vlow=1e-1;


  #c(T) trade-off values
	kappa=-0.262282;
	eta=-3.283455;

  alp.length=200; 
  avec=10^seq(-4,2,length=alp.length);
  gamvec=exp(eta+kappa*log(avec));
  cbvec=avec*exp(gamvec*5); #5 degrees    - units in 1/day
  cwvec=avec*exp(gamvec*40); #40 degrees   - units in 1/day
  
  R0vec=p*b*Utc0/(cwvec*d)
  R0tmp1=R0vec[1:which.max(R0vec)] #split R0vec into 2 to extract R0<1 on both the high and low end
  R0tmp2=R0vec[which.max(R0vec):length(R0vec)]
  alow=0.0005; ahigh=13; #avec values above/below no infection happens. found 'by hand' 
  clow=cwvec[which.min(R0tmp1<1)]
  chigh=cwvec[length(R0tmp1)-1+which.max(R0tmp2<1)]
  ccut=33; #level of cw where R0=1. Found "by hand".
  
######################################################
#Figure to show alpha-gamma decay plots
######################################################
#make a plot of c(T=5) and c(T=40) as function of alpha
#also add values from strains

graphics.off(); #close all graphics windows

ww=6.83; wh=3;
windows(width=ww, height=wh)   #open a window for plotting the results
par(mfrow=c(1,2))
par(mar=c(3, 3.5, 0.5, 0.5)) #bottom, left, top, right margins
mult=1; mult.ax=1; mult.lab=1; mult.text=0.6; #multiplier for plots
plot(resmatrix[,2],c.at.5,xlab="",ylab="",pch=18,col="white",log="xy",ylim=c(2e-3,1),xlim=c(avec[1],1),cex=mult,cex.axis=mult.ax,cex.lab=mult.lab)
text(resmatrix[,2],c.at.5,labels=strainnames,col=c(rep("black",4),"red",rep("black",7)),cex=mult.text)
lines(avec,cbvec,lty=1,lwd=2,type="l") 
mtext(expression(paste("intercept, ",alpha)),side=1,line=2,cex=1) #x-axis - change line to move in/out
mtext(expression(paste("environmental decay rate, ", c[b])),side=2,line=2,cex=1) #y-axis
text(1e-4,1,'A)',font=2);

par(mar=c(3, 3, 0.5, 0.6)) #bottom, left, top, right margins
plot(resmatrix[,2],c.at.40,xlab="",ylab="",pch=18,col="white",log="xy",ylim=c(2e-1,40),xlim=c(avec[1],tail(avec,1)),cex.lab=2,cex=mult,cex.axis=mult.ax,cex.lab=mult.lab)
text(resmatrix[,2],c.at.40,labels=strainnames,col=c(rep("black",4),"red",rep("black",7)),cex=mult.text)
lines(avec,cwvec,lty=1,lwd=2,type="l") 
lines(avec,rep(ccut,alp.length),lty=3,lwd=2)
mtext(expression(paste("intercept, ",alpha,"")),side=1,line=2,cex=1) #x-axis - change line to move in/out
mtext(expression(paste("within-host decay rate, ", c[w])),side=2,line=2,cex=1) #y-axis
text(1e-4,40,'B)',font=2);

dev.print(device=tiff,filename ="fig7.tiff",width=ww, height=wh, units="in",pointsize = 12,compression =c("lzw"), res=300)

qvec=c(0,1e-1) #with and without IR

shed1=matrix(0,nrow=alp.length,ncol=length(qvec)); #avec, all IR strength values 
shed2=shed1; shed3=shed1; #the s_i terms
#shed1a=rep(0,12); shed2a=shed1; shed3a=shed1; #the s_i terms

R0.d.s1=shed1; R0.i.s1=shed1;
R0.d.s2=shed1; R0.i.s2=shed1; #R0 values
R0.d.s3=shed1; R0.i.s3=shed1; #R0 values
fit.d.s1=shed1;  fit.i.s1=shed1;  
fit.d.s2=shed1;  fit.i.s2=shed1;  
fit.d.s3=shed1;  fit.i.s3=shed1;  

R0within=shed1; #within-host R0, just for diagnostics

killratio=shed1;  t.final=shed1;


for (ctq in 1:length(qvec))  #loop over different values of q
{
  
  for (ct1 in 1:alp.length) #loop over alpha
	{
	
    cb=cbvec[ct1]; #5 degrees    - units in 1/day
    cw=cwvec[ct1]; #40 degrees   - units in 1/day
    q=qvec[ctq]
    
    odeoutput=ode(y=Y0,times=timevec,func=odeequations,parms="",method="lsoda",rootfun=rootfun,atol=1e-5); #run and stop when rootfun happens
    tvec=odeoutput[,1];
    vvec=odeoutput[,4];
    lvvec=odeoutput[,8];
    
    max.ind=which.max(vvec); #find peak
    if (vvec[max.ind]<=1) {cat("infection did not take off!")}
    d.ind=max.ind+min(which(vvec[-(1:max.ind)]<=Vlow)); #find 1st time virus is below 1 after peak
    
    #final time for each integration
    t.final[ct1,ctq]=tvec[d.ind]
    
    #percentage of total virus removal due to IR 
    killratio[ct1,ctq]=100*odeoutput[d.ind,7]/(odeoutput[d.ind,6]+odeoutput[d.ind,7]);
     

    #different forms for total virus/shedding
    #integral of total virus load
    inf.t=tvec[1:d.ind];
    inf.v=vvec[1:d.ind];
    shed1[ct1,ctq]=max(0,sum(diff(inf.t)*(inf.v[-length(inf.v)]+inf.v[-1])/2))   
    
     #compute shedding value, i.e. function under integral
     cc1=5; cc2=5; cc3=2.5; 
     shed=inf.v * (cc1*pmax(0,log10(pmax(0,inf.v)))^cc2/(cc3^cc2+pmax(0,log10(pmax(0,inf.v)))^cc2))
     shed2[ct1,ctq]=sum(diff(inf.t)*(shed[-length(shed)]+shed[-1])/2) ;
         
     #total log virus load as area under the curve 
     shed=pmax(0,log10(pmax(0,inf.v)));
     shed3[ct1,ctq]=sum(diff(inf.t)*(shed[-length(shed)]+shed[-1])/2) ;
   
    print(sprintf('q is %f, R0 within-host is %f, final time %f, IR kill percentage %f, shedding: %f, %f, %f',q,R0vec[ct1],t.final[ct1,ctq],killratio[ct1,ctq],shed1[ct1,ctq],shed2[ct1,ctq],shed3[ct1,ctq]));
    
    #"R0" (not exactly, numerator and denominator for F) for 3 different shedding, direct and indirect transmission
    R0.d.s1[ct1,ctq]=shed1[ct1,ctq]       
    R0.i.s1[ct1,ctq]=shed1[ct1,ctq]/cbvec[ct1];       
    R0.d.s2[ct1,ctq]=shed2[ct1,ctq];
    R0.i.s2[ct1,ctq]=shed2[ct1,ctq]/cbvec[ct1];
    R0.d.s3[ct1,ctq]=shed3[ct1,ctq];
    R0.i.s3[ct1,ctq]=shed3[ct1,ctq]/cbvec[ct1];
 
  }   #finish loop over alpha


  fit.d.s1[,ctq]=R0.d.s1[,ctq]/max(R0.d.s1[,ctq]);
  fit.i.s1[,ctq]=R0.i.s1[,ctq]/max(R0.i.s1[,ctq]);   
  fit.d.s2[,ctq]=R0.d.s2[,ctq]/max(R0.d.s2[,ctq]);
  fit.i.s2[,ctq]=R0.i.s2[,ctq]/max(R0.i.s2[,ctq]);   
  fit.d.s3[,ctq]=R0.d.s3[,ctq]/max(R0.d.s3[,ctq]);
  fit.i.s3[,ctq]=R0.i.s3[,ctq]/max(R0.i.s3[,ctq]);   
  
  print(sprintf('for q = %f, min/max kill percentage is %f/%f',q,min(killratio[,ctq]),max(killratio[,ctq])));
  
} #finish loop over q


  #create  figures
  
  ######################################################
  #Fig 8 main text
  ######################################################
  ctq=1; #1st figure, no IR
  #plot the fitness for indirect/direct and 3 forms of shedding
   
  graphics.off(); #close all graphics windows
  ww=6.83; wh=3;
  windows(width=ww, height=wh)   #open a window for plotting the results
  par(mar=c(3, 3, 1, 10)) #bottom, left, top, right margins
  mult=1; mult.ax=1; mult.lab=1; #multiplier for plots
  plot(avec,fit.d.s1[,ctq],type="l",ylim=c(0.1,1),xlim=c(4e-4,20),lty=1,col="black",xlab="",ylab="",log="x",lwd=3,cex=mult,cex.axis=mult.ax,cex.lab=mult.lab,mar=c(3, 0.5, 0.5, 8))
  lines(avec,fit.i.s1[,ctq],type="l",lty=2,col="black",cex=mult,lwd=3)
  lines(avec,fit.d.s2[,ctq],type="l",lty=1,col="red",cex=mult,lwd=2)
  lines(avec,fit.i.s2[,ctq],type="l",lty=2,col="red",cex=mult,lwd=2)
  lines(avec,fit.d.s3[,ctq],type="l",lty=1,col="blue",cex=mult,lwd=2)
  lines(avec,fit.i.s3[,ctq],type="l",lty=2,col="blue",cex=mult,lwd=2)
  lines(rep(alow,100),10^seq(-5,4,length=100),lty=3,lwd=2)
  lines(rep(ahigh,100),10^seq(-5,4,length=100),lty=3,lwd=2)
  mtext(expression(paste("intercept, ",alpha)),side=1,line=2,cex=mult) #x-axis - change line to move in/out
  mtext("normalized fitness",side=2,line=2,cex=mult) #Y-axis - change line to move in/out
 
  ltext=c(expression(paste(s[1], " - direct")),expression(paste(s[1], " - environmental")),expression(paste(s[2], " - direct")),expression(paste(s[2], " - environmental")),expression(paste(s[3], " - direct")),expression(paste(s[3], " - environmental")) )
  legend(30,0.9,ltext,col=c("black","black","red","red","blue","blue"),lty=c(1,2,1,2,1,2),cex=mult,lwd=2,bty="n", xpd=TRUE)

 
  dev.print(device=tiff,filename ="fig8.tiff",width=ww, height=wh, units="in",pointsize = 12,compression =c("lzw"), res=300)
  
  
  #plot the fitness for indirect/direct and 3 forms of shedding
  ######################################################
  #Fig SM 3
  ######################################################

  ctq=2; #2nd figure, high IR
  graphics.off(); #close all graphics windows
  ww=6.83; wh=3;
  windows(width=ww, height=wh)   #open a window for plotting the results
  par(mar=c(3, 3, 1, 10)) #bottom, left, top, right margins

  plot(avec,fit.d.s1[,ctq],type="l",ylim=c(0.1,1),xlim=c(4e-4,20),lty=1,col="black",xlab="",ylab="",log="x",lwd=3,cex=mult,cex.axis=mult.ax,cex.lab=mult.lab,mar=c(3, 0.5, 0.5, 8))
  lines(avec,fit.i.s1[,ctq],type="l",lty=2,col="black",cex=mult,lwd=3)
  lines(avec,fit.d.s2[,ctq],type="l",lty=1,col="red",cex=mult,lwd=2)
  lines(avec,fit.i.s2[,ctq],type="l",lty=2,col="red",cex=mult,lwd=2)
  lines(avec,fit.d.s3[,ctq],type="l",lty=1,col="blue",cex=mult,lwd=2)
  lines(avec,fit.i.s3[,ctq],type="l",lty=2,col="blue",cex=mult,lwd=2)
  lines(rep(alow,100),10^seq(-5,4,length=100),lty=3,lwd=2)
  lines(rep(ahigh,100),10^seq(-5,4,length=100),lty=3,lwd=2)

  mtext("normalized fitness",side=2,line=2,cex=mult) #Y-axis - change line to move in/out
  mtext(expression(paste("intercept, ",alpha)),side=1,line=2,cex=mult) #x-axis - change line to move in/out

   ltext=c(expression(paste(s[1], " - direct")),expression(paste(s[1], " - environmental")),expression(paste(s[2], " - direct")),expression(paste(s[2], " - environmental")),expression(paste(s[3], " - direct")),expression(paste(s[3], " - environmental")) )
  legend(30,0.9,ltext,col=c("black","black","red","red","blue","blue"),lty=c(1,2,1,2,1,2),cex=mult,lwd=2,bty="n", xpd=TRUE)
  
  dev.print(device=tiff,filename ="S3.tiff",width=ww, height=wh, units="in",pointsize = 12,compression =c("lzw"), res=300)
  
###################################################################
#end main program
###################################################################                       


