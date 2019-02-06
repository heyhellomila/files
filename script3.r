############################################################
##run the within-host model to get D, Vpeak, and Vtot
#use that to compute fitness, both in presence and absence of virulence
############################################################
rm(list=ls()) #this clears the workspace to make sure no leftover variables are floating around. Not strictly needed
graphics.off(); #close all graphics windows
require(deSolve)  #loads ODE solver package


###################################################################
#function that specificies the ode model 
###################################################################
odeequations=function(t,y,parms) #lsoda requires parms to be there, even though it's empty
{ 
  	Utc=y[1]; Itc=y[2]; Vir=y[3];   #uninfected cells, infected cells, virus
  
	  dUtcdt=-b*Vir*Utc;
		dItcdt=b*Vir*Utc-d*Itc;
		dVirdt=p*Itc-cw*Vir;
   	
    return(list(c(dUtcdt,dItcdt,dVirdt))); 

    
} #end function specifying the ODEs

###################################################################
#main program
###################################################################
  tmax=50;              
  timevec=seq(0,tmax,0.001); #vector of times for which integration is evaluated, units of days 
	filename="fitness12strains.Rdata";
  
  #set initial conditions 
	Utc0=2.5e7; #initial number of uninfected cells  
	Itc0=0; #initial number for free virus V 
  #values for within-host model parameters, units are assumed to be 1/days
  #data comes from within-host best fit
  Vir0=1.000000e+00
  b=2.220917e-06
  p=1.178816e+00
  d=1.906049e+00
  
	Y0=c(Utc0, Itc0, Vir0);  #combine initial conditions into a vector 

	
  #within- and between-host decay rates (units of 1/day) for the 12 strains - values come from fitting decay rates and computing decay at 5C and 40C
	cwvec=c(0.9143109, 1.1466441, 0.7720865, 1.4266564, 5.3159759, 10.1449393, 0.6973284, 0.4382046, 1.9480118, 1.2631844, 7.3341671, 1.9187683);
  cbvec=c(0.030649683, 0.033049389, 0.024954676, 0.036237560, 0.120608903, 0.003221546, 0.049174557, 0.036928843, 0.026679227, 0.124843694, 0.024628383, 0.035103518);

  Duration=rep(0,12);  
  shed1=rep(0,12); shed2=shed1; shed3=shed1; #the s_i terms
  R0.d.s1=shed1; R0.i.s1=rep(0,12);
  R0.d.s2=shed1; R0.i.s2=shed1; #R0 values
  R0.d.s3=shed1; R0.i.s3=shed1; #R0 values
  
  R0star.d.s1=shed1; R0star.i.s1=shed1;
  R0star.d.s2=shed1; R0star.i.s2=shed1; #R0 values
  R0star.d.s3=shed1; R0star.i.s3=shed1; #R0 values
  
  
  R0within=shed1; #within-host R0, just for diagnostics

  for (strain in 1:12)
	{
	
    cw=cwvec[strain]
    R0within[strain]=p*b*Utc0/(cw*d)
    print(sprintf('R0 within-host is %f for strain %d',R0within[strain],strain));

    odeoutput=lsoda(Y0,timevec,odeequations,parms="",atol=1e-10,rtol=1e-10); 
    tvec=odeoutput[,1];
    vvec=odeoutput[,4];

    if (strain==1) {plot(tvec,vvec,log="y",type="l",lwd=2,xlim=c(0,tmax),ylim=c(1,Utc0)); }
    lines(tvec,vvec,type="l",lwd=2);
    lines(tvec,odeoutput[,2],type="l",lwd=2,col="blue"); 
    
    max.ind=which.max(vvec); #find peak
    d.ind=max.ind+min(which(vvec[-(1:max.ind)]<1)); #find 1st time virus is below 1 after peak

    Duration[strain]=tvec[d.ind];

    Vpeak=(max(vvec)); #peak of virus load

    #different forms for total virus/shedding
    #integral of total virus load
    inf.t=tvec[1:d.ind];
    inf.v=vvec[1:d.ind];
    shed1[strain]=sum(diff(inf.t)*(inf.v[-length(inf.v)]+inf.v[-1])/2)   
       
     #compute shedding value, i.e. function under integral
     cc1=5; cc2=5; cc3=2.5; 
     shed=inf.v * (cc1*pmax(0,log10(pmax(0,inf.v)))^cc2/(cc3^cc2+pmax(0,log10(pmax(0,inf.v)))^cc2))
     shed2[strain]=sum(diff(inf.t)*(shed[-length(shed)]+shed[-1])/2) ;
         
     #total log virus load as area under the curve 
     shed=pmax(0,log10(pmax(0,inf.v)));
     shed3[strain]=sum(diff(inf.t)*(shed[-length(shed)]+shed[-1])/2) ;
   
     
  	 #"R0" (not exactly, numerator and denominator for F) for 3 different shedding, direct and indirect transmission
     R0.d.s1[strain]=shed1[strain]       
     R0.i.s1[strain]=shed1[strain]/cbvec[strain];       
     R0.d.s2[strain]=shed2[strain];
     R0.i.s2[strain]=shed2[strain]/cbvec[strain];
     R0.d.s3[strain]=shed3[strain];
     R0.i.s3[strain]=shed3[strain]/cbvec[strain];
  
  	 #"R0" (not exactly, numerator and denominator for F) for 3 different shedding, direct and indirect transmission
     #also including virulence
     R0star.d.s1[strain]=shed1[strain]/Vpeak       
     R0star.i.s1[strain]=shed1[strain]/Vpeak/cbvec[strain];       
     R0star.d.s2[strain]=shed2[strain]/Vpeak;
     R0star.i.s2[strain]=shed2[strain]/Vpeak/cbvec[strain];
     R0star.d.s3[strain]=shed3[strain]/Vpeak;
     R0star.i.s3[strain]=shed3[strain]/Vpeak/cbvec[strain];
       
    
  }  #end loop over 12 strains
   
  #fitness values from the paper 
	fit.d.s1=R0.d.s1/R0.d.s1[1];
  fit.i.s1=R0.i.s1/R0.i.s1[1];   
  fit.d.s2=R0.d.s2/R0.d.s2[1];
  fit.i.s2=R0.i.s2/R0.i.s2[1];   
  fit.d.s3=R0.d.s3/R0.d.s3[1];
  fit.i.s3=R0.i.s3/R0.i.s3[1];   

  #fitness with virulence
	fitstar.d.s1=R0star.d.s1/R0star.d.s1[1];
  fitstar.i.s1=R0star.i.s1/R0star.i.s1[1];   
  fitstar.d.s2=R0star.d.s2/R0star.d.s2[1];
  fitstar.i.s2=R0star.i.s2/R0star.i.s2[1];   
  fitstar.d.s3=R0star.d.s3/R0star.d.s3[1];
  fitstar.i.s3=R0star.i.s3/R0star.i.s3[1];   

save(list = ls(all=TRUE), file = filename);
 
#plot ordered by clearance at 40c
cwsort=sort(cwvec,index.return=TRUE)
cw.ind=cwsort$ix
xlabels.raw=c("H1N1","H2N4","H3N2","H4N6","H5N2","H6N4","H7N6","H8N4","H9N2","H10N7","H11N6","H12N5");
xlabels=xlabels.raw[cw.ind] #resort labels according to cw

###########################
#figure for paper - figure 6 main text
###########################

#plot the fitness for different shedding - direct transmission
graphics.off()
ww=6.83; wh=6.83;
windows(width=ww,height=wh)
par(mfrow=c(2,1))
par(mar=c(2, 3.5, 0.5, 0.5)) #bottom, left, top, right margins
dx=0.1; #spacer between values
mult=1; mult.ax=1; mult.lab=1; mult.xlab=0.75; #multiplier for plots

sym.with=2;
ymax=max(c(fit.d.s1,fit.d.s2,fit.d.s3))
plot((1:12)-dx,fit.d.s1[cw.ind],type="p",ylim=c(0,ymax),pch=0,col="black",xlab="",ylab="",cex=mult,cex.axis=mult.ax,cex.lab=mult.lab,xaxt="n",lwd=sym.with)
points((1:12),fit.d.s2[cw.ind],type="p",pch=5,col="blue",cex=mult,lwd=sym.with)
points((1:12)+dx,fit.d.s3[cw.ind],type="p",pch=2,col="red",cex=mult,lwd=sym.with)
axis(1,at=1:12,lab=xlabels,cex.axis=mult.xlab);
mtext("relative fitness",side=2,line=2,cex=1.5) #y-axis
title("A) Direct transmission",line=-1) #y-axis
ltext=c(expression(paste('shedding ~ viral load, ', s[1])),expression(paste('shedding ~ total discharge, ',s[2])),expression(paste('shedding ~ log viral load, ',s[3])) )
legend(7,1.9,ltext,col=c("black","blue","red"),pch=c(0,5,2),cex=mult,bty="n",lwd=sym.with,lty=0)

par(mar=c(2, 3.5, 0.5, 0.5)) #bottom, left, top, right margins
dx=0.1; #spacer between values
ymax=max(c(fit.i.s1,fit.i.s2,fit.i.s3))
plot((1:12)-dx,fit.i.s1[cw.ind],type="p",ylim=c(0,ymax),pch=0,col="black",xlab="",ylab="",cex=mult,cex.axis=mult.ax,cex.lab=mult.lab,xaxt="n",lwd=sym.with)
points((1:12),fit.i.s2[cw.ind],type="p",pch=5,col="blue",cex=mult,lwd=sym.with)
points((1:12)+dx,fit.i.s3[cw.ind],type="p",pch=2,col="red",cex=mult,lwd=sym.with)
axis(1,at=1:12,lab=xlabels,cex.axis=mult.xlab);
mtext("relative fitness",side=2,line=2,cex=1.5) #y-axis
title("B) Environmental transmission",line=-1) #y-axis

ltext=c(expression(paste('shedding ~ viral load, ', s[1])),expression(paste('shedding ~ total discharge, ',s[2])),expression(paste('shedding ~ log viral load, ',s[3])) )
legend(1,5,ltext,col=c("black","blue","red"),pch=c(0,5,2),cex=mult,bty="n",lwd=sym.with,lty=0)

dev.print(device=tiff,filename ="fig6.tiff",width=ww, height=wh, units="in",pointsize = 12,compression =c("lzw"), res=300)

###########################
#figure for paper - virulence figure for supplementary material, S4
###########################

graphics.off()
ww=6.83; wh=6.83;
windows(width=ww,height=wh)
par(mfrow=c(2,1))
par(mar=c(2, 3.5, 0.5, 0.5)) #bottom, left, top, right margins
dx=0.1; #spacer between values
mult=1; mult.ax=1; mult.lab=1; mult.xlab=0.75; #multiplier for plots
ymax=max(c(fitstar.d.s1,fitstar.d.s2,fitstar.d.s3))
plot((1:12)-dx,fitstar.d.s1[cw.ind],type="p",ylim=c(0,ymax),pch=0,col="black",xlab="",ylab="",cex=mult,cex.axis=mult.ax,cex.lab=mult.lab,xaxt="n",lwd=sym.with)
points((1:12),fitstar.d.s2[cw.ind],type="p",pch=5,col="blue",cex=mult,lwd=sym.with)
points((1:12)+dx,fitstar.d.s3[cw.ind],type="p",pch=2,col="red",cex=mult,lwd=sym.with)
axis(1,at=1:12,lab=xlabels,cex.axis=mult.xlab);
mtext("relative fitness",side=2,line=2,cex=1.5) #y-axis
title("A) Direct transmission",line=-1) #y-axis
                     
ltext=c(expression(paste('shedding ~ viral load, ', s[1])),expression(paste('shedding ~ total discharge, ',s[2])),expression(paste('shedding ~ log viral load, ',s[3])) )
legend(1,5,ltext,col=c("black","blue","red"),pch=c(0,5,2),cex=mult,bty="n")

par(mar=c(2, 3.5, 0.5, 0.5)) #bottom, left, top, right margins
dx=0.1; #spacer between values
ymax=max(c(fitstar.i.s1,fitstar.i.s2,fitstar.i.s3))
ymax=10;
plot((1:12)-dx,fitstar.i.s1[cw.ind],type="p",ylim=c(0,ymax),pch=0,col="black",xlab="",ylab="",cex=mult,cex.axis=mult.ax,cex.lab=mult.lab,xaxt="n",lwd=sym.with)
points((1:12),fitstar.i.s2[cw.ind],type="p",pch=5,col="blue",cex=mult,lwd=sym.with)
points((1:12)+dx,fitstar.i.s3[cw.ind],type="p",pch=2,col="red",cex=mult,lwd=sym.with)
axis(1,at=1:12,lab=xlabels,cex.axis=mult.xlab);
mtext("relative fitness",side=2,line=2,cex=1.5) #y-axis
title("B) Environmental transmission",line=-1) #y-axis
ltext=c(expression(paste('shedding ~ viral load, ', s[1])),expression(paste('shedding ~ total discharge, ',s[2])),expression(paste('shedding ~ log viral load, ',s[3])) )
legend(1,9.5,ltext,col=c("black","blue","red"),pch=c(0,5,2),cex=mult,bty="n")

dev.print(device=tiff,filename ="S4.tiff",width=ww, height=wh, units="in",pointsize = 12,compression =c("lzw"), res=300)


  
###################################################################
#end main program
###################################################################                      
