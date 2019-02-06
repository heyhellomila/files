############################################################
##run the within-host model to get D and Vtot
#use that to compute fitness
#same as in script 3, but this one includes an immune response in the within-host model
############################################################
rm(list=ls()) #this clears the workspace to make sure no leftover variables are floating around. Not strictly needed
graphics.off(); #close all graphics windows
require(deSolve)  #loads ODE solver package


###################################################################
#function that specificies the ode model 
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
	

  #set initial conditions 
	Utc0=2.5e7; #initial number of uninfected cells  
	Itc0=0; #initial number for free virus V 
  #values for within-host model parameters, units are assumed to be 1/days
  #data comes from within-host best fit
  Vir0=1.000000e+00
  b=2.220917e-06
  p=1.178816e+00
  d=1.906049e+00
  
  #immune response components
	r=1;
	B0=1; #initial B-cell/AB immune response
	Y0=c(Utc0, Itc0, Vir0,B0,0,0,log10(Vir0));  #combine initial conditions into a vector 

  #level of virus at which integration stops
	Vlow=1e-1;
	
  #within- and between-host decay rates (units of 1/day) for the 12 strains
	cwvec=c(0.9143109, 1.1466441, 0.7720865, 1.4266564, 5.3159759, 10.1449393, 0.6973284, 0.4382046, 1.9480118, 1.2631844, 7.3341671, 1.9187683);
  cbvec=c(0.030649683, 0.033049389, 0.024954676, 0.036237560, 0.120608903, 0.003221546, 0.049174557, 0.036928843, 0.026679227, 0.124843694, 0.024628383, 0.035103518);

  qvec=c(1e-5,1e-4,1e-3,1e-2,1e-1)

  shed1=matrix(0,nrow=12,ncol=length(qvec)); #12 strains, all IR strength values 
  shed2=shed1; shed3=shed1; #the s_i terms
  
  R0.d.s1=shed1; R0.i.s1=shed1;
  R0.d.s2=shed1; R0.i.s2=shed1; #R0 values
  R0.d.s3=shed1; R0.i.s3=shed1; #R0 values
  fit.d.s1=shed1;  fit.i.s1=shed1;  
  fit.d.s2=shed1;  fit.i.s2=shed1;  
  fit.d.s3=shed1;  fit.i.s3=shed1;  

  R0within=shed1; #within-host R0, just for diagnostics

  killratio=shed1;  t.final=shed1;

  for (ct1 in 1:length(qvec))  #loop over different values of q (strength of IR)
  {
	  for (strain in 1:12) #loop over all 12 strains
  	{
  	  cw=cwvec[strain]
      R0=p*b*Utc0/(cw*d)
  	  q=qvec[ct1]
  	  
      odeoutput=ode(y=Y0,times=timevec,func=odeequations,parms="",method="lsoda",rootfun=rootfun,atol=1e-5); #run and stop when rootfun happens
      tvec=odeoutput[,1];
      vvec=odeoutput[,4];
      lvvec=odeoutput[,8];
    
      max.ind=which.max(vvec); #find peak
      if (vvec[max.ind]<=1) {cat("infection did not take off!")}
      d.ind=max.ind+min(which(vvec[-(1:max.ind)]<=Vlow)); #find 1st time virus is below 1 after peak
  
      #final time for each integration
      t.final[strain,ct1]=tvec[d.ind]
      
      #percentage of total virus removal due to IR 
      killratio[strain,ct1]=100*odeoutput[d.ind,7]/(odeoutput[d.ind,6]+odeoutput[d.ind,7]);
  
      print(sprintf('q is %f, R0 within-host is %f, final time %f, IR kill percentage %f',q,R0,t.final[strain,ct1],killratio[strain,ct1]));
  
    #different forms for total virus/shedding
    #integral of total virus load
     inf.t=tvec[1:d.ind];
     inf.v=vvec[1:d.ind];
     shed1[strain,ct1]=sum(diff(inf.t)*(inf.v[-length(inf.v)]+inf.v[-1])/2)   
         
     #compute shedding value, i.e. function under integral
     cc1=5; cc2=5; cc3=2.5; 
     shed=inf.v * (cc1*pmax(0,log10(pmax(0,inf.v)))^cc2/(cc3^cc2+pmax(0,log10(pmax(0,inf.v)))^cc2))
     shed2[strain,ct1]=sum(diff(inf.t)*(shed[-length(shed)]+shed[-1])/2) ;
         
     #total log virus load as area under the curve 
     shed=pmax(0,log10(pmax(0,inf.v)));
     shed3[strain,ct1]=sum(diff(inf.t)*(shed[-length(shed)]+shed[-1])/2) ;
    
             
     
  	 #"R0" (not exactly, numerator and denominator for F) for 3 different shedding, direct and indirect transmission
     R0.d.s1[strain,ct1]=shed1[strain,ct1]       
     R0.i.s1[strain,ct1]=shed1[strain,ct1]/cbvec[strain];       
     R0.d.s2[strain,ct1]=shed2[strain,ct1];
     R0.i.s2[strain,ct1]=shed2[strain,ct1]/cbvec[strain];
     R0.d.s3[strain,ct1]=shed3[strain,ct1];
     R0.i.s3[strain,ct1]=shed3[strain,ct1]/cbvec[strain];

  
  } #finish loop over 12 strains   

  #fitness
	fit.d.s1[,ct1]=R0.d.s1[,ct1]/R0.d.s1[1,ct1];
  fit.i.s1[,ct1]=R0.i.s1[,ct1]/R0.i.s1[1,ct1];   
  fit.d.s2[,ct1]=R0.d.s2[,ct1]/R0.d.s2[1,ct1];
  fit.i.s2[,ct1]=R0.i.s2[,ct1]/R0.i.s2[1,ct1];   
  fit.d.s3[,ct1]=R0.d.s3[,ct1]/R0.d.s3[1,ct1];
  fit.i.s3[,ct1]=R0.i.s3[,ct1]/R0.i.s3[1,ct1];   
  
	print(sprintf('for q = %f, min/max kill percentage is %f/%f',q,min(killratio[,ct1]),max(killratio[,ct1])));
	  
  ct1=ct1+1;
} #finish loop over q


###########################
#figure for paper - supplementary material, S2
###########################
#plot ordered by clearance at 40c
cwsort=sort(cwvec,index.return=TRUE)
cw.ind=cwsort$ix
xlabels.raw=c("H1N1","H2N4","H3N2","H4N6","H5N2","H6N4","H7N6","H8N4","H9N2","H10N7","H11N6","H12N5");
xlabels=xlabels.raw[cw.ind] #resort labels according to cw

dummyct=0;

graphics.off(); #close all graphics windows
ww=6.83; wh=9;
windows(width=ww, height=wh)   #open a window for plotting the results
par(mfrow=c(3,2))
mult=1; mult.ax=1; mult.xax=1; mult.lab=0.8; #multiplier for plots
dx=0.0; #spacer between values


#create all figures

for (ct2 in 1:3) #loop over shedding scenarios
{
  if (ct2==1) {fit.d=fit.d.s1; fit.i=fit.i.s1; usecol="black"; lab1="s[1]"; }  
  if (ct2==2) {fit.d=fit.d.s2; fit.i=fit.i.s2; usecol="blue"; lab1="s[2]"; }  
  if (ct2==3) {fit.d=fit.d.s3; fit.i=fit.i.s3; usecol="red"; lab1="s[3]"; }  

  
  filenamepdf.direct=paste("fig-fitness-s",as.character(ct2),"-direct.pdf",sep="");
  filenamepdf.indirect=paste("fig-fitness-s",as.character(ct2),"-indirect.pdf",sep="");
  
  #plot the fitness for different shedding - direct transmission
  par(mar=c(3, 3.5, 2, 0.5)) #bottom, left, top, right margins
  plot((1:12)-dx,fit.d[cw.ind,1],type="p",ylim=c(0.01,2.1),pch=21,col=usecol,log="",xlab="",ylab="",cex.axis=mult.xax,xaxt="n")
  points((1:12),fit.d[cw.ind,2],type="p",pch=22,col=usecol,cex=mult)
  points((1:12)+dx,fit.d[cw.ind,3],type="p",pch=23,col=usecol,cex=mult)
  points((1:12)+dx,fit.d[cw.ind,4],type="p",pch=24,col=usecol,cex=mult)
  points((1:12)+dx,fit.d[cw.ind,5],type="p",pch=25,col=usecol,cex=mult)
  axis(1,at=1:12,lab=xlabels,cex.axis=mult.lab,las=2);
  ltext=c((paste("q=",qvec[1],sep="")),(paste("q=",qvec[2],sep="")),paste("q=",qvec[3],sep=""),paste("q=",qvec[4],sep=""),paste("q=",qvec[5],sep=""))
  if (dummyct==0) {legend("topright",ltext,col=usecol,pch=c(21:25),cex=mult,bty="n"); dummyct=1;}
  if (ct2==1) {mtext(expression(paste(s[1], " - direct transmission")),side=3,line=0)} #y-axis
  if (ct2==2) {mtext(expression(paste(s[2], " - direct transmission")),side=3,line=0)} #y-axis
  if (ct2==3) {mtext(expression(paste(s[3], " - direct transmission")),side=3,line=0)} #y-axis
  mtext("relative fitness",side=2,line=2) #y-axis

  #plot the fitness for different shedding - indirect transmission
  par(mar=c(3, 3, 2, 0.5)) #bottom, left, top, right margins
  plot((1:12)-dx,fit.i[cw.ind,1],type="p",ylim=c(0.04,10.5),pch=21,col=usecol,xlab="",ylab="",log="y",cex.axis=mult.xax,xaxt="n")
  points((1:12),fit.i[cw.ind,2],type="p",pch=22,col=usecol,cex=mult)
  points((1:12)+dx,fit.i[cw.ind,3],type="p",pch=23,col=usecol,cex=mult)
  points((1:12)+dx,fit.i[cw.ind,4],type="p",pch=24,col=usecol,cex=mult)
  points((1:12)+dx,fit.i[cw.ind,5],type="p",pch=25,col=usecol,cex=mult)
  axis(1,at=1:12,lab=xlabels,cex.axis=mult.lab,las=2);
  ltext=c((paste("q=",qvec[1],sep="")),(paste("q=",qvec[2],sep="")),paste("q=",qvec[3],sep=""),paste("q=",qvec[4],sep=""),paste("q=",qvec[5],sep=""))
  if (ct2==1) {mtext(expression(paste(s[1], " - environmental transmission")),side=3,line=0)} #y-axis
  if (ct2==2) {mtext(expression(paste(s[2], " - environmental transmission")),side=3,line=0)} #y-axis
  if (ct2==3) {mtext(expression(paste(s[3], " - environmental transmission")),side=3,line=0)} #y-axis
  
  mtext("relative fitness",side=2,line=2) #y-axis
} #finish plotting loop  

dev.print(device=tiff,filename ="S2.tiff",width=ww, height=wh, units="in",pointsize = 12,compression =c("lzw"), res=300)
 

  
###################################################################
#end main program
###################################################################                       
