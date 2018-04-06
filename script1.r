##################################################################
#script used to fit the virus-decay data from Brown et al 2009
#also does simple stats in the results
##################################################################
rm(list=ls());
graphics.off();
require(xtable)  #use xtable to get format of results useable in latex manuscript file
require(nloptr)  #package for optimization routine


#load Brown et al 2009 data from spreadsheet
#data from Brown et al 2009 is reported as RT values, which is the time (in days) required
#for each virus to reduce infectivity by 90%, as evidenced by
#a decrease in titer by 1 log10 TCID50/ml
#i.e. we have 0.1=V(t=RT)/V(t=0)
AIdata=as.matrix(read.csv('decay-data.csv'))
strain=AIdata[1,2:13] #strain name
group=AIdata[2,2:13]   #genogroup of strain
Rtdata=matrix(as.numeric(AIdata[3:9,2:13]),nrow=7,ncol=12)

#we want decay rate,c, according to exp(-ct) (t in days). 
#we have 0.1=V(t=RT)/V(t=0), with V(t=0)=V0 and V(t=RT)=V0*exp(-c * RT) we get
# 0.1 = exp(-c * RT) or c= -log(0.1)/RT 
cw=-log(0.1)/Rtdata; #decay data to fit
Tvec=as.numeric(AIdata[3:9,1]); #temperature vector

Tplot=seq(1,40,length=100); #temperature vector for plotting purposes
cwvec=matrix(0,nrow=12,ncol=100); #store all best fit plots

#we fit decay rate c as function of temperature for each strain
bestfitmat=matrix(0,12,2); bestfitmat.lin=bestfitmat;
c.at.5=rep(0,12);        #decay at 5 degrees, units of 1/day
c.at.40=rep(0,12);

strainnames=c("H1N1","H2N4","H3N2","H4N6","H5N2","H6N4","H7N6","H8N4","H9N2","H10N7","H11N6","H12N5");

names(c.at.5)<-strainnames
names(c.at.40)<-strainnames

resmatrix=matrix(0,12,5); #combine all results
rownames(resmatrix) <- strainnames
colnames(resmatrix) <- c('group','alpha','gamma','c.at.5','c.at.40')
resmatrix[,1]=1; resmatrix[c(3,4,7,10),1]=2; #assign group type, H3, H4, H7 and H10 belong to group 2, rest to group 1

###################### Function Called For Fitting Data ########################
fitfunction=function(pars)
  { 
    a=pars[1]; g=pars[2] 
    cwfit=a*exp(g*(Tvec))
    #Returns the objective value to be minimizied by optim function
    return(sum((cwfit-cw[,x])^2))
  }

############################# Main Program ####################################
for (x in 1:12) #repeat for each of the 12 strains
{   
   a0=1e-2;  g0=0.05;  
   p0=c(a=a0,g=g0); 
   lb=c(1e-12,1e-12); ub=c(1e5,1e5); 

   xerr.local=1e-12;    
   stps=10000;
   #call optim to optimize the proposed function for fitting
   fitresult <- nloptr(x0=p0,eval_f=fitfunction,lb=lb,ub=ub,opts=list("algorithm"="NLOPT_LN_PRAXIS",xtol_rel=xerr.local,maxeval=stps,print_level=0)); 
   finalparams=fitresult$solution
   bestfitmat[x,]=finalparams 
    
   #produce best fit curve for each strain
   cwvec[x,]=finalparams[1]*exp(finalparams[2]*Tplot)                
   c.at.5[x]=finalparams[1]*exp(finalparams[2]*5)
   c.at.40[x]=finalparams[1]*exp(finalparams[2]*40)
   resmatrix[x,2:5]=c(finalparams[1],finalparams[2],c.at.5[x],c.at.40[x])
   
}

##############################################################################################################
#plot data and clearance as function of temperature for all strains over the whole temperature range
#all 12 strains in one figure
cxl=1; cxa=1;
graphics.off(); #close all graphics windows
ww=6.83; wh=3;
windows(width=ww, height=wh)   #open a window for plotting the results
par(mfrow=c(2,6))
par(mar=c(2, 2, 0.5, 0.5),oma=c(3,3,0.5,0.5)) #bottom, left, top, right margins
for (j in 1:12)
{
      plot(Tplot,cwvec[j,],type="l",ylab="",xlab="",xlim=c(2,40),ylim=c(0,5.2),cex.axis=cxa,cex.lab=cxl)
      points(Tvec,cw[,j],pch=j)
      text(20,5,labels=strainnames[j],cex=1.5)
     
}
mtext("Temperature [C]",side=1,line=1,outer=TRUE) #x-axis - change line to move in/out
mtext("Virus Decay Rate [1/day]",side=2,line=1,outer=TRUE) #y-axis

dev.print(device=tiff,filename ="fig3.tiff",width=ww, height=wh, units="in",pointsize = 12,compression =c("lzw"), res=300)


###############################################################
#compute correlation between ranks of alpha and gamma (using spearman rank correlation coefficient)
ag.cor=cor(resmatrix[,2],resmatrix[,3],method="spearman") #get correlation between alpha and gamma
ag.pval=cor.test(resmatrix[,2],resmatrix[,3],method="spearman",exact=TRUE) #get p-value for correlation between alpha and gamma
print(sprintf('correlation is %f, p-value is %e',ag.cor,ag.pval$p.value))

###############################################################
#compute regression between log of alpha and log of gamma
reg=lm(log(resmatrix[,3])~log(resmatrix[,2]))
eta=reg$coefficients[1]
kappa=reg$coefficients[2]
summary(reg) #write result to screen


####################
#3-panel Figure for manuscript
####################


graphics.off(); #close all graphics windows
ww=6.83; wh=2.5;
windows(width=ww, height=wh)   #open a window for plotting the results
cxl=1; cxa=1;
par(mfrow=c(1,3))
par(mar=c(3.5, 3.5, 0.5, 0.5)) #bottom, left, top, right margins

strainvec=c(8,9,10)

dct=1; #dummy counter
cvec=c("black","black","black")
xmax=40; ymax=2;
for (j in strainvec)
{
      if (dct==1)
      {
      plot(Tplot,cwvec[j,],type="l",ylab="",xlab="",log="y",xlim=c(2,xmax),ylim=c(0.01,ymax),cex.axis=cxa,cex.lab=cxl,lty=dct,col=cvec[dct],lwd=2)

      }
      if (dct>1)
      {
        lines(Tplot,cwvec[j,],type="l",lty=dct,col=cvec[dct],lwd=2)
      }  
      points(Tvec,cw[,j],pch=j,col=cvec[dct])
      dct=dct+1;
}
mtext("Temperature [C]",side=1,line=2.5) #x-axis - change line to move in/out
mtext("Virus Decay Rate [1/day]",side=2,line=2) #y-axis
text(3,1.8,"A)",font=2)
legend("bottomright",strainnames[strainvec],lty=c(1,2,3),lwd=2,bty="n",col=cvec)


                                               
#plot alpha and gamma (on log scales) with regression line
par(mar=c(3.5, 3.5, 0.5, 0.5)) #bottom, left, top, right margins
mult=1; mult.ax=1; mult.lab=1; #multiplier for plots
xmin=6e-4; xmax=0.15; ymin=6.5e-2; ymax=0.235;
mult=1;
plot(resmatrix[,2],resmatrix[,3],xlab="",ylab="",xlim=c(xmin,xmax),ylim=c(ymin,ymax),pch=18,col="white",log="xy",cex=mult,cex.axis=mult.ax,cex.lab=mult.lab)
text(resmatrix[,2],resmatrix[,3],labels=strainnames,cex=mult)
aa=seq(from=min(resmatrix[,2]),to=max(resmatrix[,2]),length=100)
lines(aa,exp(eta)*aa^kappa,lwd=2)
abline(reg,col="black",lwd=2)
mtext(expression(paste("intercept, ",alpha)),side=1,line=2.5,cex=1) #x-axis - change line to move in/out
mtext(expression(paste("temperature dependence, ",gamma)),side=2,line=2,cex=1) #y-axis
text(xmax-0.02*xmax,ymax-0.02*ymax,"B)",font=2)


#plot ranks for alpha and gamma
par(mar=c(3.5, 3.5, 0.5, 0.5)) #bottom, left, top, right margins
xmin=0.1; ymin=0.2; xmax=13; ymax=12.5;
plot(rank(resmatrix[,2]),rank(resmatrix[,3]),xlim=c(xmin,xmax),ylim=c(ymin,ymax),xlab="",ylab="",pch=18,cex=mult,col="white",cex.axis=mult.ax,cex.lab=mult.lab)
text(rank(resmatrix[,2]),rank(resmatrix[,3]),labels=strainnames,cex=mult)
fitline=lm(rank(resmatrix[,2])~rank(resmatrix[,3]))
abline(fitline,lwd=2)
mtext(expression(paste("intercept, ",alpha,", (rank)")),side=1,line=2.5,cex=1) #x-axis - change line to move in/out
mtext(expression(paste("temp. dep., ",gamma,", (rank)")),side=2,line=2,cex=1) #y-axis
text(xmax-0.02*xmax,ymax-0.02*ymax,"C)",font=2)

 
dev.print(device=tiff,filename ="fig4.tiff",width=ww, height=wh, units="in",pointsize = 12,compression =c("lzw"), res=300)

filename="decay-fitting-results.Rdata"
save(list = ls(all=TRUE), file = filename);

