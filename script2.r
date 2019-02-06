##################################################################################
#fitting duck flu data from Justin
##################################################################################
rm(list=ls()) #this clears the workspace to make sure no leftover variables are floating around. Not strictly needed
graphics.off(); #close all graphics windows
require(deSolve)  #loads ODE solver package
require(nloptr) #for fitting

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#load experimental data from Justin
#for details, see their paper
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
duckdata=read.csv("duck-data.csv",header=TRUE)

#use fecal/excrement measurements
#data is already in log10 EID units
fecal.data.full=duckdata[duckdata$Sample=="fecal",] 

#remove entries where no sample was collected or no feces could be obtained (NA entries)
fecal.data=fecal.data.full[!is.na(fecal.data.full$Titer),]

#convert Day & AM/PM information (coded as 1/2) into hours post infection. T=0 is at day 0 at 4pm.
time.data=24*fecal.data$Day+9+7*((fecal.data$Time+1) %% 2); 
virus.data=fecal.data[,5]

#any 0 in the dataset is no virus isolation, a 1 is virus isolated but below limit of titration/quantification
#set both to limit of detection and in fitting, treat as censored
Vmin=0.6; #level of detection for virus, log10 EID units 
virus.data[virus.data==0]=Vmin; 
virus.data[virus.data==1]=Vmin;


###################################################################
#function specifying the ode model
###################################################################
odeequations=function(t,y,x) 
{ 

    
      b=exp(x[1]); p=exp(x[2]); d=exp(x[3]); clear=cw; #clearance as mean of other strains 
   
    
    #these are the differential equations
		dUtc=-b*y[3]*y[1];
		dItc=b*y[3]*y[1]-d*y[2]; 
    dV=p*y[2]-clear*y[3]; #virus 
   
    dydt=c(dUtc,dItc,dV)
                                         
    return(list(dydt))
 } #end function withinode
    

###################################################################
#function that fits the ODE model to data 
###################################################################
fitfunction <- function(pars.all)
{
        
        #integrate ODE 
        tvec=seq(0,480,by=1);
        Y0=c(U0, I0, V0); pars.ode=pars.all

        odesol=try(lsoda(Y0,tvec,odeequations,pars.ode,atol=atolv,rtol=rtolv)); #
        #try command catches error from lsoda. If error occurs and lsoda "breaks", we exit the whole optimizer routine with a high objective function value
        if (length(odesol)==1) {cat('!!unresolvable lsoda error - triggering early return from optimizer!!'); return(1e10) }
                
        #extract model values used for fitting
        vir.model.lin=odesol[match(time.data,odesol[,1]),4]; #extract values for virus load at time points corresponding to experimental measurements
        vir.model=log10(pmax(eps,vir.model.lin)); #convert to log scale. prevent potential negative entries by setting them to a small number (eps). 
        
        #for censored data, set difference to zero if model prediction is lower than censored value/level of detection
        vir.model.cens=vir.model
        vir.model.cens[(vir.data==Vmin & (vir.data-vir.model)>0)]=Vmin
        
        #fit in SSR approach
        Fobject=sum((vir.model.cens-vir.data)^2); 
        if (is.na(Fobject)) {Fobject=1e10}; #if something went wrong, e.g. because ODE solver messed up and returned nonsense leading to NA, we set Fobject to a really large number, otherwise optimizer might stop
               
         
       	return(Fobject)
} #end of fitting function, main program starts below


#################################
#main program
#################################
tstart=proc.time(); #capture current time to measure speed

eps=1e-12;
atolv=1e-6; rtolv=1e-6; #tolerances for ODE solver

finalplot=1;

U0=2.5E7; #number of target cells in duck - (Uni et al 1998)
I0=0;

cw=2.78/24; #mean clearance rate per hour of other strains at 40C
              
#bounds on parameters, all is in units of hours
blow=1e-12; bhigh=1e6;
plow=1e-8; phigh=1e5;
dlow=1e-2; dhigh=1e3;  #duration of infectious period in 1/hours 

clearlow=1e-2; clearhigh=1000;
V0low=1e-12; V0high=1e7;

vir.data=virus.data; 

    
     V0=1;
     p0=c(9.253823e-08, 4.911733e-02, 7.941870e-02); Obj.Fct=3.504143e+02; AICc=1.359814e+02; 
     lb=c(blow,plow,dlow); ub=c(bhigh,phigh,dhigh); #units are per hour
    
    p0=log(p0); lb=log(lb); ub=log(ub);
   
    
    stps=5000; xrel=1e-12; mtime=0.5*60*60; #hours runtime
    xerr.global=1e-8;
    xerr.local=1e-12;
     
    
    fres <- nloptr(x0=p0,eval_f=fitfunction,lb=lb,ub=ub,opts=list("algorithm"="NLOPT_LN_SBPLX",xtol_rel=xerr.local,maxeval=stps,maxtime=mtime,print_level=0)); 
    
    fpar=exp(fres$solution); Obj.final=fres$objective; Obj.improve=(Obj.Fct-Obj.final)/Obj.final*100; 
          
    tend=proc.time(); #capture current time
    tdiff=tend-tstart;
    print(sprintf('Optimization took %f seconds',tdiff[3]));
    #code that checks and warns the user if the fit has not converged and instead stopped for another reason (e.g. because it reached maxsteps iterations)
    print(sprintf('Solver status is %d (between 1-4 is good, see nlopt webpage for details)',fres$status)); 

        
        
    R0=fpar[1]*fpar[2]*U0/(fpar[3]*cw)
    cat(sprintf('V0=%e, b=%e, p=%e, d=%e, c=%e, R0=%e; Obj.Fct=%e; Obj.Improve=%e; AICc=%e;\n',V0,fpar[1]*24,fpar[2]*24,fpar[3]*24,cw*24,R0,Obj.final,Obj.improve,AICc));
    
    print(sprintf('parameters are in units of 1/day'));

    #print final/best fit plot into a file
    if (finalplot==1)
    {
      ww=3.27; wh=ww;
      graphics.off()
      windows(width=ww,height=wh)
      par(mar=c(3, 3, 0.5, 0.5)) #bottom, left, top, right margins
      
        tvec=seq(0,480,by=1); 
        if (scenario<3) 
        {
          Y0=c(U0, I0, fpar[1]); #fit initial virus load 
          pars.ode=log(fpar[-1]);
        }
        if (scenario==3) 
        {
          Y0=c(U0, I0, V0); #fit initial virus load 
          pars.ode=log(fpar);
        }
        odesol=lsoda(Y0,tvec,odeequations,pars.ode,atol=atolv,rtol=rtolv); #
        mult=1; mult.ax=1; mult.lab=1; #multiplier for plots
        
        plot(time.data,10^vir.data,col='black',type='p',ylab="",xlab="",log="y",xlim=c(0,360),ylim=c(1e0,1e9),pch=19,cex.axis=mult.ax,cex.lab=mult.lab);  
        lines(odesol[,1],odesol[,4],col='black',lwd=2)  #virus
        lines(odesol[,1],rep(10^Vmin,length(odesol[,1])),col='black',lwd=2,lty=2)  #level of detection
        legend('topright',c("data","model"),col=c('black','black'),lwd=2,lty=c(-1,1),pch=c(19,-1),merge=TRUE,cex=mult,bty="n")
        mtext("hours post infection",side=1,line=2) #x-axis - change line to move in/out
        mtext("virus load [EID50/g]",side=2,line=2) #y-axis
      
      dev.print(device=tiff,filename ="fig5.tiff",width=ww, height=wh, units="in",pointsize = 12,compression =c("lzw"), res=300)
 
    }


###################################################################
#end main program
###################################################################
