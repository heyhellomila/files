############################################################
##a simple model for an acute virus infection
##written by Andreas Handel (ahandel@uga.edu), last change 5/29/10
##this short program is a companion to the YaRI tutorial, 
##available at http://ahandel.myweb.uga.edu/resources.htm 
############################################################
rm(list=ls()) #this clears the workspace to make sure no leftover variables are floating around. Not strictly needed but often a good idea
graphics.off(); #close all graphics windows, in case there are still some open from previous stuff we did
library(deSolve)  #loads ODE solver package. You need to have the package installed, for instance by going to "Packages" -> "Install Packages"  

###################################################################
#first, we specify the function that describes the differential equation model for the simulated virus infection 
#this function is called by lsoda (the ode solver) in the main program
###################################################################
odeequations=function(t,y,parameters) #the function needs to have a certain form, dictated by lsoda
{ 
  	Utc=y[1]; Itc=y[2]; Vir=y[3];  #uninfected target cells, infected target cells, virus
		b=parameters[1]; delta=parameters[2]; p=parameters[3]; clear=parameters[4]; #four model parameters, passed into function by main program
	 
	  #these are the 3 differential equations which describe an acute viral infection
		#see for instance: Perelson (2002) Nature Reviews Immunology for a brief, easy introduction to such models
    dUtcdt=-b*Vir*Utc;
		dItcdt=b*Vir*Utc-delta*Itc;
		dVirdt=p*Itc-clear*Vir;

		return(list(c(dUtcdt,dItcdt,dVirdt))); #this command returns the result, which is the right side of the ODEs as a list, back to the solver (i.e. lsoda). 
} #end function specifying the ODEs

                                                                                                         
###################################################################
#main program
###################################################################

  timevec=seq(0,10,by=0.1); #this creates a vector of times for which integration is evaluated (from 0 to 10 days in steps of 0.1)
	
  #assign numerical values to the parameters and initial conditions of the model  

  #initial conditions             
  Utc0=1e8; #initial number of uninfected cells 
	Itc0=0; #initial number of infected cells
	Vir0=10; #initial condition for free virus V
  Y0=c(Utc0, Itc0, Vir0);  #combine initial conditions into a vector 
  
  #values for model parameters, units are assumed to be 1/days
	b=1e-8; #rate at which virus infects uninfected cells
	delta=2; #rate at which infected cells die
	p=1e2; #rate at which infected cells produce new virus
	clear=10; #rate at which free virus is cleared 
  #Note: we don't call the parameter "c" since c is a reserved function in R that creates vectors (see next line for a use of it)
	#always make sure that your variable names do not clash with built-in definitions, or you might get unpredictable results
  parvec=c(b,delta,p,clear); #this combines all parameters into a vector called parvec which is sent to the ODE function
		
	#call ode-solver to integrate ODEs. this solver is part of the deSolve package. 
  #It starts the simulation with initial conditions Y0, integrates and returns values for the time and the different variables (here: uninfected cells, infected cells, virus)
  #the function odeequations specifies the actual set of differential equations describing our model - see below. 
  #parvec is the vector of parameters we created above. It is sent into the function odeequations and used there - see below
  odeoutput=lsoda(Y0,timevec,odeequations,parvec);

  #plot results
  #odeoutput, which is returned from calling lsoda, contains a matrix with rows that contain values for a given time, and columns with the 1st column being the time and the remaining columns the variables
  #here, column 2 is uninfected cells, column 3 contains infected cell counts, and column 4 contains the virus load
  plot(odeoutput[,1],odeoutput[,2],type="l",xlab="time (days)",ylab="",col="green",lwd=2,log="y",xlim=c(0,10),ylim=c(100,1e9))
  lines(odeoutput[,1],odeoutput[,3],type="l",col="red",lwd=2)
  lines(odeoutput[,1],odeoutput[,4],type="l",col="blue",lwd=2)
  legend(7,1e8, c("Uninfected cells","Infected cells","Virus"),col = c("green","red","blue"),lwd=2)
  
###################################################################
#end main program
###################################################################                           


