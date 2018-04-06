############################################################
##a simple model for a singe outbreak of an infectious disease on a population
##written by Andreas Handel (ahandel@uga.edu), last change 10/6/11
##this short program is a companion to the YaRI tutorial, 
##available at http://ahandel.myweb.uga.edu/resources.htm 
############################################################
rm(list=ls()) #this clears the workspace to make sure no leftover variables are floating around. Not strictly needed but often a good idea
graphics.off(); #close all graphics windows, in case there are still some open from previous stuff we did
library(deSolve)  #loads ODE solver package. You need to have the package installed, for instance by going to "Packages" -> "Install Packages"  

###################################################################
#first, we specify the function that describes the differential equation model for the simulated outbreak 
#this function is called by lsoda (the ode solver) in the main program
###################################################################
odeequations=function(t,y,parameters) #the function needs to have a certain form, dictated by lsoda
{ 
  	S=y[1]; In=y[2]; R=y[3];  #uninfected, infected and recovered people are being simulated. 
    #Note that I'm calling the infecteds neither I nor Inf since both are reserved words in R. I could use I instead if In and R would be ok, but it's cleaner this way  
		b=parameters[1]; d=parameters[2]; 
		
	  #these are the 3 differential equations which describe the outbreak dynamics
		#see for instance: Keeling and Rohani (2007) "Modeling Infectious Diseases in Humans and Animals" for background on those kind of models
    dSdt=-b*S*In;
		dIndt=b*S*In-d*In;
		dRdt=d*In;

		return(list(c(dSdt,dIndt,dRdt))); #this command returns the result, which is the right side of the ODEs as a list, back to the solver (i.e. lsoda). 
} #end function specifying the ODEs

                                                                                                         
###################################################################
#main program
###################################################################

  tmax=100; #maximum time/days for which to run the integration/simulation 
  dt=0.1; #steps for the integration. Note that this only affects the times for which the solution is returned, not the actual steps the ODE solver takes
  timevec=seq(0,tmax,by=dt); #this creates a vector of times for which integration is evaluated (from 0 to tmax days in steps of dt)
	
  #assign numerical values to the parameters and initial conditions of the model  

  #initial conditions             
  S0=1e5; #initial number of susceptible hosts
  In0=1; #initial number of infected hosts
	R0=0; #initial number of recovered hosts
  Y0=c(S0, In0, R0);  #combine initial conditions into a vector 
  
  #values for model parameters, units are assumed to be 1/days
	b=1e-5; #rate at which an infected person infects an uninfected person
  d=0.2; #rate at which a person recovers (1/d is the average duration of the infectious period)
  parvec=c(b,d); #this combines all parameters into a vector called parvec which is sent to the ODE function
		
	#call ode-solver to integrate ODEs. this solver is part of the deSolve package. 
  #It starts the simulation with initial conditions Y0, integrates and returns values for the time and the different variables (here: uninfected/infected/recovered hosts)
  #the function odeequations specifies the actual set of differential equations describing our model - see below. 
  #parvec is the vector of parameters we created above. It is sent into the function odeequations and used there - see below
  odeoutput=lsoda(Y0,timevec,odeequations,parvec);

  #plot results
  #odeoutput, which is returned from calling lsoda, contains a matrix with rows that contain values for a given time, and columns with the 1st column being the time and the remaining columns the variables
  #here, column 2 is uninfected hosts, column 3 contains infected hosts counts, and column 4 contains the recoverd hosts
  #note the log-scale for the plot
  plot(odeoutput[,1],odeoutput[,2],type="l",xlab="time (days)",ylab="",col="green",lwd=2,log="y",xlim=c(0,tmax),ylim=c(1,S0),main="ID outbreak")
  lines(odeoutput[,1],odeoutput[,3],type="l",col="red",lwd=2)
  lines(odeoutput[,1],odeoutput[,4],type="l",col="blue",lwd=2)
  legend("right", c("Uninfected","Infected","Recovered"),col = c("green","red","blue"),lwd=2)
  
###################################################################
#end main program
###################################################################                           


