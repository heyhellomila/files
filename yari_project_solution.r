############################################################
##example solution to the project described in Chapter 6 of the YaRI tutorial
##written by Andreas Handel (ahandel@uga.edu), last change 9/3/12
##this short program is a companion to the YaRI tutorial, 
##available at http://ahandel.myweb.uga.edu/resources.htm 
############################################################
rm(list=ls()) #this clears the workspace to make sure no leftover variables are floating around. Not strictly needed but often a good idea
graphics.off(); #close all graphics windows, in case there are still some open from previous stuff we did

d=0.1; 					#dispersal rate of individuals
L=20; 					#number of patches
tmax=50; 				#maximum number of timesteps to run simulation
lambdavec=c(rep(0.9,L/2),rep(1.2,L/2)); #create a vector of length L with growth rates for every patch
mvec=matrix(0,L,1); #note that mvec is a vector, but I created it as a matrix with only a single column = a standing vector
njt=matrix(0,tmax,L); #create the matrix/array njt with tmax rows and L columns that will store the population. 
          					#Each row will contain the population at time t for the L patches in the L different columns 
					         #to initialize, set all values to zero.
njt[1,]=5;				#every patch at time t=1 gets 5 individuals


#########################################################
#1st way of doing it - an outer loop over time, and an inner loop over all patches for every time-step
#########################################################
#start the main simulation loop, running from time t=1 to t=tmax
#note that the loop only goes to tmax-1, since we address the matrix njt at njt[t+1,], and the matrix only has the size tmax
for (t in 1:(tmax-1))
{
	for (p in 1:L) #inner loop over patches - growth step 1st
	{
	  mvec[p]=lambdavec[p]*njt[t,p]; 
  }
  for (p in 1:L) #dispersal step next
  {
	  if (p==1) { njt[t+1,p]=(1-d)*mvec[p] + d*mvec[p+1] }     #dispersal for the 2 edges and the middle
	  if (p==L) {  njt[t+1,p]=(1-d)*mvec[p]+d*mvec[p-1]  }
	  if (p>1 & p<L) { njt[t+1,p]=(1-2*d)*mvec[p]+d*mvec[p-1]+d*mvec[p+1]} 
	  
	}
}

plot(1:tmax,log(apply(njt,1,sum))); #plot total population size as a function of time
#since the command sum(njt) would sum not only over the patches at every time, but also over all times and return only one value,
#one can't use it straight. Instead, one can use the command apply, which applies another command (here sum) to every row or column
#see also table 5 in Ellner's tutorial. Also note that it's perfectly ok to stack commands into each other, like plot(log(apply())). 
#But at some point, it can get too messy and it's better to split it and do one command, such as njtsum=apply(njt,1,sum) 
#and then njtlog=log(njtsum), etc.

browser(); #this command stops the execution of the program and shows you the console. To continue the rest of the program, press c.
           #this command is very useful if you want to find mistakes in your code and don't know where it goes wrong. 
           #just place the command in certain places, let the program run, and see if it gets this far. 
           #also, when you get the console, you can check your variables to see if they all 'behave'
           #I have used the command here not to find a mistake, but to make the first plot, then stop. 
           #this way, one can look at the plot, and then by pressing 'c' at the console, continue to the next plot below
           #of course one could also create two separate plots, or 2 plots in the same window, etc. 
                                                                                                                          
image(1:tmax,1:L,log(njt)) #plot a 2-d image of the population log-densities as a function of patch position and time

browser();
	
#########################################################
#2nd way of doing the same thing - create a dispersal matrix that applies all the 
#dispersal steps at one point, then only do one loop over time
#########################################################
njt2=matrix(0,tmax,L); #create the matrix/array njt with tmax rows and L columns that will store the population. 
njt2[1,]=5;				#every patch at time t=1 gets 5 individuals

#create of matrix that contains dispersal rates/coefficients. Put 1-2d on the diagonal and d on the off-diagonals.
dispmatrix=diag(1-2*d,L,L)+cbind(rep(0,L),diag(d,L,L-1))+rbind(rep(0,L),diag(d,L-1,L))
#we still need to deal with the boundaries and give them special values
dispmatrix[1,1]=1-d; dispmatrix[L,L]=1-d; 
                                                                           
#start the main simulation loop, running from time t=1 to t=tmax
for (t in 1:(tmax-1))
{
	mvec=lambdavec*njt2[t,];		#multiply the population in every patch at time with the growth rate for each patch
	njt2[t+1,]=dispmatrix%*%mvec;	#multipy the vector mvec with the matrix dispmatrix, which contains the rates of dispersal
					#note that this is a matrix-vector multiplication and therefore requires the %*% symbol
          #also note that mvec is a lying vector of dimension 1 x L. so technically, one can't multiply it with an L x L matrix
          #but R is 'smart' and switches mvec into a standing L x 1 vector, then does the matrix multiplication, which results
          #in another L x 1 vector. R switches this result around again so it matches njt[t+1,], which is again a 1 x L vector.
          #Conclusion: R is rather powerful, and it's sometimes convenient. But if you are not sure you are doing the right thing, 
          #stuff can go wrong a lot. So be careful.   
}
plot(1:tmax,log(apply(njt2,1,sum))); #plot total population size as a function of time
browser();
image(1:tmax,1:L,log(njt2)) #plot a 2-d image of the population log-densities as a function of patch position and time
browser();
