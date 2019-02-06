function dummy=populationevolution
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%This code simulates a population of bacteria that undergo exponential growth and serial dilution cycles.
%The bacteria can mutate to increase their fitness and the program tracks the fitness of all the different subpopulations.
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%While I tried to document the code as well as possible, I'm aware that it is rather terse since I make heavy use of Matlab's vectorization abilities.
%I believe the code works as promised and doesn't contain any serious bugs -- nevertheless use at your own risk :)
%Please feel free to contact me (ahandel/at|emory.edu) if certain aspects are unclear.
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%Also note that this code uses (in addition to Matlab and the statistics toolbox) some external functions
%provided in the lightspeed package, see: http://research.microsoft.com/~minka/software/lightspeed/
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

clear all;

%add lightspeed functions to path
addpath /data/current/matlabstuff/lightspeed

appendix='s'; %will be added to filename under which data is saved
randstate=110; %initialize random number generators. gamrnd uses randn and unidrnd uses rand. threfore need to set both seeds.
rand('twister',randstate); randn('state',randstate);
randomseed([1 1 1]); %set seed for lightspeed functions

avmax=100; %number of 'experiments'
L=50; %number of local, possible (beneficial) mutations 
p=1*10^-5; %mutation rate per genome/bacteria per replication. From Perfeito et al, Science 2007 
spmax=5000; %maximum number of serial passages (about log2(dilfac) generations per serial passage)
ms=12; %number of mutation steps to track, for detailed analysis

rfrac=1.0;

changesites=rfrac*L; %number of neighbor fitness sites that are changed when mutation occurs (0=no epistasis, L=full epistasis)

alp=30; %alp*exp(-alp*s), alpha from gerrish98/rozen02, fitness exponentially distributed 
%alp=20; %alp*exp(-alp*s),  
gam_a=1; %shape parameter. setting it to one makes it an exponential distribution
gam_b=1/30; %scale parameter (smaller=steeper), for gam_a=1, gam_b=1/alpha from gerrish98/rozen02 

dilfac=1000; %bottleneck/factor by which pops are allowed to grow before being diluted (1/D)

fit0=1; %fitness value for initial strain
delfrac=0.0; %fraction of deleterious mutations. both + and - mutations are chosen from exp distr, deleterious mutations simply get a sign switch to become negative.

inisize=500; %initial size for fitness neighbourhood array
FA=nan(L+2+ms,inisize); %preallocate matrix of fitness neighborhood

%fitvals=exprnd(1/alp,1,L);
fitvals=gamrnd(gam_a,gam_b,1,L);

if (delfrac>0),
	changevals=(rand(1,L)<delfrac);
	fitvals(changevals)=-fitvals(changevals); %change some fitness values to negative, according to delfrac proportion
end

FA(:,1)=[fit0 0 zeros(1,ms) fit0+fitvals]'; %fitness of current, number of mutation steps that have occured, ms values record these mutation steps to keep track of every clone. the rest contains the fitness neighbourhood (mutants 1 step away from current mutant), drawn from a gamma distribution. the gamma dist. only gives increase in fitness, not total value
ma=1; %number of currently existing fitness neighborhoods

ni=1+ms+1; ti=1+ms+2; fi=1+ms+3; %indices of internal time, numbers, fitness

localpeak=max(FA(ms+3:end,1)); %save fitness of local peak in this variable
localmin=min(FA(ms+3:end,1)); %fitness of local minimum 

%loop over different population sizes. For each population size, avmax experiments are performed. For the PLoS paper, only two population sizes, corresponding to the experiment, are considered. However, for our current work, we consider more.
for ct2=[3:-1:1]
	%define name of file that will contain data, as well as size of SA and FA arrays (specified below)
	if (ct2==1)
		N0=1*10^2; filename=strcat('dataN',num2str(ct2),appendix,'.mat'); sainisize=5*10^2; %p=1e-3;
	elseif (ct2==2)
		N0=1*10^4; filename=strcat('dataN',num2str(ct2),appendix,'.mat'); sainisize=1*10^3; %p=1e-5;
	elseif (ct2==3)
		N0=1*10^6; filename=strcat('dataN',num2str(ct2),appendix,'.mat'); sainisize=5*10^3; %p=1e-7;
	%elseif (ct2==4)
	%	N0=1*10^6; filename=strcat('dataN',num2str(ct2),appendix,'.mat'); sainisize=5*10^3;
	end
	
	Nmax=dilfac*N0; %population size at end of exponential growth phase
    Fav=[];  %initialize array that will contain the main data, the fitness of the population as a function of generation time
    for ct=1:avmax	%loop over number of experiments/runs
		SA=nan(4+ms,sainisize); %initialize array of different clones
        	SA(:,1)=[0 zeros(1,ms) N0 1/fit0 fit0]'; %number of mutations, ms entries containing the consecutive sites of mutations, frequency of clone, internal time (=1/fitness),  fitness
		x=1; %initial number of different clones
		spct=1; %counter for serial passages (start at 1 instead of zero to allow better indexing of arrays)
		clmaxini=0; clmaxf=0; %max number of clones before and after serial passage, just two other variables to track the dynamics 
		if (ct==1), disp(sprintf('randstate %d, L=%d, min/max fitness available %f/%f',randstate,L,localmin,localpeak)); end
		Fmean=SA(fi,1)/fit0; Fmax=SA(fi,1)/fit0; Fmin=SA(fi,1)/fit0; %set arrays for mean/max/min fitness of the population to 1
		while (spct<=spmax)	%simulate the exponential growth and dilution cycle for spmax serial passages
			Ndiv=sum(SA(ni,1:x)); %current total number of bacteria
			[SAvals,SAind]=sort(SA(ti,1:x)'); SAar=[SAvals SAind]; %sort times at which next division happens. Also make column vector and save as array, is faster.
			%Division works such that a clone of fitness 'fit' divides after time 1/fit. Every clone has an internal time attached. Whichever value is the lowest is the clone that divides next. Then this clone's internal clock is increased by 1/fit. Then the clone with the next lowest time is picked to divide, etc. This approach is more realistic than having clones divide in sync and create something like 2.2 offspring for a clone with increase fitness. 
			while (Ndiv<Nmax)	%do exponential growth until maximum size is reached
				mut3=0; cl=SAar(1,2); %get the position of clone that will divide next (first entry=lowest time in SAind)
			    if (SA(ni,cl)+Ndiv>1.1*Nmax); break; end %if the clonal expansion would overshoot Nmax, leave growth loop and do serial passage/dilution
				mut0=randbinom(p,SA(ni,cl)); %draw from a binomial distribution the number of bacteria that mutate. Mutation rate is p. (We assume that only one mutation per individual/genome can appear). mut0 is the number of mutations occuring when offpsring is created, given that probability per sequence is p. Note that this is a lightspeed function. for purposes of speed, it is included in this fuction at the end, but it still calls external lightspeed routines. 
				if (mut0>0) %check if the currently dividing clone will produce any mutants
					nms=SA(1,cl)+1; %keep track of the number of mutations that the newly created mutants will have undergone
					if (nms>ms), disp(sprintf('too many mutation steps')); end %check if more mutations occcur than we keep track. If this is the case, ms should be increased to keep track of an adequate number of steps 
					%keyboard;
					mutsites=unidrnd(L,mut0,1); %for each of the mut mutations, randomly pick the new mutant 
				    unimutsites=unique(mutsites); %find unique new mutants. if for instance mutant 123 (i.e mutsites=123) occurs 10 times, we need to make sure that we don't produce 10 distinct new mutants but instead 10 of mutant 123
				    if (length(unimutsites)>1) %if multiple mutants occur, check how many and keep track of them
		            	            mutf=hist(mutsites,unimutsites); %record how many mutations occured at each unique site mutated
   				    else 
                		    	    mutf=length(mutsites);
		                    end
               			    mutsites=unimutsites; %record location of unique mutations (together with their frequency mutf) for further use
					mut2=length(mutsites); %recompute number of unique mutations that occur 

					mutids=[nms*ones(1,mut2); repmat(SA(2:nms,cl),1,mut2); mutsites']; %mutation history array for all new clones. this allows us to give each clone a distinct identifier, which allows checking if a currently created mutant already exists.
					premut=ismember(mutids',SA(1:nms+1,1:x)','rows'); %check which newly created mutants already exist
					if (sum(premut)>0) %if some of the created mutants already exist, add new ones to existing ones, don't create a new clone
									mutidspre=mutids(:,premut); %index of already existing mutants
									mutfpre=mutf(premut); %frequency of preexisting mutant					
									premut2=ismember(SA(1:nms+1,1:x)',mutidspre','rows'); %get index for SA which correspond to preexisting mutants
									SA(ni,premut2)=SA(ni,premut2)+mutfpre; %increase frequency of preexisting mutants. note that the newly added mutants could have a different time until next division. this is ignored, instead they are assumed to divide at the same time as the currently already existing mutant
									mutsites=mutsites(~premut); %remove those mutation sites that would lead to an already existing mutant
									mutids=mutids(:,~premut); %remove already existing mutants from new mutant array
									mutf=mutf(~premut); %remove already existing mutants from new mutant frequency array
					end
					mut3=size(mutids,2); %recompute number of unique mutations that occur 

					if (mut3>0) %check if novel mutants are still created or if all occuring mutations where already existing mutants

						mn=~sum(abs(FA(3:ms+2,1:ma)-repmat(SA(2:ms+1,cl),1,ma)),1); %check if mutating clone already has a mutant neighborhood. the fitness landscape as described by the mutant neighborhood is initially 'empty' and is created dynamically and randomly whenever a specific clone undergoes mutation for the first time. however, for all consecutive mutations of this clone, as well as for runs at different population sizes and the avmax different realizations, the fitness landscape remains preserved, corresponding to fixed experimental coditions (note that the epistasis and non-epistasis experiments don't have that restriction, these are created in distinct runs of the code, therefore having independent fitness neighborhoods).
						if (sum(mn)>1), disp(sprintf('something is wrong')); keyboard;  %can't have more than one preexisting mutant neighborhood
						elseif (sum(mn)==0), %if no neighborhood exists, create one
							%find index for ancestor of current clone
							mss=SA(1,cl); %number of mutations of current clone
							anc=find(sum(abs(FA(2:mss+1,1:ma)-repmat([mss-1; SA(2:mss,cl)],1,ma)),1)==0); 
							%use fitness neighbourhood from ancestor, choose new values for a certain fraction of sites. for no epistasis, the new mutant will have the same fitness neighborhood as the ancestor. for complete epistasis, the new 1-mutant fitness neighborhood will be chosen completely at random, drawn from the gamma/exponential distribution as was done initially.
							if (length(anc)>1 || sum(anc)==0), disp(sprintf('anc is wrong')); keyboard; end %internal check to make sure the ancestor exists and is unique
							newfit=FA(ms+3:end,anc); %old fitness neighbourhood from ancestor
        				                newfitsites=unidrnd(L,changesites,1); %pick randomly changesites sites that will get fitness changed
							fitdiff=gamrnd(gam_a,gam_b,1,changesites);
							%fitdiff=exprnd(1/alp,1,changesites); %generate new fitness value differentials for a fraction of the neighbourhood
							if (delfrac>0)
								changevals=(rand(1,changesites)<delfrac);
								fitdiff(changevals)=-fitdiff(changevals); %change some fitness value differentials to negative, according to delfrac proportion
							end
							newfit(newfitsites)=fit0+fitdiff; %absolute fitness values for new sites		

							if (size(newfit)<L), disp(sprintf('newfit is wrong')); keyboard; end %another internal check
							FA(:,ma+1)=[SA(fi,cl); SA(1:ms+1,cl); newfit]; ma=ma+1; faind=ma; %new fitness neighbourhood
						elseif (sum(mn)==1), faind=find(mn==1);	end %if mutant neighborhood exists, get index to location 

						newfitness=FA(ms+2+mutsites,faind)'; %get new fitness for each mutant from neighbor mutant array. 
    
                        if (x+mut3>sainisize), disp(sprintf('SA array too small')); keyboard; end %check if preallocated array is too small. if needed increase sainisize (slows down the code).
			SA(:,x+1:x+mut3)=[mutids; zeros(ms-nms,mut3); mutf; (SA(ti,cl)*ones(1,mut3))+1./newfitness; newfitness]; %clone old population. first row contains new sequence identifiers, 2nd row contains number of mutant steps away from wt, 3rd row contains frequency (=1), fourth row contains new internal time clock (old value + 1/new fitness), fifth row contains new fitnesses
                        x=x+mut3; %update counter that keeps track of number of sequences in sequence array 
                        if (x>sainisize), disp(sprintf('array size wrong')); keyboard; end
			%disp(sprintf('%d clones present, created %d new mutants from clone with frequency %d and fitness %f',x,mut3,SA(ni,cl),SA(fi,cl)));
                        %keyboard;
                       end %end mut3>0
                end %end mut>0

  		        Ndiv=Ndiv+SA(ni,cl); %current total number of sequences. has increased by SA(ni,cl) due to division of this clone    
				SA(ni,cl)=2*SA(ni,cl)-mut0; %replicate clone, minus mutants 
				SA(ti,cl)=SA(ti,cl)+1/SA(fi,cl); %increase internal clock of clone cl that just divided by 1/fitness
		
				if (mut3>0) %check if unique mutations occured. If yes, we need to resort the time array that shows which clone divides next
			                [SAvals,SAind]=sort(SA(ti,1:x)'); SAar=[SAvals SAind]; %if mutation occurs, resort/recreate array
    				else                
					newind=find(SAar(:,1)<SA(ti,cl),1,'last');  %if no mutation occurs, no resorting needed. 
					SAar=[SAar(2:newind,:); [SA(ti,cl) cl]; SAar(newind+1:end,:)]; %put clone that has divided into ordered array, together with its new timestamp
					%keyboard;
				end 

  
            end %end while loop over exponential growth 
            
			%if stationary phase is reached, do serial passage/dilution
			clini(spct)=x; %number of clones present before serial passage
			varini(spct)=var(SA(fi,1:x)/fit0-1); %variance in fitness of clones present before serial passage
			mutmeanini(spct)=mean(SA(1,1:x)); %mean number of mutations undergone by population
			Fmeanini(spct)=(sum((SA(fi,1:x)/fit0).*SA(ni,1:x))/sum(SA(ni,1:x))); %mean population fitness before serial passage
			
			pp=SA(ni,1:x)/sum(SA(ni,1:x));  %probability of each clone to make it through serial passage (to be drawn from a multinomial sampling)
			newclones=sample_hist(pp',N0);  %draw N0 random numbers from a multinomial distribution to find how many clones per type are transferred. Returns a vector of numbers of length size(pp), which is the number of clones present after transfer. This is a lightspeed function. Matlab's internal function is way too slow to be of any use.
			cloc=find(newclones>0); %find clones that get transferred
			SA(ni,cloc)=newclones(cloc); %assign those clones their new frequency
			x=length(cloc); %number of new clones
			SA=[SA(:,cloc) nan(4+ms,sainisize-x)]; %create new array of clones that get transferred, set everything to zero for clones that did not survive the transfer
			%disp(sprintf('Serial passage %d,  current size %e, projected size %e, fitness and clones before: %f, %d after: %d',spct,Ndiv,(2^Fmean(end))*Ndiv,Fmean(end),length(pp),x));

			clfinal(spct)=x; %number of clones present after serial passage
			varfinal(spct)=var(SA(fi,1:x)/fit0-1);; %variance in fitness of clones present after serial passage
			mutmeanfinal(spct)=mean(SA(1,1:x)); %mean number of mutations undergone by population
			Fmeanfinal(spct)=(sum((SA(fi,1:x)/fit0).*SA(ni,1:x))/sum(SA(ni,1:x))); %population fitness after serial passage
			Fmax(spct)=max(SA(fi,1:x)/fit0); %max fitness after serial passage
			Fmin(spct)=min(SA(fi,1:x)/fit0); %min fitness after serial passage
			spct=spct+1; %count number of serial passages
  			
			%serial passage/dilution is over, back to next growth cycle
			%keyboard;			
	end %ends loop over all growth/dilution cycles which constitute one simulation/experiment 

	%keyboard;

	%everything below are variables that keep track of certain aspects of the dynamics, to allow tracking of what happens during the simulation
	%for the PLoS manuscript, only the mean and variance of the Fav array is used.

	Favini(ct,:)=Fmeanini; %for each experiment, save mean fitness before serial passage during experiment
	Fav(ct,:)=Fmeanfinal; %for each experiment, save mean fitness after serial passage during experiment
	Fminar(ct,:)=Fmin; %also save max and min fitness at every serial passage
	Fmaxar(ct,:)=Fmax; %also save max and min fitness at every serial passage
	cliniar(ct,:)=clini; %save number of clones present before every serial passage
	clfinalar(ct,:)=clfinal; %save number of clones present after every serial passage
	variniar(ct,:)=varini; %variance in fitness before every serial passage
	varfinalar(ct,:)=varfinal; %variance in fitness before every serial passage
	mutiniar(ct,:)=mutmeanini; %mean number of mutation steps before serial passage	
	mutfinalar(ct,:)=mutmeanfinal; %mean number of mutation steps after serial passage	

    maxmut(ct)=max(SA(1,:)); %maximum number of mutations that occured
    
    [dummy,cind]=max(SA(ni,:));	%find most frequent clone  
	maxnifreq(ct)=SA(ni,cind)/sum(SA(ni,1:x)); %frequency of most frequent clone
	maxnifit(ct)=SA(fi,cind)/fit0;  %fitness of most frequent clone	
    maxnisteps(ct)=SA(1,cind); %number of mutation steps of most frequent clone
	mn=~sum(abs(FA(3:ms+2,1:ma)-repmat(SA(2:ms+1,cl),1,ma)),1); %check mutant neighborhood of most frequent clone
	mind=find(mn==1);
    %keyboard;
	if (sum(mind)==0), maxnirank(ct)=NaN; %find fitness rank of most frequent clone
	else, maxnirank(ct)=sum(FA(ms+3:end,mind)>SA(fi,cind)); end
    
    [dummy,cind]=max(SA(fi,:));	%find fittest clone
    maxfifreq(ct)=SA(ni,cind)/sum(SA(ni,1:x)); %frequency of fittest clone
	maxfifit(ct)=SA(fi,cind)/fit0;  %fitness of fittest clone	
    maxfisteps(ct)=SA(1,cind); %number of mutation steps of fittest clone
	mn=~sum(abs(FA(3:ms+2,1:ma)-repmat(SA(2:ms+1,cl),1,ma)),1); %check mutant neighborhood of fittest clone
	mind=find(mn==1);
    if (length(mind)>1), keyboard; end
	if (sum(mind)==0), maxfirank(ct)=NaN; %find fitness rank to fittest clone
    else maxfirank(ct)=sum(FA(ms+3:end,mind)>SA(fi,cind)); end
	    

    dvec=datestr(now);	
    disp(sprintf('end run %d, experiment %d/%d, time %s, mean/N1/F1 fitness: %f/%f/%f, N1/F1/max mutations %d/%d/%d, N1/F1 rank %d/%d, max clones before/after passage %d/%d, ma %d',ct2,ct,avmax,dvec(12:end),Fav(ct,end),maxnifit(ct),maxfifit(ct),maxnisteps(ct),maxfisteps(ct),maxmut(ct),maxnirank(ct),maxfirank(ct),max(clini),max(clfinal),ma));
	if (mod(ct,10)==0), save(filename); end %every 10 experiments, save workspace. All important data for each population size will be saved in the file 'filename'. A different matlab script then uses that data to create the figures.
	%profile viewer;
	%keyboard;
	if (ct2==1), sy='b'; else sy='r'; end
	%plot(Fmean,sy); hold on; drawnow;
	%plot(Fmax,'r'); drawnow;
end %ends loop over avmax experiments/simulations
%keyboard;
end %end loop over different population sizes
%profile viewer; keyboard;
dummy=1;

function r = randbinom(p, n)

if n < 15
  % coin flip method this takes O(n) time
  r = 0;
  for i = 1:n
    if rand < p
      r = r + 1;
    end
  end
elseif n*p < 150
  % waiting time method
  % this takes O(np) time
  q = -log(1-p);
  r = n;
  e = -log(rand);
  s = e/r;
  while(s <= q)
    r = r - 1;
    if r == 0
      break
    end
    e = -log(rand);
    s = s + e/r;
  end
  r = n - r;
else
  % recursive method, this makes O(log(log(n))) recursive calls
  i = floor(p*(n+1));
  b = randbeta(i, n+1-i);
  if b <= p
    r = i + randbinom((p-b)/(1-b), n-i);
  else
    r = i - 1 - randbinom((b-p)/b, i-1);
  end
end