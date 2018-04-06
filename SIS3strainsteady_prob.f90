!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Gillespie algorithm to simulate dynamics of drug resistance emergence through compensatory mutations
!produces data for probability of emergence 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

PROGRAM main
	USE ifport !to use sort ifort RNG
	IMPLICIT NONE
	INTEGER*8 :: nmax, maxsteps 
	INTEGER*8 :: dummy,count1,m,j 
	INTEGER*4 :: randseed
	REAL*8 :: t_ini, t_end
	REAL*8 :: S,Iu,It,I1,I2,I3,N
	REAL*8 :: lam,b, minTinvade, maxstep 
    	REAL*8 :: nuu,nut,nu1,nu2,nu3
	REAL*8 :: f,fstart,fend,fstep
	REAL*8 :: mut,mu1,mu2,mutend,mutstart
	REAL*8 :: a(8)
	REAL*8 :: av, invct, invmax
	REAL*8 :: a0, asum, mu_tstart,comptime,invade,extinct
	REAL*8 :: r2a0,tau,r1,r2,r3,a12,t,threshold, thlevel
	REAL*8 :: b1,b2,b3,f0,PN1,PN10, cu, ct, c1, c2, c3 
	REAL*8 :: S0, R0u, R0t, Rwt, R01, R02, R03, R1hat, R2hat, R3hat, Seq, Iueq, Iteq
	REAL*8 :: lamnu1, lamnu2, lamnu3, a1a, a1b, a2a, a3a
	CHARACTER :: resultfile*50, junkinput	

	!nut=12.0; fend=0.06; fstart=0.6; mut=1.0/10.0; R0u=2.0; R01=0.75*R0u; R02=0.85*R0u; R03=0.95*R0u; resultfile='figstoch1steadyprob.txt'; invmax=1E7;
	nut=12.0; fend=0.06; fstart=0.6; mut=1.0/1000.0; R0u=2.0; R01=0.75*R0u; R02=0.85*R0u; R03=0.95*R0u; resultfile='figstoch2steadyprob.txt'; invmax=1E9;
	!nut=12.0; fend=0.12; fstart=0.6; mut=1.0/10.0; R0u=2.0; R01=0.6*R0u; R02=0.75*R0u; R03=0.9*R0u; resultfile='figstoch3steadyprob.txt'; invmax=1E7;
	!nut=12.0; fend=0.12; fstart=0.6; mut=1.0/1000.0; R0u=2.0; R01=0.6*R0u; R02=0.75*R0u; R03=0.9*R0u; resultfile='figstoch4steadyprob.txt'; invmax=1E9;
	

	randseed=123;
	nmax=1E16;
	lam=1.0/25.0; 	
	mu1=mut; 
	mu2=mut;
	N=10000.0; thlevel=0.05;

	nuu=2.0; nu1=2.0; nu2=2.0; nu3=2.0;
	maxsteps=20;
	b=R0u*(lam+nuu);
	b1=R01*(lam+nu1);
	b2=R02*(lam+nu2);
	b3=R03*(lam+nu3);
	cu=lam+nuu; ct=lam+nut+mut; c1=lam+nu1+mu1; c2=lam+nu2+mu2; c3=lam+nu3;

	!mutstart=-2.0; maxsteps=20; mutend=-6.0; f=0.4
	
	CALL seed(randseed) !ifort RNG
        CALL CPU_TIME(t_ini)   !some function specific to ifort to find out time the simulation took
 	 
	DO count1=1,maxsteps !loop over treatment level f or mu
	
	f=fstart-(count1-1)*(fstart-fend)/(maxsteps-1);
	!mu=10**(mutstart-(count1-1)*(mutstart-mutend)/(maxsteps-1));
	
   	invade=0.0 
	extinct=0.0 
	maxstep=0.0
	invct=1.0
	av=1.0		
   	S0=N;
  	Seq=floor(N*ct*cu/(b*(ct-f*ct+f*cu))); 
  	Iueq=floor((1-f)*N*lam*ct*(ct*b-ct*cu-b*ct*f+b*cu*f)/(b*(ct*f-ct-cu*f)*(ct*f*lam-ct*lam-cu*f*lam-cu*f*mut)));
	Iteq=floor(f*N*lam*cu*(ct*b-ct*cu-b*ct*f+b*cu*f)/(b*(ct*f-ct-cu*f)*(ct*f*lam-ct*lam-cu*f*lam-cu*f*mut)));

  	Rwt=b/N*S0*((1-f)/cu+f/ct);

	   Print '(A,7F10.2)', 'f, mut, S, Iu, It, R01, Rwt', f, mut, Seq, Iueq, Iteq, R01, Rwt
	   !Print '(A,3F10.2)', 'br, b2, b3', br, b2, b3

	   DO WHILE (av.LE.invmax)
       		t=0.0
		S=Seq; Iu=Iueq; It=1; I1=0.0; I2=0.0; I3=0.0
	    	lamnu1=lam+nu1; lamnu2=lam+nu2; lamnu3=lam+nu3; 
		a1a=mut*It; a1b=b1/N*S; a2a=b2/N*S; a3a=b3/N*S;	
		threshold=ceiling(thlevel*(Iueq+Iteq));
		DO dummy=1,nmax
		!if (f.le.1) then
		!	Rwt=S*((1-f)*bu/cu+f*bt/ct);
			!Print '(A,F8.3,E12.4,6F12.1)', 'current values:', t, tau, S, Iu, It, Ir, I2, I3 
		!end if
			a(1)=a1a+a1b*I1;
			a(2)=lamnu1*I1;
			a(3)=mu1*I1;

			a(4)=a2a*I2;
			a(5)=lamnu2*I2;
			a(6)=mu2*I2;

			a(7)=a3a*I3;
			a(8)=lamnu3*I3;

			a0=a(1)+a(2)+a(3)+a(4)+a(5)+a(6)+a(7)+a(8)
			tau=LOG(1.0/random(0))/a0 !ifort RNG
			r2a0=random(0)*a0  !ifort RNG			
			asum=0.0
			DO m=1,8
			   asum=asum+a(m)
			   IF (asum.ge.r2a0) THEN
				EXIT
    			   END IF
			ENDDO
        	        select case(m)
				case(1)
					I1=I1+1 
				case(2)
				 	I1=max(I1-1,0.0);  
				case(3)
				 	I1=max(I1-1,0.0); 
				        I2=I2+1 
				case(4)
					I2=I2+1 
				case(5)
				 	I2=max(I2-1,0.0);  
				case(6)
				 	I2=max(I2-1,0.0); 
				        I3=I3+1 
				case(7)
					I3=I3+1 
				case(8)
				    I3=max(I3-1,0.0);
		        end select        	
			t=t+tau
			!mc(m)=mc(m)+1
			if  (((I1+I2+I3).eq.0)) then	
				extinct=extinct+1
				!Tinvade(av)=0.0
				!Print '(A,2F15.7)', 'complete extinction # at time ', extinct, t
				av=av+1;
				EXIT
			end if
			!Iusum=Iusum+tau*Iu
			!Itsum=Itsum+tau*It
		   	if ((I1.gt.threshold).OR.(I2.gt.threshold).OR.(I3.gt.threshold)) then
				invade=invade+1
				!Tinvade(invct)=t
				!minTinvade1=MIN(minTinvade1,t)		
				!Print '(A,1F12.7,I7)', 'invasion at time t, run #', t, av
				av=av+1;
				!invct=invct+1;
				EXIT
			end if			
		ENDDO !finish one simulation
		if (dummy.eq.nmax) then
			Print *,'INCREASE MAXIMUM SIMULATION STEPS' 
		end if
	  ENDDO	!finish averaging
	  !print '(A,15I11)', 'reaction numbers', (mc(j),j=1,15)
	  !print '(A)', 'test'
	  PN1=1.0-(1.0-invade/av)**(mut*Iteq);
	  PN10=1.0-(1.0-invade/av)**(10.0*mut*Iteq);
	  if (count1.eq.1) then 
	  	open(unit=11,file=resultfile,status='replace',access='sequential')
	  else
		open(unit=11,file=resultfile,status='old',access='sequential',POSITION='append')
	  end if
	  Write (11,'(6F22.15)',ADVANCE='YES') f,mut,extinct/av,invade/av,PN1,PN10
	  close (unit=11,status='keep') 
	  f=f+fstep
      ENDDO 	!finish loop over f 
      CALL CPU_TIME(t_end)      
      Print *, 'done in', (t_end-t_ini)/60.0, 'minutes'	
END PROGRAM main