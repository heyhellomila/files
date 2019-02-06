function [varargout]=compensatorymutation

clear all;

%this loops over all the different situations/parameter settings to create data for all the figures
for ct2=1:1:16

ct1=1; %counter over different values of f

%setting some parameters that don't change
lam=1/25; nuu=2; nu1=2; nu2=2; nu3=2; nu4=2; N=10000; 


%for each scenario, set the parameters accordingly
if (ct2==1)
nut=12; fvec=linspace(0.0,0.6,200); mut=10^(-1); R0u=2; R01=0.75*R0u; R02=0.85*R0u; R03=0.95*R0u; R04=0.1;  mu1=mut; mu2=mut; mu3=0; filename='figdata1.mat';
end
if (ct2==2)
nut=12; fvec=linspace(0.0,0.6,200); mut=10^(-3); R0u=2; R01=0.75*R0u; R02=0.85*R0u; R03=0.95*R0u; R04=0.1; mu1=mut; mu2=mut; mu3=0; filename='figdata2.mat';
end
if (ct2==3)
nut=12; fvec=linspace(0.0,0.6,200); mut=10^(-1); mu1=mut; mu2=mut; R0u=2; R01=0.6*R0u; R02=0.75*R0u; R03=0.9*R0u; R04=0.1; mu3=0; filename='figdata3.mat';
end
if (ct2==4)
nut=12; fvec=linspace(0.0,0.6,200); mut=10^(-3); mu1=mut; mu2=mut; R0u=2; R01=0.6*R0u; R02=0.75*R0u; R03=0.9*R0u; R04=0.1; mu3=0; filename='figdata4.mat';
end
if (ct2==5)
nut=12; fvec=linspace(0.0,0.6,200); mut=10^(-3); mu1=mut; R0u=2; R01=0.6*R0u; R02=0.9*R0u; R03=0.1*R0u; R04=0.1*R0u; mu2=0; mu3=0; filename='figdata5.mat';
end
if (ct2==6)
nut=12; fvec=linspace(0.0,0.6,200); mut=10^(-3); mu1=mut; R0u=2; R01=0.6*R0u; R02=0.75*R0u; R03=0.9*R0u; R04=0.1*R0u; mu2=mu1; mu3=0; filename='figdata6.mat';
end
if (ct2==7)
nut=12; fvec=linspace(0.0,0.6,200); mut=10^(-3); mu1=mut; R0u=2; R01=0.6*R0u; R02=0.7*R0u; R03=0.8*R0u; R04=0.9*R0u; mu2=mu1; mu3=mu1; filename='figdata7.mat';
end
if (ct2==8)
nut=12; fvec=linspace(0.0,0.6,200); mut=10^(-2); mu1=10^(-6); R0u=2; R01=0.6*R0u; R02=0.9*R0u; R03=0.1*R0u; R04=0.1*R0u; mu2=0; mu3=0; filename='figdata8.mat';
end
if (ct2==9)
nut=12; fvec=linspace(0.0,0.6,200); mut=10^(-2); mu1=10^(-3); R0u=2; R01=0.6*R0u; R02=0.75*R0u; R03=0.9*R0u; R04=0.1*R0u; mu2=mu1; mu3=0; filename='figdata9.mat';
end
if (ct2==10)
nut=12; fvec=linspace(0.0,0.6,200); mut=10^(-2); mu1=10^(-2); R0u=2; R01=0.6*R0u; R02=0.7*R0u; R03=0.8*R0u; R04=0.9*R0u; mu2=mu1; mu3=mut; filename='figdata10.mat';
end
if (ct2==11)
nut=12; fvec=linspace(0.0,0.6,200); mut=10^(-1); mu1=mut; R0u=2; R01=0.6*R0u; R02=0.75*R0u; R03=0.9*R0u; R04=0.1*R0u; mu2=mu1; mu3=0; filename='figdata11.mat';
end
if (ct2==12)
nut=12; fvec=linspace(0.0,0.6,200); mut=10^(-2); mu1=mut; R0u=2; R01=0.6*R0u; R02=0.75*R0u; R03=0.9*R0u; R04=0.1*R0u; mu2=mu1; mu3=0; filename='figdata12.mat';
end
if (ct2==13)
nut=12; fvec=linspace(0.0,0.6,200); mut=10^(-3); mu1=mut; R0u=2; R01=0.6*R0u; R02=0.75*R0u; R03=0.9*R0u; R04=0.1*R0u; mu2=mu1; mu3=0; filename='figdata13.mat';
end
if (ct2==14)
nut=12; fvec=linspace(0.0,0.6,200); mut=10^(-2); mu1=10^(-3); R0u=2; R01=0.6*R0u; R02=0.75*R0u; R03=0.9*R0u; R04=0.1*R0u; mu2=mu1; mu3=0; filename='figdata14.mat';
end
if (ct2==15)
nut=12; fvec=linspace(0.0,0.6,200); mut=10^(-4); mu1=10^(-2); R0u=2; R01=0.6*R0u; R02=0.75*R0u; R03=0.9*R0u; R04=0.1*R0u; mu2=mu1; mu3=0; filename='figdata15.mat';
end
if (ct2==16)
nut=12; fvec=linspace(0.0,0.6,200); mut=10^(-3); mu1=10^(-2); R0u=2; R01=0.6*R0u; R02=0.75*R0u; R03=0.9*R0u; R04=0.1*R0u; mu2=10^(-3); mu3=0; filename='figdata16.mat';
end

thlevel=0.05; %threshold level

%main loop over different treatment levels f 
for f=fvec
	%compute all the different quantities as described in the paper
	cu=lam+nuu; ct=lam+nut+mut; c1=lam+nu1+mu1; c2=lam+nu2+mu2; c3=lam+nu3+mu3; c4=lam+nu4;
	b=R0u*cu; b1=R01*c1; b2=R02*c2; b3=R03*c3; b4=R04*c4; 	%choose b1, b2 and b3 such that fitness is a fraction of wt fitness R0u
	R0t=b/ct; 
	S0=N;	
  	S_hatvec(ct1)=floor(N*ct*cu/(b*(ct-f*ct+f*cu))); 
  	Iu_vec(ct1)=floor((1-f)*N*lam*ct*(ct*b-ct*cu-b*ct*f+b*cu*f)/(b*(ct*f-ct-cu*f)*(ct*f*lam-ct*lam-cu*f*lam-cu*f*mut)));
  	It_vec(ct1)=floor(f*N*lam*cu*(ct*b-ct*cu-b*ct*f+b*cu*f)/(b*(ct*f-ct-cu*f)*(ct*f*lam-ct*lam-cu*f*lam-cu*f*mut)));
  	S_hat=S_hatvec(ct1); Iu_hat=Iu_vec(ct1); It_hat=It_vec(ct1);
	Ith=ceil(thlevel*(Iu_hat+It_hat));
  	Rwt0vec(ct1)=b/N*S0*((1-f)/cu+f/ct);
	Rwtvec(ct1)=b/N*S_hat*((1-f)/cu+f/ct);
	Rwt0=Rwt0vec(ct1);
	Rwt=Rwtvec(ct1);
	R1=b1/N*S_hat/c1;
	R2=b2/N*S_hat/c2;
	R3=b3/N*S_hat/c3;
	R4=b4/N*S_hat/c4;
  	R1hatvec(ct1)=R1/Rwt;
  	R2hatvec(ct1)=R2/Rwt;
  	R3hatvec(ct1)=R3/Rwt;
  	R4hatvec(ct1)=R4/Rwt;
	R1hat=R1hatvec(ct1);
  	R2hat=R2hatvec(ct1);
  	R3hat=R3hatvec(ct1);
	R4hat=R4hatvec(ct1);
    	Qvec(ct1)=mut*It_hat;
	r1=(R01/Rwt0-1)*c1;
	r2=(R02/Rwt0-1)*c2;
	r3=(R03/Rwt0-1)*c3;
	r4=(R04/Rwt0-1)*c4;
	Q=Qvec(ct1); Q1=mu1; Q2=mu2; Q3=mu3;	
	I1eq=0; I2eq=0; I3eq=0;
	if (r1<0) I1eq=floor(Q/abs(r1)); end 
	if (r2<0) I2eq=floor(Q*Q1/abs(r1*r2)); end
	if (r3<0) I3eq=floor(Q*Q1*Q2/abs(r1*r2*r3)); end
	disp(sprintf('STARTING: treatment level %f, R1hat %f, R2hat %f, R3hat %f',f, R1hat,R2hat,R3hat));
	%disp(sprintf('Rwt0 %f S %d, Iu %d, It %d, I1eq %d, I2eq %d, I3eq %d, Ith %d',Rwt0,S_hat,Iu_hat,It_hat,I1eq,I2eq,I3eq,Ith));
	%disp(sprintf('R0u %f, R0t %f, R0r %f, R02 %f, R03 %f',R0u,R0t,R01,R02,R03));
	mu2d=mu2/c2; mu1d=mu1/c1; mu3d=mu3/c3;
	%4th strain invasion probability
	aa=R3hat*(mu3d-1); cc=mu3d*R3hat*(1-1/R4hat);
	pint=max((1+aa+cc-sqrt((1+aa)^2-2*(aa-1)*cc+cc^2))/(2*aa),eps); 
	aa=R2hat*(mu2d-1); cc=mu2d*R2hat*pint;
	pint=max((1+aa+cc-sqrt((1+aa)^2-2*(aa-1)*cc+cc^2))/(2*aa),eps); 
	aa=R1hat*(mu1d-1); cc=mu1d*R1hat*pint;
	px4(ct1)=max((1+aa+cc-sqrt((1+aa)^2-2*(aa-1)*cc+cc^2))/(2*aa),eps); 
	%3rd strain invasion probability
	aa=R2hat*(mu2d-1); cc=mu2d*R2hat*(1-1/R3hat);
	%pint=max((1+aa+cc-sqrt((1+aa)^2-2*(aa-1)*cc+cc^2))/(2*aa),eps); 
	pint=(1+aa+cc-sqrt((1+aa)^2-2*(aa-1)*cc+cc^2))/(2*aa); 
	aa=R1hat*(mu1d-1); cc=mu1d*R1hat*pint;
	%px3(ct1)=max((1+aa+cc-sqrt((1+aa)^2-2*(aa-1)*cc+cc^2))/(2*aa),eps); 
	px3(ct1)=(1+aa+cc-sqrt((1+aa)^2-2*(aa-1)*cc+cc^2))/(2*aa); 
	papp3(ct1)=max((1-1/R3hat)*mu1*(R1hat/(1-R1hat))*mu2*(R2hat/(1-R2hat)),eps);
	%2nd strain invasion probability
	aa=R1hat*(mu1d-1); cc=mu1d*R1hat*(1-1/R2hat);
	px2(ct1)=max((1+aa+cc-sqrt((1+aa)^2-2*(aa-1)*cc+cc^2))/(2*aa),eps); 
	papp2(ct1)=max((1-1/R2hat)*mu1*(R1hat/(1-R1hat)),eps);
	%1st srain invasion probability
	px1(ct1)=max(1-1/R1hat,0); papp1(ct1)=max(1-1/R1hat,eps);  
	%px3 gives the total probability for invasion, no matter how many strains are used
	P1S1(ct1)=max([1-(1-px1(ct1))^(Q*1) , 1e-12]); 
	P1S2(ct1)=max([1-(1-px2(ct1))^(Q*1) , 1e-12]); 
	P1S3(ct1)=max([1-(1-px3(ct1))^(Q*1) , 1e-12]); 
	P0=0.95;
	%T1(ct1)=(1/Q)*(log(1-P0)/log(1-px3(ct1)));	
	Ts2(ct1)=1/max(Q*px3(ct1),eps);	
	Ts2a(ct1)=1/max(Q*px2(ct1),eps);	
	Ts2b(ct1)=1/max(Q*px1(ct1),eps);	
		

	tfinaldet=10000/(f+0.1);
	tspan=[0 tfinaldet];
	%if (I1eq>Ith || I2eq>Ith || I3eq>Ith) disp(sprintf('equilibrium too large')); keyboard; end
    	Y0=[S_hat Iu_hat It_hat 0 0 0 0];
	options = odeset('OutputFcn',@myoutputfcn2,'RelTol',1e-4,'AbsTol',1e-6);
	[tvec,yvec]=ode15s(@sirode,tspan,Y0,options,lam,b,b1,b2,b3,b4,f,nuu,nut,nu1,nu2,nu3,nu4,mut,mu1,mu2,mu3,N,thlevel); %integrate ODE
	tfvec(ct1)=tvec(end);
	if (tvec(end)==tfinaldet) tfvec(ct1)=1e50; end
	
	%solution for deterministic T by root finding
	fct1=@(t)(Q/r1*(exp(r1*t)-1)-Ith);
	fct2=@(t)(Q1*Q/(r1*r2*(r1-r2))*(r1*(1-exp(r2*t))-r2*(1-exp(r1*t)))-Ith);
	fct3=@(t)(Q*Q1*Q2*((exp(r1*t)-1)*r2*(r2-r3)*r3+r1^2*((exp(r3*t)-1)*r2+r3-exp(r2*t)*r3)+r1*((1-exp(r3*t))*r2^2+(exp(r2*t)-1)*r3^2))/(r1*(r1-r2)*r2*(r1-r3)*(r2-r3)*r3)-Ith);
	tinv1(ct1)=1e50; tinv2(ct1)=1e50; tinv3(ct1)=1e50;
	if (R1hat>1) 
		tm=tfvec(ct1); while (fct1(tm)<0) tm=tm+10; end
		tinv1(ct1)=fzero(fct1,[0 tm]); 
	end
	if (R2hat>1) 		
		tm=tfvec(ct1); while (fct2(tm)<0) tm=tm+10; end
		tinv2(ct1)=fzero(fct2,[0 tm]); 
	end
	if (R3hat>1) 
		tm=tfvec(ct1)-100; while (fct3(tm)<0) tm=tm+10; end
		tinv3(ct1)=fzero(fct3,[0 tm]); 
	end
	Tinvdet(ct1)=min([tinv1(ct1) tinv2(ct1) tinv3(ct1)]);
	T2(ct1)=Ts2(ct1)+tvec(end);  	
	T2a(ct1)=Ts2a(ct1)+tvec(end);  	
	T2b(ct1)=Ts2b(ct1)+tvec(end);  	
	T4(ct1)=1/max(Q*px4(ct1),eps)+tvec(end);	
	ct1=ct1+1;
end %finish loop over f	
%keyboard;
save(filename);
disp(sprintf('DONE INITIAL CONDITION %d',ct2));
end %finish loop over ct2 (different ini. conditions)

%keyboard;
varargout{1}=ct;

%this function provides the RHS of the ODE 
function dydt=sirode(t,y,lam,b,b1,b2,b3,b4,f,nuu,nut,nu1,nu2,nu3,nu4,mut,mu1,mu2,mu3,N,thlevel)
%three resistant strains
S=y(1); Iu=y(2); It=y(3); I1=y(4); I2=y(5); I3=y(6); I4=y(7);
dydt=[lam*N-lam*S-S/N*(b*Iu+b*It+b1*I1+b2*I2+b3*I3+b4*I4)+nuu*Iu+nut*It+nu1*I1+nu2*I2+nu3*I3+nu4*I4; 
      (1-f)*S*b/N*(Iu+It)-(nuu+lam)*Iu;
      f*S*b/N*(Iu+It)-(nut+lam+mut)*It; 
      mut*It+(b1*S/N-nu1-lam-mu1)*I1;
      mu1*I1+(b2*S/N-nu2-lam-mu2)*I2;
      mu2*I2+(b3*S/N-nu3-lam-mu3)*I3;
      mu3*I3+(b4*S/N-nu4-lam)*I4];

% stop integration when y4 or y5 or y6 become too large 
function status = myoutputfcn2(t,y,flag,lam,b,b1,b2,b3,b4,f,nuu,nut,nu1,nu2,nu3,nu4,mut,mu1,mu2,mu3,N,thlevel)
status=0;
if (isempty(flag)) 
	if (y(4)>thlevel*(y(2)+y(3)) || y(5)>thlevel*(y(2)+y(3)) || y(6)>thlevel*(y(2)+y(3)) || y(7)>thlevel*(y(2)+y(3))) 
	status=1; 
	end; 
end;