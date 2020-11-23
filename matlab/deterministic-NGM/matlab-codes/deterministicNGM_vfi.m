%Deterministic Neoclassical Growth Model - Value function iteration        %
%          Program written by: Alessandro Ruggieri - UAB and BGSE         %
%                          Version 11/05/2016                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear
close all
fprintf('\nSolve Deterministic Neoclassical Growth Model using VFI\n');

%% Parametrization
par.alfa=0.33;                                                            %Capital share of income
par.sigma=1;                                                               %Relative risk aversion
par.delta=0.1;                                                             %Depreciation rate
par.beta=0.95;                                                             %Discount factor
par.n=1000;                                                                %Capital grid points
par.toler=1e-10;                                                           %Tolerance level
par.T=300;                                                                 %Length for transition path

%% Steady state
kstar  = (par.alfa/(1/par.beta -1+par.delta))^(1/(1-par.alfa));               
ystar  = kstar^par.alfa;
cstar  = ystar - par.delta*kstar ;
kystar = (par.alfa/(1/par.beta -1+par.delta)) ;

%% Capital Grid
kmin=1; 
kmax=5; 
kstep=(kmax-kmin)/(par.n-1);                                               %Bounds and steps for Capital Grid
k=(kmin:kstep:kmax)';                                                      %Capital Grid

%% Consumption streams and utility flows
C = max(ones(par.n,1)*k'.^par.alfa - k*ones(1,par.n) + (1-par.delta).*ones(par.n,1)*k',0);   %Compute Consumption streams
if par.sigma==1
    U=log(C);                                                              %Utility levels 
else
    U=(C.^(1-par.sigma))/(1-par.sigma);
end

%% Value Function iteration
iter=0;                                                                    %Count iterations
diff=1;                                                                    %Norm 
tme = cputime;                                                             %Time
vinitial=zeros(1,par.n);                                                   %Initialize guessed value function

while diff>par.toler 
   [vrevised]=max(U+par.beta.*vinitial'*ones(1,par.n));                    %Maximize utility 
   diff=norm(vinitial-vrevised)/norm(vrevised);                            %Check convergence
   vinitial=vrevised;                                                      %Update guess
   iter = iter+1;
end
fprintf('\nFixed point solved via value function iteration took %d iterations and %d seconds\n', iter,cputime-tme);

%% Policy Functions
[vfinal,decision]=max(U+par.beta.*vrevised'*ones(1,par.n));                %Functional fixed point
kprime=k(decision);                                                        %Recover policy function for capital
c = k.^par.alfa - kprime + (1-par.delta).*k;                              %Recover policy function for consumption

%% Simulation of transition dynamics
kind=ones(1,par.T);                                                        %Initialize vector of indexes
ktran=ones(1,par.T+1);                                                     %Initialize vector of capital stocks 
ctran=ones(1,par.T);
ytran=ones(1,par.T);
kytran=ones(1,par.T);

kind(1)=1000;                                                              %Arbitrary starting point in terms of the index
ktran(1)= k(kind(1));                                                      %Arbitrary starting point in terms of capital stock

for t=2:par.T+1
    kind(t)=(decision(kind(t-1)));                                         %Evolution in index space
    ktran(t)=k(kind(t));                                                   %Evolution in capital space
    ytran(t-1)=ktran(t-1).^par.alfa;
    ctran(t-1)=ytran(t-1) - ktran(t) + (1-par.delta).*ktran(t-1);
    kytran(t-1)=ktran(t-1)/ytran(t-1);
end

%% Figures
figure(1); 
plot(k,vfinal,'b','linewidth',2);hold on; title('Value function'); grid on;

figure(2); 
subplot(121)
plot(k,kprime,'b','linewidth',2); hold on; plot(k,k,'r','linewidth',2); grid on
xlabel('capital stock t');ylabel('capital stock t+1'); title('policy function for capital')
subplot(122)
plot(k,c,'b','linewidth',2); hold on;  grid on
xlabel('capital stock t');ylabel('consumption t'); title('policy function for consumption')

figure(3); subplot(221);plot(ktran(1:par.T),'b','linewidth',2); 
hold on; plot(kstar*ones(par.T,1),'--r','linewidth',2); grid on;
xlabel('time period');ylabel('capital stock'); title('Transition path')
subplot(222);plot(ctran,'b','linewidth',2); grid on;
hold on; plot(cstar*ones(par.T,1),'--r','linewidth',2); grid on;
xlabel('time period');ylabel('consumption'); title('Transition path')
subplot(223);plot(ytran,'b','linewidth',2); grid on;
hold on; plot(ystar*ones(par.T,1),'--r','linewidth',2); grid on;
xlabel('time period');ylabel('output'); title('Transition path')
subplot(224);plot(kytran,'b','linewidth',2); grid on;
xlabel('time period');ylabel('capital-output ratio'); title('Transition path')
hold on; plot(kystar*ones(par.T,1),'--r','linewidth',2); grid on;








