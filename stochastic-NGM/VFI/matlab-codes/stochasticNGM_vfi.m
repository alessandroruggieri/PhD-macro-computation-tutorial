%Stochastic Neoclassical Growth Model - Value function iteration          %
%          Program written by: Alessandro Ruggieri - UAB and BGSE         %
%                          Version 11/05/2016                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clear
close all
fprintf('\nSolve Stochastic Neoclassical Growth Model using VFI\n');

%% Parametrization
par.alfa=0.33;                                                             %Capital share of income
par.sigma=1;                                                               %Relative risk aversion
par.delta=0.1;                                                             %Depreciation rate
par.beta=0.95;                                                             %Discount factor
par.n=1000;                                                                %number of capital grid points
par.s=2;                                                                   %number of states
par.toler=1e-10;                                                           %Tolerance level

%% Productivity grid and transition matrix
Ahigh  = 1.5;                % high value for productivity
Alow   = 0.5;                % low value for productivity
probA   = [ .5 .5; .5 .5];   % probA(i,j) = probAability (A(t+1)=Aj | A(t) = Ai)

%% Capital grid 
kmin=0.1; 
kmax=10; 
kstep=(kmax-kmin)/(par.n-1);                                               %Bounds and steps for Capital Grid
k=(kmin:kstep:kmax)';                                                      %Capital Grid

%% Consumption streams and utility flows
Chigh=max(Ahigh.*ones(par.n,1)*k'.^par.alfa - k*ones(1,par.n) + (1-par.delta).*ones(par.n,1)*k',0);   %Consumption streams
Clow =max( Alow.*ones(par.n,1)*k'.^par.alfa - k*ones(1,par.n) + (1-par.delta).*ones(par.n,1)*k',0);   

if par.sigma==1
    Ulow=log(Clow); 
    Uhigh=log(Chigh);                                                             %Utility levels 
else
    Ulow=(Clow.^(1-par.sigma))/(1-par.sigma); 
    Uhigh=(Chigh.^(1-par.sigma))/(1-par.sigma);
end


%%  Value function iteration
vinitial= zeros(par.n,par.s);
vrevised= zeros(par.n,par.s);
diff=1;                                                                           %Norm 
iter=0;                                                                           %Count iterations
tme = cputime;

while diff > par.toler
  [vrevisedhigh]=max(Uhigh + par.beta*repmat(vinitial*probA(1,:)',1,par.n));
  [vrevisedlow] =max(Ulow  + par.beta*repmat(vinitial*probA(2,:)',1,par.n));
  vrevised=[vrevisedhigh' vrevisedlow'];
  diff=norm(norm(abs(vinitial-vrevised)));
  vinitial=vrevised;
  iter = iter+1;
end
fprintf('\nFixed point solved via VFI took %d iterations and %d seconds\n', iter,cputime-tme);


%% Policy Function
[vfinalhigh,decisionhigh]=max(Uhigh + par.beta*repmat(vinitial*probA(1,:)',1,par.n));
[vfinallow,decisionlow]=max(Ulow + par.beta*repmat(vinitial*probA(2,:)',1,par.n));
vfinal=[vfinalhigh',vfinallow'];
decision=[decisionhigh',decisionlow'];

kprime(:,1) =k(decision(:,1));                                             %Recover policy function for capital
kprime(:,2) =k(decision(:,2));

c(:,1)=Ahigh.*k.^par.alfa - kprime(:,1) + (1-par.delta).*k;                %Recover policy function for consumption
c(:,2)=Alow.*k.^par.alfa  - kprime(:,2) + (1-par.delta).*k;                


%% Stationary distribution
% Form transition matrix from state at t (row) to the state at t+1 (column) 
% The eigenvector associated with the unit eigenvalue of transition matrix (transposed) is  the stationary distribution. 

g2=sparse(par.n,par.n);
g1=sparse(par.n,par.n);
for i=1:par.n
    g1(i,decisionhigh(i))=1;
    g2(i,decisionlow(i))=1;
end
trans=[ probA(1,1)*g1 probA(1,2)*g1; probA(2,1)*g2 probA(2,2)*g2];
trans= trans';
probAst = (1/(2*par.n))*ones(2*par.n,1);
test = 1;
while test > par.toler
   probAst1 = trans*probAst;
   test=max(abs(probAst1-probAst));
   probAst = probAst1;
end

% Vectorize the decision rule to calculate mean level of capital
kk=k(decision(:));
meanK=probAst'*kk;

% Calculate measure over (k,A) pairs lambda has same dimensions as decis
lambda=zeros(par.n,2);
lambda(:)=probAst;

% Calculate stationary distribution of capital 
probAk=sum(lambda');     
probAk=probAk';


%% Figures
figure(1); 
plot(k,vfinal,'b','linewidth',2);hold on; title('Value functions'); grid on;

figure(2); 
subplot(121)
plot(k,kprime(:,1),'b-.','linewidth',2); hold on; plot(k,kprime(:,2),'b--','linewidth',2); hold on; plot(k,k,'r','linewidth',2); grid on
xlabel('capital stock t');ylabel('capital stock t+1'); title('policy function for capital')
subplot(122)
plot(k,c(:,1),'b-.','linewidth',2); hold on; plot(k,c(:,2),'b--','linewidth',2); grid on
xlabel('capital stock t');ylabel('consumption t'); title('policy function for consumption')


figure(3)
plot(k,probAk); title('Stationary Distribution of capital');
xlabel('capital');
ylabel('fraction of agents');


