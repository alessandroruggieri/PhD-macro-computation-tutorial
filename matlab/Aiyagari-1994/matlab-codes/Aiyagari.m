%      Aiyagari (1994) Incomplete Market Model - Value Function Iteration    %
%          Program written by: Alessandro Ruggieri - UAB and BGSE         %
%                          Version 11/02/2016                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clear
close all
fprintf('\nSolve Incomplete Market Model a la Aiyagari 1994\n');


%% Calibration
beta   = 0.96;              % discount factor 
gamma  = 1;                 % risk aversion 
alfa   = 1/3;               % capital share
delta  = 0.035;             % depreciation
prob   = [ .7 .3; .5 .5];   % prob(i,j) = probability (A(t+1)=Aj | A(t) = Ai)
s      = size(prob,1);      % number of states
b      = 0.1;               % unemployment benefit (share of equilibrium wage)
F      = -0.1;              % Borrowing constraints
A      = 0.1;                 % Technology parameter

%% Initial value for capial
n=100;                      % number of grid points
kmin=0; 
kmax=1;  
kstep=(kmax-kmin)/(n-1);    
k=(kmin:kstep:kmax)';


%% Stationary labor supply (onle extensive margin accounted for)
probss= stationary(prob);
N=probss(1);

%% Initialize equilibrium price iteration
vinitial= zeros(n,s);
toler=1e-07;
maxiter=70;
iter1=0;
iter2=0;
diff1=1;
diff2=1;


% homotopy and initial condition
r=0.01;
g=0.05;

tic;
while  (diff1> toler) && (iter1 <= maxiter) 

   %Compute Capital stock from firm FOCS 
   K = ((alfa*A*N^(1-alfa))/r)^(1/(1-alfa));
   
   % Compute wage bill
   w = (1-alfa)*A* K^(alfa)* N^(-alfa);
   
   % Compute consumption 
   Cemp= max(w*ones(n,n) + (r+1-delta)*ones(n,1)*k' - k*ones(1,n),0) ;
   Cunemp = max(b*w*ones(n,n) + (r+1-delta)*ones(n,1)*k' - k*ones(1,n),0) ;
   % Compute utility
    if gamma==1
        Uemp=log(Cemp); 
        Uunemp=log(Cunemp);                                     %Utility levels 
    else
        Uemp=(Cemp.^(1-gamma))/(1-gamma); 
        Uunemp=(Cunemp.^(1-gamma))/(1-gamma);
    end

tme = cputime;
while diff2 > toler
  [vrevisedemp]   =max(Uemp    + beta*repmat(vinitial*prob(1,:)',1,n));
  [vrevisedunemp] =max(Uunemp  + beta*repmat(vinitial*prob(2,:)',1,n));
  vrevised=[vrevisedemp' vrevisedunemp'];
  diff2=norm(norm(abs(vinitial-vrevised)));
  vinitial=vrevised;
  iter2 = iter2+1;
end
disp(['the fixed point via VFI took ' num2str(cputime-tme) ' seconds and ' num2str(iter2) ' iterations to complete'])

%% Policy Function
[vfinalemp,decisionemp]   =max(Uemp    + beta*repmat(vinitial*prob(1,:)',1,n));
[vfinalunemp,decisionunemp] =max(Uunemp  + beta*repmat(vinitial*prob(2,:)',1,n));
vfinal=[vfinalemp',vfinalunemp'];
decision=[decisionemp',decisionunemp'];


%% Stationary distribution
% Form transition matrix from state at t (row) to the state at t+1 (column) 
% The eigenvector associated with the unit eigenvalue of transition matrix (transposed) is  the stationary distribution. 

g2=sparse(n,n);
g1=sparse(n,n);
for i=1:n
    g1(i,decisionemp(i))=1;
    g2(i,decisionunemp(i))=1;
end
trans=[ prob(1,1)*g1 prob(1,2)*g1; prob(2,1)*g2 prob(2,2)*g2];
trans= trans';
probst = (1/(2*n))*ones(2*n,1);
test = 1;
while test > toler
   probst1 = trans*probst;
   test=max(abs(probst1-probst));
   probst = probst1;
end

% Vectorize the decision rule to calculate mean level of capital
kk=k(decision(:));
meanK=probst'*kk;

% Calculate measure over (k,A) pairs lambda has same dimensions as decis
lambda=zeros(n,2);
lambda(:)=probst;

% Calculate stationary distribution of capital 
probk=sum(lambda,2);     

%Compute Capital stock from firm FOCS 
rnew = alfa*A* meanK^(alfa-1)* N^(1-alfa);

% update interest rate
diff1 = abs((r-rnew));
r = g*rnew + (1-g)*r;
iter1 = iter1+1;
diff2 = 1;

end
tElapsed = toc;
disp(['the fixed point via bisection took ' num2str(tElapsed) ' seconds and ' num2str(iter1) ' iterations to complete'])

%% Policy Functions
derulekemp=k(decision(:,1));                                    %Recover policy function for capital
derulekunemp =k(decision(:,2)); 

figure(1); 
plot(k,vfinal,'b','linewidth',2);hold on; title('Value functions'); grid on;

figure(2); 
plot(k,derulekemp,'b-.','linewidth',2); hold on; plot(k,derulekunemp,'b--','linewidth',2); hold on; plot(k,k,'r','linewidth',2); grid on
xlabel('capital stock t');ylabel('capital stock t+1'); title('policy function for capital')


figure(3)
plot(k,probk); title('Stationary Distribution of capital');
xlabel('capital'); ylabel('fraction of agents');




