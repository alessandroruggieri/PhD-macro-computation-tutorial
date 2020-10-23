%      Hugget (1993) Incomplete Market Model - Value Function Iteration    %
%          Program written by: Alessandro Ruggieri - UAB and BGSE         %
%                          Version 11/02/2016                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clear
close all
fprintf('\nSolve Incomplete Market Model a la Hugget 1993\n');


%% Calibration
beta   = 0.96;              % discount factor 
gamma  = 3;                 % risk aversion 
prob   = [ .95 .05; .2 .8]; % prob(i,j) = probability (A(t+1)=Aj | A(t) = Ai)
s      = size(prob,1);      % number of states
ehigh  = 1.0;               % high value for income
elow   = 0.8;               % low value for income
B      = 0;                 % Borrowing constraints


%% Asset Grid
n=100;                                           %Grid points
amin=B; 
amax=0.2;  
astep=(amax-amin)/(n-1);    %Bounds and steps for Asset Grid
a=(amin:astep:amax)';                             %Asset Grid

%% Initialize equilibrium price iteration
tic;
meanA=1;
toler=1e-10;
maxiter=15;
iter1=0;

qmin=1/2*beta;
qmax=4*beta;

while abs(meanA)>toler && iter1<maxiter
%% Set equilibrium price
q=qmin/2+qmax/2;

%%  Value function iteration
vinitial= zeros(n,s);
diff=1;                                              %Norm 
toler=1e-10;                                         %Tolerance level
iter2=0;                                              %Count iterations

Chigh= max(ehigh - q*a*ones(1,n) + ones(n,1)*a',0) ;
Clow = max(elow  - q*a*ones(1,n) + ones(n,1)*a',0) ;

if gamma==1
    Ulow=log(Clow); 
    Uhigh=log(Chigh);                                     %Utility levels 
else
    Ulow=(Clow.^(1-gamma))/(1-gamma); 
    Uhigh=(Chigh.^(1-gamma))/(1-gamma);
end

tme = cputime;
while diff > toler
  [vrevisedhigh]=max(Uhigh + beta*repmat(vinitial*prob(1,:)',1,n));
  [vrevisedlow] =max(Ulow  + beta*repmat(vinitial*prob(2,:)',1,n));
  vrevised=[vrevisedhigh' vrevisedlow'];
  diff=norm(norm(abs(vinitial-vrevised)));
  vinitial=vrevised;
  iter2 = iter2+1;
end
disp(['the fixed point via VFI took ' num2str(cputime-tme) ' seconds and ' num2str(iter2) ' iterations to complete'])


%% Policy Function
[vfinalhigh,decisionhigh]=max(Uhigh + beta*repmat(vinitial*prob(1,:)',1,n));
[vfinallow,decisionlow]=max(Ulow + beta*repmat(vinitial*prob(2,:)',1,n));
vfinal=[vfinalhigh',vfinallow'];
decision=[decisionhigh',decisionlow'];

%% Stationary distribution
% Form transition matrix from state at t (row) to the state at t+1 (column) 
% The eigenvector associated with the unit eigenvalue of transition matrix (transposed) is  the stationary distribution. 

g2=sparse(n,n);
g1=sparse(n,n);
for i=1:n
    g1(i,decisionhigh(i))=1;
    g2(i,decisionlow(i))=1;
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
aa=a(decision(:));
meanA=probst'*aa;

% Calculate measure over (k,A) pairs lambda has same dimensions as decis
lambda=zeros(n,2);
lambda(:)=probst;

% Calculate stationary distribution of capital 
proba=sum(lambda,2);     

% Update prices
if meanA>0
    qmin=q;
else
    qmax=q;
end
iter1=iter1+1;

end
tElapsed = toc;
disp(['the fixed point via bisection took ' num2str(tElapsed) ' seconds and ' num2str(iter1) ' iterations to complete'])


%% Policy Functions
deruleahigh=a(decision(:,1));                                    %Recover policy function for capital
derulealow =a(decision(:,2)); 
derulechigh=ehigh - q*deruleahigh +a;                            %Recover policy function for consumption
deruleclow=elow - q*derulealow + a ;                             %Recover policy function for consumption


disp('Parameters');
disp('');
disp('    gamma      beta          B   '); 
disp([ gamma        beta        B ]);
disp(''); 
disp('Results ');
disp('');
disp('    q              A     ');
disp([ q  meanA ]);


%% Plot Figures
figure(1); 
plot(a,vfinal,'b','linewidth',2);hold on; title('Value functions'); grid on;

figure(2); 
subplot(121)
plot(a,deruleahigh,'b-.','linewidth',2); hold on; plot(a,derulealow,'b--','linewidth',2); hold on; plot(a,a,'r','linewidth',2); grid on
xlabel('asset t');ylabel('asset t+1'); title('policy function for asset')
subplot(122)
plot(a,derulechigh,'b-.','linewidth',2); hold on; plot(a,deruleclow,'b--','linewidth',2); grid on
xlabel('asset t');ylabel('consumption t'); title('policy function for consumption')

figure(3)
plot(a,proba); title('Stationary Distribution of asset');
xlabel('asset'); ylabel('fraction of agents');

