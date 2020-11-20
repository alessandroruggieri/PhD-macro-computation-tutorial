%Stochastic Neoclassical Growth Model - Parametrized Expectation Algorithm%
%          Program written by: Alessandro Ruggieri - UAB and BGSE         %
%                          Version 11/05/2016                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear
close all
fprintf('\nSolve Stochastic Neoclassical Growth Model using Parametrized Expectation Algorithm\n');


%% Parametrization
alfa    = 0.33;     % Capital share 
beta    = 0.95;     % Discount factor
gamma   = 1;        % Risk aversion parameter
delta   = 0.10;     % Depreciation rate
sigma   = 0.01;     % Standard deviation for log noise
rho     = 0.95;     % Persistence of log technology shock
T       = 1000;     % Length of simulation

%% Steady State
kss=(alfa/(1/beta -1+delta))^(1/(1-alfa));  
css= kss^alfa - delta*kss;	
yss= kss^alfa;
ess= css^(-gamma)*( 1-delta+alfa*kss^(alfa-1) );


%% Allocate memory for the simulates series
z= zeros(T,1);      % Technology Shocks
k= zeros(T+1,1);    % Capital Stock
c= zeros(T,1);      % Consumption
e= zeros(T-1,1);    % Conditional expectation


%% Initial values of capital and technology shock
k(1)= kss;          % Initial value of capital 
z(1)= 1;            % Initial value of shock

%% Simulate AR(1) technology process
epsi = randn(T,1)*sigma;
for t = 2:T
    z(t) = z(t-1)^rho*exp(epsi(t)); 
end

%% Set-up
coeff  = [log(ess); 0.00001; 0.00001];    % Initial polinomial coefficients  
toler  = 1e-10;            	              % Convergence criterion
update = 1;            	                  % Updating parameter for homotopy 
iter   = 0;                               % Set iteration counter
diff   = 1; 					          % Set convergence criterion
up_bound  = 2*kss;                        % Upper bound for capital
low_bound = 0.5*kss;                      % Lower bound for capital

tme = cputime;

%% Iteration
while (diff > toler)  || (hit==1)
hit= 0;                                   % Indicator, 1= bound is hit ;  
for t = 1:T
      % parametrize expectation
      uprime = exp( coeff(1) + coeff(2)*log(k(t)) + coeff(3)*log(z(t)));
      % compute consumption choice for euler equation
      c(t) = ( beta*uprime )^(-1/gamma);
      % compute capital choice from budget constraints
      k(t+1)=z(t)*k(t)^alfa - c(t) + (1-delta)*k(t);
      % taking care of boundaries
      if      k(t+1)>up_bound ;  k(t+1)=up_bound; hit=1; 
      elseif  k(t+1)<low_bound;  k(t+1)=low_bound; hit=1;
      end
  
end
% Construct expectation from simulated data
e(1:T-1) = c(2:T).^(-gamma).*(z(2:T).*alfa.*k(2:T).^(alfa-1)+1-delta);
    
% Compute 'coeffout' by using non-linear least square regression
x        = [ones(T-1,1) log( k(1:T-1) ) log( z(1:T-1) )];   % Regressors 
coeffnlls= nlinfit(x,e,'objective',coeff);                % NLLS regression
diff     = norm(coeff-coeffnlls)/norm(coeff);            % Display difference between 
coeff    = update*coeffnlls + (1-update)*coeff;            % Update the coefficients (homotopy)
iter     = iter+1;                                          % Update number of iteration     
end
disp(['the fixed point via PEA took ' num2str(cputime-tme) ' seconds and ' num2str(iter) ' iterations to complete'])


% Policy Functions
coeffout = x\log(c(1:end-1));

% Plot the consumption policy function as a function of k
% Define grid for capital 
klow  = 0.5*kss;
khigh = 2*kss;
kgridfine = (klow:(khigh-klow)/1000:khigh)';



%% Simulations 
T=1000;

% Initialize value productivity
zsim=zeros(T,1);
zstart=0;
zsim(1)=zstart;

% Draw random shock
for t=2:T
    zsim(t) = rho*zsim(t-1) + sigma*randn(1,1);
end
zsim=exp(zsim);

% Initialize value capital
kstart = 0.1*low_bound;
ksim=zeros(T+1,1);
ksim(1)=kstart;

% Allocate memory
csim=zeros(T,1);
ysim=zeros(T,1);

for t=2:T+1
    csim(t-1) = consfun(ksim(t-1),zsim(t-1),coeffout);
    ysim(t-1) = zsim(t-1).*ksim(t-1).^alfa;
    ksim(t)   = (1-delta).*ksim(t-1) + zsim(t-1).*ksim(t-1).^alfa - consfun(ksim(t-1),zsim(t-1),coeffout);
end


%% Figures
figure(1)
subplot(221)
plot(kgridfine,consfun(kgridfine,0.9,coeffout), ...
     kgridfine,consfun(kgridfine,1.0,coeffout), ...
     kgridfine,consfun(kgridfine,1.1,coeffout)); 
xlabel('capital t');  grid on;
title('Consumption t')

subplot(222)
plot(kgridfine, 0.9.*kgridfine.^alfa - consfun(kgridfine,0.9,coeffout), ...
     kgridfine, 1.0.*kgridfine.^alfa - consfun(kgridfine,1.0,coeffout), ...
     kgridfine, 1.1.*kgridfine.^alfa - consfun(kgridfine,1.1,coeffout));
xlabel('capital t');  grid on;
title('Investment t')

subplot(223)
plot(kgridfine, 0.9.*kgridfine.^alfa, ...
     kgridfine, 1.0.*kgridfine.^alfa,  ...
     kgridfine, 1.1.*kgridfine.^alfa);
xlabel('capital t');  grid on;
title('Output t')

subplot(224)
plot(kgridfine, 0.9.*kgridfine.^alfa + (1-delta).*kgridfine - consfun(kgridfine,0.9,coeffout), ...
     kgridfine, 1.0.*kgridfine.^alfa + (1-delta).*kgridfine - consfun(kgridfine,1.0,coeffout), ...
     kgridfine, 1.1.*kgridfine.^alfa + (1-delta).*kgridfine - consfun(kgridfine,1.1,coeffout)); 
xlabel('capital t'); grid on;
title('Capital t+1')


figure(2)
subplot(221)
plot(zsim); xlabel('time'); ylabel('zt'); title('Productivity')
grid on
subplot(222)
plot(csim);hold on;plot(css*ones(T,1),'--r');xlabel('time');ylabel('ct');
title('Consumption')
grid on
subplot(223);plot(ksim(1:T));hold on;plot(kss*ones(T,1),'--r');xlabel('time');ylabel('kt');title('Capital')
grid on
subplot(224);plot(ysim);hold on;plot(yss*ones(T,1),'--r');xlabel('time');ylabel('yt');title('Production')
grid on





