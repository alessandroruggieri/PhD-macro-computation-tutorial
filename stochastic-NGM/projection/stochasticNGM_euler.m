%        Stochastic Neoclassical Growth Model - Projection Method         %
%          Program written by: Alessandro Ruggieri - UAB and BGSE         %
%                          Version 11/05/2016                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear
close all
fprintf('\nSolve Stochastic Neoclassical Growth Model using Local Projection Method\n');


%% Parametrization
alfa    = 0.33;     % Capital share 
beta    = 0.95;     % Discount factor
gamma   = 3.00;     % Risk aversion parameter
delta   = 0.07;     % Depreciation rate
sigma   = 0.01;     % Standard deviation for log noise
rho     = 0.95;     % Persistence of log technology shock
param   = [alfa,beta,gamma,delta,sigma,rho];

%% Steady State
kss= (alfa/(1/beta -1+delta))^(1/(1-alfa)); 
css= kss^alfa - delta*kss;	
yss= kss^alfa;	

%% Define grid for capital and shock
klow  = 0.5*kss;
khigh = 2*kss;
nk    = 10;
kstep = (khigh-klow)/(nk-1);
kgrid = (klow:kstep:khigh)';

zlow  =-3*sqrt(sigma^2/(1-rho^2));
zhigh = 3*sqrt(sigma^2/(1-rho^2));
nz    = 10;
zstep = (zhigh-zlow)/(nz-1);
zgrid = exp(zlow:zstep:zhigh)';

% Construct nodes and weights for numerical integration
nq = 5; [qnodes,qweights]= hernodes(nq);


%% Minimization routine
coeffin = [log(1-alfa*beta); 0.001; 0.001];
options = optimset('MaxFunEvals',1e+06,'MaxIter',1e+06,'TolFun',1e-10,'TolX',1e-10,'Display','iter','Algorithm','trust-region-reflective');
tme = cputime;
coeffout = fminunc(@(x) griderror(x,param, qnodes, qweights, kgrid, zgrid), coeffin,options);
disp(['the fixed point via PM took ' num2str(cputime-tme) ' seconds to complete'])

%% Simulations 
T=500;

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
kstart = 0.1*klow;
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

% capital grid for consumption policy
kgridfine = (klow:(khigh-klow)/1000:khigh)';

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




