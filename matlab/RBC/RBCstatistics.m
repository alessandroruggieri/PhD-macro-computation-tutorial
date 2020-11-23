%% Code to perform filtering analysis of US aggregate series
%  Author: Alessandro Ruggieri  (UAB and Barcelona GSE)
%  This version: December 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clear
close all;

%% Parameters
par.delta =0.05; %Depreciation
par.alfa  =0.33; %Employment elasticity
par.lambda=1600; %Smoothing coefficient, quarterly data

%% Load data
fprintf('\nLoading data ...\n')

% working age population
pop=xlsread('pop');
fprintf('\nPopulation data loaded\n')
logpop=log(pop);

% Employment
e=xlsread('emp');
fprintf('\nEmployment data loaded\n')
logn=log(e);

% average weekly hours
h=xlsread('hours');
fprintf('\nWeekly hours data loaded\n')
logh=log(h);

% total hours worked
he=e.*(h*4.3);
loghe=log(he);

% national account data
[data]= xlsread('data') ;
fprintf('\nNational account (GDP, consumption, investment, net export) loaded\n')

y=data(:,1);
c=data(:,2);
i=data(:,3);
nx = y-c-i;

logy=log(y);
logc=log(c);
logi=log(i);
lognx=log(nx);
N=size(logy,1);


%% Compute capital stoc
fprintf('\nCompute capital stock using perpetual invetory method\n')
ygrowth=(y(2:end)-y(1:end-1))./y(1:end-1); 
k=zeros(N,1);
k(1) = i(1)./(par.delta+mean(ygrowth));
for j=2:N
   k(j)=k(j-1)*(1-par.delta) + i(j); 
end
logk=log(k);

%% Compute TFP
fprintf('\nConstruct TFP\n')
tfp = y./(k.^par.alfa.*he.^(1-par.alfa));
logtfp = log(tfp); 
tfpgrowth=(tfp(2:end)-tfp(1:end-1))./tfp(1:end-1);

%% Filter
fprintf('\nFiltering...\n')

ytrend  =hpfilter(logy, par.lambda);
ctrend  =hpfilter(logc, par.lambda);
itrend  =hpfilter(logi, par.lambda);
nxtrend =hpfilter(lognx, par.lambda);
hetrend =hpfilter(loghe, par.lambda);
tfptrend=hpfilter(logtfp,par.lambda); 

ycycle=logy-ytrend;
ccycle=logc-ctrend;
icycle=logi-itrend;
nxcycle=lognx-nxtrend;
hecycle=loghe-hetrend;
tfpcycle=logtfp-tfptrend;

figure(1)
subplot(321)
plot(exp(logy)/exp(logy(1)));hold on; plot(exp(ytrend)/exp(ytrend(1)));grid on;title('GDP')
axis([0 200 1 5])
subplot(322)
plot(exp(logc)/exp(logc(1)));hold on; plot(exp(ctrend)/exp(ctrend(1)));grid on; title('Consumption')
axis([0 200 1 5])
subplot(323)
plot(exp(logi)/exp(logi(1)));hold on; plot(exp(itrend)/exp(itrend(1)));hold on;  grid on; title('Investment')
axis([0 200 1 8])
subplot(324)
plot(exp(loghe)/exp(loghe(1)));hold on; plot(exp(hetrend)/exp(hetrend(1)));grid on; title('Hours')
axis([0 200 1 3])
subplot(325)
plot(exp(lognx)/exp(lognx(1)));hold on; plot(exp(nxtrend)/exp(nxtrend(1)));hold on;  grid on; title('Net-Export')
axis([0 200 1 5])
subplot(326)
plot(exp(logtfp)/exp(logtfp(1)));hold on; plot(exp(tfptrend)/exp(tfptrend(1)));grid on; title('TFP')
axis([0 200 1 1.5])

%% Compute standard deviation of cyclical components
ystd  = 100*std(ycycle);
cstd  = 100*std(ccycle);
istd  = 100*std(icycle);
nxstd = 100*std(nxcycle);
hestd  = 100*std(hecycle);
tfpstd= 100*std(tfpcycle);

fprintf('\nStandard Deviation of Cyclical Components (percent)\n');
disp('');
disp(['GDP: ', num2str(ystd)])
disp(['Cons: ', num2str(cstd)])
disp(['Inv: ', num2str(istd)])
disp(['NX: ', num2str(nxstd)])
disp(['Hours: ', num2str(hestd)])
disp(['TFP: ', num2str(tfpstd)])
disp('');

%% Compute correlation with output
X=[ycycle,ccycle,icycle,nxcycle,hecycle,tfpcycle];
rho=corr(X);
fprintf('\nCorrelation of Cyclical Components with GDP\n');
disp('');
disp(['GDP: ', num2str(rho(1,1))])
disp(['Cons: ', num2str(rho(2,1))])
disp(['Inv: ', num2str(rho(3,1))])
disp(['NX: ', num2str(rho(4,1))])
disp(['Hours: ', num2str(rho(5,1))])
disp(['TFP: ',   num2str(rho(6,1))])

%% Compute autocorrelation 
fprintf('\nFirst-order Autocorrelation of Cyclical Components\n');
disp('');
X=[ycycle(2:end),ycycle(1:end-1),];
rho=corr(X);
disp(['GDP: ', num2str(rho(1,2))])
X=[ccycle(2:end),ccycle(1:end-1),];
rho=corr(X);
disp(['Cons: ', num2str(rho(1,2))])
X=[icycle(2:end),icycle(1:end-1),];
rho=corr(X);
disp(['Inv: ', num2str(rho(1,2))])
X=[nxcycle(2:end),nxcycle(1:end-1),];
rho=corr(X);
disp(['NX: ', num2str(rho(1,2))])
X=[hecycle(2:end),hecycle(1:end-1),];
rho=corr(X);
disp(['Hours: ', num2str(rho(1,2))])
X=[tfpcycle(2:end),tfpcycle(1:end-1),];
rho=corr(X);
disp(['TFP: ', num2str(rho(1,2))])



