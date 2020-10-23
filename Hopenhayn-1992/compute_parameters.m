%% Parameters
beta=0.8;                 %discounting rate
theta=0.64;               %production function parameter
cf=15;                    %fixed cost (approximated by mean employment)
ce=100;                   %entry cost
D = 300;                  %size of the market (exogeneous)

a=0.37;                   %autoregressive intercept (approximated by five year exit rate)
ro=0.93;                  %autoregressive parameter

vare=0.53*(1-theta)^2;    %vare/(1-theta)^2=0.53
varz=vare/(1-ro^2);       %simplify z=log(s)
sigma=sqrt(varz);         %standard deviation of z
stde=sqrt(vare);

Z=21;                     % Grid point for productivity
N=251;                    % Grid point for employment

toler=1e-06;