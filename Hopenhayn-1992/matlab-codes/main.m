%      Hopenhayn (1992) Firm Dynamics Market Model - VFI                  %
%          Program written by: Alessandro Ruggieri - UAB and BGSE         %
%                          Version 11/02/2016                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clear
close all
fprintf('\nSolve Firm Dynamics Model as in Hopenhayn (1992)\n');

 
%% Parameters
compute_parameters;

%% State Space
compute_statespace;

%% Iterate over price of goods
d=1;
% Boundaries for price
pmin=0.01; 
pmax=100;

while d>toler

    % Guess prices    
    price=(pmin+pmax)/2;

    % Solve firm value function iteration
    [vinitial,dr,exit] = solve_vfi(price,z,Z,n,N,theta,beta,cf,p);

    %Compute the decision rule for labor
    decrule=zeros(1,Z);
    for i=1:Z
        decrule(i)=n(dr(i));
    end

    %Define expected value of entrant
    value=inidis*vinitial';

    % Update price till EV=ce
    if value<ce
     pmin=price;
    else
     pmax=price;
    end

    % Check convergence
    d=abs(value-ce)/ce;
end


%Given the value function and policy function, iterate on the industry structure
%until it converges
d=1; 
muinitial=inidis;
while d>toler
    muexit=muinitial.*exit;                 %exit decision
    mustay=muexit*p;                        %update for the incumbents stay
    muentry=mustay+inidis;                   %entry 
    murevised=muentry./sum(muentry);
    d=norm(murevised-muinitial)/norm(murevised);
    muinitial=murevised;
end

%% Calculating the entry mass M 
% Using equilibrium condition in goods market
y=D-price;  
Xstar=z(Z-sum(exit));
Pstar=price;
Size = (decrule)*murevised'; 
Y=(decrule.^theta.*z)*murevised';
Mstar=y/[Y+(decrule.^theta.*z)*inidis'];
Exrate=sum(murevised(1:Z-sum(exit)))*100;

disp('Results ');
disp('');
disp('    Price     Firms    Avg.size       ');
disp([ Pstar  Mstar    Size   ]);
 