%      Hopenhayn and Rogerson (1993) Firm Dynamics Market Model - VFI     %
%          Program written by: Alessandro Ruggieri - UAB and BGSE         %
%                          Version 11/02/2016                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clear
close all
fprintf('\nSolve Firm Dynamics Model as in Hopenhayn and Rogerson (1993)\n');

clear;
clc;
tic;

% Parameters
compute_parameters;

% State Space
compute_statespace;

% Price boundaries
pmin=0.01;
pmax=10;

while d>tol
    
    % Price of goods    
    price = (pmin+pmax)/2;

    % Value function iteration
    [v,PolicyFunctions] = solve_vfi(price,n_grid,n_size,z_grid,z_size,z_prob, r,theta,co,cf);

    %Value of entry    
    value = z_ergprob*v(:,1);     

    % Update price till EV=ce
    if value<ce
     pmin=price;
    else
     pmax=price;
    end
    fprintf('\nEquilibrium price %d\n', price)
    
    % Check convergence
    d=abs(value-ce)/ce;
    
end

%Given the value function and policy function, iterate on the industry structure
%until it converges
d=1;
muinitial =ones(z_size,n_size)./(n_size*z_size);
murevised =zeros(z_size,n_size);

while d>tol
    for j=1:z_size
            murevised(j,1) =  (PolicyFunctions.indic_exit(j,1)==1) .* sum( ((muinitial'*z_prob(:,j))' + z_ergprob(j)).* (PolicyFunctions.pol(j,:)==1) );
        for i=2:n_size
            murevised(j,i) =  (PolicyFunctions.indic_exit(j,i)==1) .* sum( ((muinitial'*z_prob(:,j))' ).* sum(PolicyFunctions.pol(j,:)==i) ) ;
        end
    end
    murevised=murevised./sum(sum(murevised));
    d=norm(murevised-muinitial);
    muinitial=murevised;
end
    

%% Calculating the entry mass M 
% Using equilibrium condition in goods market
Pstar=price;
yd =D-P;  
ys =sum(sum((repmat(z_grid',1,n_size).*n_grid(PolicyFunctions.pol).^theta).*murevised ));
Mstar = yd/ys;
N = sum(sum(n_grid(PolicyFunctions.pol).*murevised ));
 
disp('Results ');
disp('');
disp('    Price     Firms     Avg.size       ');
disp([ Pstar  Mstar     N   ]);
 