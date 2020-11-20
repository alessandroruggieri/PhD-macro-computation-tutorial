% Parameters
n_min      = 1;       % Min number of workers
n_max      = 1000;    % Max number of workers
n_size     = 300;     % Number of workers grid points
z_size     = 10;       % Number of productivity grid points

theta      = 0.85;     % Capital share
delta      = 0.15;    % Separation rate

muZ        = 0;       % Productivity scalar
roZ        = 0.95;    % Productivity persistence
sigmaZ     = 0.10;    % Productivity volatiltiy

r          = 0.04;    % Interest rate
beta       = 1/(1+r); % Discount factor
cf         = 3;       % Adjustment cost (firing cost)
co         = 4;       % Fixed cost
ce         = 25;      % Entry cost
D          = 400;     % Demand shifter

tol        = 1e-5;    % Tolerance
d          = 1;       % Convergence criteria
