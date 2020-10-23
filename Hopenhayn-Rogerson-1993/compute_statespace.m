% State Space

% Grid for capital
n_grid = linspace(n_min,n_max,n_size);
n_mat  = ones(n_size,1)*n_grid; % Rows are next period, columns are current period



% Grid for productivity
[z_grid , z_prob]=tauchen(z_size,muZ,roZ,sigmaZ);
z_grid = exp(z_grid);        
     
% Calculate the ergodic distribution
Trans= z_prob';
probst = (1/z_size)*ones(z_size,1); % initial distribution of states
test = 1;

while test > 10^(-8)
              probst1 = Trans*probst;
              test=max(abs(probst1-probst));
              probst = probst1;   
end
z_ergprob=probst';
          