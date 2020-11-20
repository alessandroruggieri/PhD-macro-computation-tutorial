function [ValueFunctions,PolicyFunctions] = solve_vfi(p,n_grid,n_size,z_grid,z_size,z_prob,r,theta,co,cf)


% Revenue functions                                      
revenue = p.*repmat(reshape(z_grid,z_size,1),1,n_size).* ...
             repmat(reshape(n_grid,1,n_size),z_size,1).^theta;
 
% Profit functions 
profit = revenue - kron(n_grid,ones(1,z_size)') - co*ones(z_size,n_size);        
  
profit_v = single(zeros(z_size*n_size,n_size));    
for i=1:z_size
     profit_v((i-1)*n_size+1:i*n_size,:) = ones(n_size,1)*profit(i,:);
end

%Hiring and firing costs 
C = cf*abs(tril(kron(n_grid,ones(n_size,1))-kron(n_grid',ones(1,n_size)))) ; 
C = tril(C) - diag(diag(C));


% Allocate value and policy function
v_ini  = zeros(z_size,n_size);
v_upd  = zeros(z_size,n_size);
policy = zeros(z_size,n_size); 

d = 1;
tol =1e-5;

while d> tol

    vstay = z_prob*v_ini;    
          
    for i=1:z_size
        vcont= ones(n_size,1)*vstay(i,:) ;
        
        [v_upd(i,:), policy(i,:)] = max((1/(1+r))*(...
            (profit_v((i-1)*n_size+1:i*n_size,:)- C + max(vcont,zeros(n_size,n_size)))'));                             
    end
   
    d =  norm(v_upd-v_ini)/norm(v_ini);
    v_ini   = v_upd;
           
end

% Value functions
ValueFunctions = v_ini;

% Policy functions
PolicyFunctions = struct();

% indicator function for firms that exits
PolicyFunctions.indic_exit = zeros(z_size,n_size); 
PolicyFunctions.indic_exit (vstay > zeros(z_size,n_size)) = 1; 
    
% policy function for employment
PolicyFunctions.pol = policy;

% indicator function for firms that expand
PolicyFunctions.indic_hire = zeros(z_size,n_size);
PolicyFunctions.indic_hire (PolicyFunctions.pol > (ones(z_size,1)*linspace(1,n_size,n_size))) = 1;
 
% indicator function for firms that retain their size
PolicyFunctions.indic_rest = zeros (z_size,n_size);
PolicyFunctions.indic_rest (PolicyFunctions.pol == (ones(z_size,1)*linspace(1,n_size,n_size)))  = 1;

% indicator function for firms that fire
PolicyFunctions.indic_fire = zeros (z_size,n_size);
PolicyFunctions.indic_fire (PolicyFunctions.pol < (ones(z_size,1)*linspace(1,n_size,n_size)))  = 1;
 
