% function [x,w]	= hernodes(n);
% 
% Purpose:
% Compute the weights and nodes for Hermite-Polynomials
% 
% Inputs:
% n	The number of nodes to compute
% 
% Outputs:
% x	The nodes
% w	the associated weights
%
% Source:
% Taken from Numerical Recipes in Fortran77, 2ed
%
% -----------------------------------------------------------------

function [x,w]	= hernodes(n);

TOL	= sqrt(eps);
MAXIT 	= 30;
PIM4	= pi^(-1/4);

m	= (n+1)/2;
x	= zeros(n,1);
w	= zeros(n,1);

for i	= 1:m;

	% Initialize the first four roots
	% -------------------------------
	if i 	== 1; 
		z	= sqrt(2*n+1) - 1.85575 * (2*n+1)^(-1/6);
	elseif i== 2;
		z	= z - 1.14 * (n^.426)/z;
	elseif i== 3;
		z	= 1.86 * z - .86 * x(1);
	elseif i== 4;
		z	= 1.91 * z - .91 * x(2);
	else;
		z	= 2 * z - x(i-2);
	end;
	
	for its	= 1:MAXIT;
		p1 	= PIM4;
		p2 	= 0;
		for j	= 1:n;
			p3	= p2;
			p2	= p1;
			p1	= z*sqrt(2/j)*p2 - sqrt((j-1)/j)*p3;
		end;
		
		pp	= p2 * sqrt(2*n);
		z1	= z;
		z	= z1 - p1/pp;
		if abs(z-z1) < TOL; break; end;
	end;
	
	x(i)	= z;
	x(n+1-i)= -z;
	w(i)	= 2/(pp*pp);
	w(n+1-i)= w(i);
end;

% **********************************************************************

% **********************************************************************
