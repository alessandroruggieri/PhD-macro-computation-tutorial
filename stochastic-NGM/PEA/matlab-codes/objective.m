% The objective function for the program
% PEAmbound.m 
%
% ksiz   initial coefficients 
% x      regressors
% y      explanatory variable 
% ---------------------------

function y = objective(ksiz,x)
y=exp(x*ksiz);