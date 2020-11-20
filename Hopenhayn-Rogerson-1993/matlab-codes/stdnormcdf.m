function f=stdnormcdf(X)

f=zeros(size(X));

f= 0.5 * erfc(-(X)./sqrt(2));
    





