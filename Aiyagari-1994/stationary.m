function [y2]=stationary(Q)
%Q is the transition matrix of the Markov Chain. We use the method of state
%space reduction to compute the stationary distributions of Q.

P=Q;
[ns,~]=size(P);
n=ns;
n8=n;
while n8>1 
   for n1=1:(n8-1)
          for n2=1:(n8-1)
P(n1,n2)=P(n1,n2)+P(n1,n8)*P(n8,n2)/(1-P(n8,n8));
          end
   end
n8=n8-1;
end

y=(ones(n,1));
   n3=1;
while n3<n
    x=(ones(n3,1));
    for n5=1:n3
    x(n5)=y(n5);
    end        
    z=sum(x.*P(1:n3,n3+1))/(1-P(n3+1,n3+1));
        y(n3+1)=z/(1+z);
    for n6=1:n3
        y(n6)=y(n6)*(1-y(n3+1));
    end
    n3=n3+1;
end
y2=y;
end