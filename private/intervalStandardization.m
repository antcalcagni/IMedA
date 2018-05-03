function [ETA,C,R,PHI] = intervalStandardization(A)

n=size(A,1);
mean_A = (ones(2,1)'*A'*ones(n,1))*(1/(2*n));
var_A = ((((A.*A)*ones(2,1) + (A*[1;0]).*(A*[0;1]))'*ones(n,1))*(1/(3*n))) - (1/(4*n^2))*((ones(2,1)'*A'*ones(n,1)))^2;

ETA = (A-ones(n,2)*mean_A).*(ones(n,2)*sqrt(var_A)^-1);
C = median(ETA,2);
R = range(ETA,2);

PHI = ones(n,2)*mean_A + (ETA.*(ones(n,2)*sqrt(var_A))); 

end