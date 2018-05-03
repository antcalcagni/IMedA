function [cv] = covar(x,y)

a = cov(x,y);
cv=a(2,1);

end