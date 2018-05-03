function [x] = adjBoxPlotOUT(x,k)

nOrig = length(x);

%% medcouple computing (robust measure of skewness)
x(isnan(x))=[];
n = size(x,1);
x = sort(x);
z = x - repmat(median(x),n,1);
ip = find(z >= 0);
im = find(z <= 0);
p = size(ip,1); q = size(im,1);
[mi, mj] = ind2sub([p,q],1:p*q);
zp = z(ip);
zm = z(im);
h = (zp(mi)+zm(mj))./(zp(mi)-zm(mj));
[ipz]= find(zp==0); % row numbers of z+ = 0, i.e., x_{i} = median(x)
[imz]= find(zm==0); % row numbers of z- = 0, i.e., x_{i} = median(x)
piz = ismember(mi,ipz); % positions in mi where z+=0
pjz = ismember(mj,imz); % positions in mi where z-=0
zmember = piz+pjz; % same size as mi, mj
pijz = find(zmember == 2); % positions where z+ = z- = 0, i.e., x_{i} =
% = x_{j} = median(x)
[indi,indj] = ind2sub([p,q],pijz); % pxq matrix position of the zero entries
indi = indi - min(indi) + 1; % row position of the zero entries as if they
% were in a separated matrix
indj = indj - min(indj) + 1; % column position of the zero entries as if they were in a separated matrix
for i=1:size(pijz,2),
    if (indi(i) + indj(i) - 1) > size(find(z==0),1),
    h(pijz(i)) = 1;
    elseif (indi(i) + indj(i) - 1) < size(find(z==0),1),
    h(pijz(i)) = -1;
    else
    h(pijz(i)) = 0;
    end
end
MC = median(h);


%% Computing adjusted interquartile measures
if MC >=p 0
    int = [quantile(x, 0.25) - (k*exp(-3.5*MC))*(quantile(x, 0.75)-quantile(x, 0.25))  quantile(x, 0.75) + (k*exp(4*MC))*(quantile(x, 0.75)-quantile(x, 0.25))];
else
    int = [quantile(x, 0.25) - (k*exp(-4*MC))*(quantile(x, 0.75)-quantile(x, 0.25))  quantile(x, 0.75) + (k*exp(3.5*MC))*(quantile(x, 0.75)-quantile(x, 0.25))];
end

%% Detecting and removing outiers (values of x which fallen outside the lower and upper int)
pos_low = (find(x <= int(1)));
pos_upper = (find(x >= int(2)));
values_low = x(find(x <= int(1)));
values_upper = x(find(x >= int(2)));
x(find(x <= int(1))) = NaN;
x(find(x >= int(2))) = NaN;

if isempty(x(find(x <= int(1)))) && isempty(x(find(x >= int(2))))
    x = [x;zeros(nOrig-length(x),1)*NaN];
end




%results.outliers.pos.low = pos_low;
%results.outliers.pos.upper = pos_upper;
%results.outliers.values.low = values_low;
%results.outliers.values.upper = values_upper;
%results.adjMed = MC;
%results.x = x;

end

