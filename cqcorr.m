function [ dcor ] = cqcorr( d, pilot )
% Correlate uncorrelated data with pilot trace and truncate
% 
% input
% -----
% d = uncorrelated data
% pilot = pilot trace
%
% output
% ------
% dcor = correlated data

pilot = pilot(:);
[s,lag] = xcorr(d(:,1),pilot);
n0 = find(lag==0);
nend = n0 + size(d,1) - 1;
dcor = zeros(size(d));
dcor(:,1) = s(n0:nend);

for k = 2:size(d,2)
    s = xcorr(d(:,k),pilot);
    dcor(:,k) = s(n0:nend);
end


end

