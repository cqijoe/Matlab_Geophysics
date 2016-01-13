function [ p, dp ] = cqfktp_dp( dx, fmax, n, prange )
% Return p that could be used by cqfktp and cqfktpinv functions
%
% input
% ------
% dx = x sampling rate
% fmax = maximum frequency in Hz
% n = number of traces in your original xt domain data
% prange = [pmin,pmax]
%
% output
% ------
% p = vector of p
% dp = p sampling rate

n = 2^nextpow2(n)*4;
dk = 1/n/dx;
dp = dk/fmax;
p = min(prange):dp:max(prange);
np = length(p);
if mod(np,2)==0
    np = np + 1;
end
p = linspace(min(prange),max(prange),np);
dp = p(2) - p(1);


end

