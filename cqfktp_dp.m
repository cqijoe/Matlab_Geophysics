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




end

