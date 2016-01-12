function [ xt ] = cqfktpinv( tp, dt, dp, pleft, dx, xleft, nx, ...
    fk_background,term )
% Inverse Tau-P transform using fk domain method
% 
% input
% ------
% tp = tau - p domain data matrix
% dt = sampling rate in time
% dp = sampling rate in p
% pleft = what is the left-most p in your tp
% dx = sampling rate in x
% xleft = what is the left-most x in your xt data
% nx = how many x traces you have? (refer to xt_fk's column number)
% fk_background = (optional for better reconstruction)background fk domain 
%         data generated when doing foward tau-p transform using 
%         cqfktp function
%         ----------------
%         column number of fk_background must = nx
% term = interpolation term used
%
% output
% ------
% xt = reconstructed x-t domain data matrix

% 1D fft in t direction
[nt0,np0] = size(tp);
nt = 2^nextpow2(nt0);
fp = fft(tp,nt,1);
fp = fp(1:nt/2+1);

fnq = 1/2/dt;
knq = 1/2/dx;
axis_f = linspace(0,fnq,nt/2+1);
axis_k = linspace(0,2*fnq,nx); % k is in [0,2*knq)
index_f = 1:(nt/2+1);

% composing fk data
xt_fk = zeros(nt/2+1,nx);
% loop for each k
for n = 1:nx
    p_search = -axis_k(n)./axis_f;
    n_in_range = find( p_search )
    
    
end

end

