function [ xt, xt_fk ] = cqfktpinv( tp, dt, dp, pmin, dx, xleft, nx, ...
    fk_background,term )
% Inverse Tau-P transform using fk domain method
% 
% input
% ------
% tp = tau - p domain data matrix
% dt = sampling rate in time
% dp = p sampling rate
% pmin = minimum p
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
fp = fp(1:nt/2+1,:);
% recompose fp data in p direction
fp = [fp(:,(np0+1)/2:np0),fp(:,2:(np0-1)/2)];
% now pmin is pmin + dp
min_p = pmin + dp;
max_p = pmin + dp * (np0-1); % pmax doesn't change
% dpinv = 1/2/( (np0+1)/2 )/dp;
dpinv = 1/(np0-1)/dp;

fnq = 1/2/dt;
knq = 1/2/dx;
axis_f = linspace(0,fnq,nt/2+1);
axis_f = axis_f(2:end);
axis_k = linspace(0,knq,nx/2+1);
axis_k = [axis_k,-fliplr(axis_k(2:end-1))];


% composing fk data
xt_fk = zeros(nt/2+1,nx);
% loop for each k
for n = 1:nx
    p_search = -axis_k(n)./axis_f;
    n_in_range = find( p_search>=min_p & p_search<=max_p );
    p_search = p_search(n_in_range);
    for m = 1:length(p_search)
        xt_fk(n_in_range(m),n) = ...
            cqcomplex_interp(fp(n_in_range(m),:),dpinv, p_search(m),term);
    end
end

% if user has background, added back ground at locations calculated
for n = 1:nx
    if axis_k(n)>=0
        fbound = -axis_k(n)/min_p;
    else
        fbound = -axis_k(n)/max_p;
    end
    n_add = find(axis_f<fbound);
    xt_fk(n_add,n) = fk_background(n_add,n);
end
% then come back to xt domain
kmat = repmat(axis_k,nt/2+1,1);
xt_fk = xt_fk .* exp( 1i*2*pi*kmat*xleft );
xt_fk = [real(xt_fk);flipud(real(xt_fk(2:end-1,:)))] + ...
    1i * [imag(xt_fk);-flipud(imag(xt_fk(2:end-1,:)))];
xt = real(ifft2(xt_fk));

end

