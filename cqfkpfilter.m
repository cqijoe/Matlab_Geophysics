function [ xt_filter, fk_filter, fk, f, k ] = ...
    cqfkpfilter( xt, dt, dx,us, f0, f1, tswp, pshift_x0,pmin, pmax  )
% Remove doppler effects in fk domain by searching for k - p pair
%
% input
% -----
% xt = xt domain data, either cdp, crg or csg
% dt = time sampling rate
% dx = space sampling rate
% us = boat speed
% f0 = start sweeping frequency
% f1 = end sweeping frequency (assuming linear sweeping)
% tswp = sweeping duration in seconds
% pshift_x0 = (optional) [pshift, x0],such that preal = p - pshift. pshift > 0 down shift
%             LMO pshift < 0 up shift LMO. x0 term is for calculatinf how much time shift is needed.
%             x0 indicates which trace is with offset = 0.
% pmin = (optional) minimum p to filter
% pmax = (optional) maximum p to filter
%
% output
% ------
% xt_filter = xt data after filtering
% fk_filter = fk data after filtering (positive f only) [o,2*pi)
% fk = fk data before filtering (positive f only) [0,2*p)
% f = frequency sampling points
% k = wavenumber sampling points

% Before everything else, check if to linear shift
if exist('pshift_x0','var') && ~isempty(pshift_x0)
    pshift = pshift_x0(1);
    x0 = pshift_x0(2);
    xvector = ([1:size(xt,2)] - x0) * dx;
    shift = -pshift * xvector;
    % decide how much zeros need to be padded
    npad_below = floor(max(shift)/dt);
    if npad_below > 0
        % pad below
        xt = [xt;zeros(npad_below,size(xt,2))];
    else 
        npad_below = 0;  %  will be used at the end of program
    end
    npad_up = floor(min(shift)/dt);
    if npad_up < 0
        npad_up = abs(npad_up); %  will be used at the end of program
        % pad above
        xt = [zeros(npad_up,size(xt,2));xt];
    else
        npad_up = 0; %  will be used at the end of program
    end
    [ xt ] = cqlinshift( xt, dt, shift, 'spline' );
    linear_shift = true;
else
    linear_shift = false;
end

% fk tranform
% -----------


[nt0, nx0] = size(xt);
nt = 2^nextpow2(nt0);
nx = 2^nextpow2(nx0);
fk_two_side = fft2(xt, nt, nx);
fk = fk_two_side(1:nt/2+1,:);
f = linspace(0, 1/2/dt, nt/2+1); % original frequency sampling points
f = f(:);

% generate k such that k is the wavenumber sampling points from [0,knyq] + (-knqy,0)
k = linspace(0, 1/2/dx, nx/2+1); k = [k,-fliplr(k(2:end-1))];


% do phase correction
% --------------------------------------------
p_table = -repmat(k, nt/2+1, 1)./repmat(f, 1, nx); 
p_table(f==0,:) = 0;
if f0 < f1
    p_table(f<f0 | f>f1,:) = 0; 
else
    p_table(f<f1 | f>f0,:) = 0;
end
if exist('pmin','var') && ~isempty(pmin)
    p_table(p_table<pmin) = 0;
end
if exist('pmax','var') && ~isempty(pmax)
    p_table(p_table>pmax) = 0;
end
if linear_shift
    p_table = p_table + pshift;
end
doppler_factor = -us * p_table;

% Dragoset method WRONG
% phase_correction = -2 * pi * doppler_factor * tswp .* repmat(f, 1, nx).^2 ./ (f1 - f0);

% Chen method
f_mat = repmat(f, 1, nx);
phase_correction = -2 * pi * doppler_factor * tswp .* f_mat.* (f_mat-f0) ./ (f1 - f0);

fk_filter = fk .* exp(1i*phase_correction); 

% inverse fk transform to xt domain
% ---------------------------------
fk_filter_negative = rot90(conj(fk_filter(2:end-1,:)),2);
fk_filter_negative = [fk_filter_negative(:,end),fk_filter_negative(:,1:end-1)];
fk_filter_two_side = [fk_filter; fk_filter_negative];


xt_filter = real(ifft2(fk_filter_two_side));
xt_filter = xt_filter(1:nt0,1:nx0);

% final truncating
if linear_shift
    xt_filter = cqlinshift( xt_filter, dt, -shift, 'spline' );
    xt_filter = xt_filter(npad_up+1:end-npad_below,:);
end





end

