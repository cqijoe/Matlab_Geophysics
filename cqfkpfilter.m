function [ xt_filter, fk_filter, fk, f, k ] = ...
    cqfkpfilter( xt, dt, dx,us, f0, f1, tswp, pmin, pmax  )
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

% fk tranform
% -----------

% ------------ TEST: Padding zero traces around to reduce edge effect -----------
% N0 = 100;
% xt = [zeros(size(xt,1),N0),xt,zeros(size(xt,1),N0)];
% Conclusion: No effect on the output result
% ===============================================================================

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

doppler_factor = -us * p_table;

% Dragoset method WRONG
% phase_correction = -2 * pi * doppler_factor * tswp .* repmat(f, 1, nx).^2 ./ (f1 - f0);

% Chen method
f_mat = repmat(f, 1, nx);
phase_correction = -2 * pi * doppler_factor * tswp .* f_mat.* (f_mat-f0) ./ (f1 - f0);

fk_filter = fk .* exp(1i*phase_correction); % we have problem here!!

% inverse fk transform to xt domain
% ---------------------------------
fk_filter_negative = rot90(conj(fk_filter(2:end-1,:)),2);
fk_filter_negative = [fk_filter_negative(:,end),fk_filter_negative(:,1:end-1)];
fk_filter_two_side = [fk_filter; fk_filter_negative];


xt_filter = real(ifft2(fk_filter_two_side));
xt_filter = xt_filter(1:nt0,1:nx0);

% ---- test ---
% =============





end

