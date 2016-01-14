function [ tp, fk_bg ] = cqfktp( xt, dt, dx, xleft, p, n_pad,term, frange )
% Implement tau-p transform in f-k domain. Reference: Wade's Thesis in
% UH, 1989
%
% input
% -----
% xt = 2D matrix with x-t domain data
% dt = sampling rate in t direction
% dx = sampling rate in x direction
% xleft = what is the x coordinate for the left-most input trace?
% p = vector for slope (dt/dx)
% n_pad = how many times of 0 traces will be padded in xt data
% term = number of term in complex number interpolation
% frange = [0,1/2/dt] (default) frequency range,
%          this gives you the ability of doing frequency filtering in
%          fk domain.
%
% output
% ------
% tp = tau - p data
% fk_bg = fk data background, [0,2*pi)

if ~exist('frange','var')||isempty(frange)
    frange = [0,1/2/dt];
end

min_f = min(frange);
max_f = max(frange);

% F-K Transform
[nt0,nx0] = size(xt);

nx = 2^nextpow2(nx0)*n_pad;
nt = 2^nextpow2(nt0);

xt_fk = fft2(xt,nt,nx);
xt_fk = xt_fk(1:nt/2+1,:);

fnyq = 1/2/dt;
knyq = 1/2/dx;

f = linspace(0,fnyq,nt/2+1);
f_index_range = find(f>=min_f & f<=max_f);
f = f(f_index_range);

k = linspace(0,knyq,nx/2+1);
k = [-fliplr(k(2:end-1)),k];

kmat = repmat(k,nt/2+1,1);
xt_fk = [xt_fk(:,nx/2+2:end),xt_fk(:,1:nx/2+1)];

xt_fk = xt_fk .* exp(1i*2*pi*kmat*(-xleft));
fk_bg = [xt_fk(:,nx/2:end),xt_fk(:,1:nx/2-1)];


fp = zeros(nt/2+1,length(p));

for m = 1:length(p)
    index_of_f = f_index_range;
    k_interp = -p(m) * f; % each f we find a k to interpolate
    % filter out of outsider
    n_rm = find(k_interp<min(k) | k_interp>max(k));
    k_interp(n_rm) = [];
    index_of_f(n_rm) = [];
    
    % do interpolation 
    for iter = 1:length(k_interp)
        this_ind = index_of_f(iter);
        fp(this_ind,m) = ...
            cqcomplex_interp( xt_fk(this_ind,:), k, k_interp(iter), term );
    end
end

fp = [real(fp);flipud(real(fp(2:end-1,:)))] + ...
    1i*[imag(fp);-flipud(imag(fp(2:end-1,:)))];
tp = real(ifft(fp,[],1));


end

