function [ tp, xt_fk, fp ] = cqfktp( xt, dt, dx, xleft, p, term )
% Implement tau-p transform in f-k domain. Reference: Wade's Thesis in
% UH, 1989
%
% input
% -----
% xt = 2D matrix with x-t domain data
% dt = sampling rate in t direction
% dx = sampling rate in x direction
% p = vector for slope (dt/dx)
% xleft = what is the x coordinate for the left-most input trace?
% term = number of term in complex number interpolation
%
% output
% ------
% tp = tau - p data
% xt_fk = fk data of the input xt data. f is in [0,fnqy] and k is in
%         [0,2*fnqy). This data set can be pushed back in cqfktpinv
%         program to be the background information in order to better
%         inverse to the original xt data.

% F-K Transform
[nt0,nx0] = size(xt);

nx = 2^nextpow2(nx0)*4; % 2 is required by the complex number interpolation
nt = 2^nextpow2(nt0);

xt_fk = fft2(xt,nt,nx);
xt_fk = xt_fk(1:nt/2+1,:);

knq = 1/2/dx;
fnq = 1/2/dt;
axis_f = linspace(0,fnq,nt/2+1);


% shift theorom for fk data
axis_k = linspace(0,2*knq,nx+1);
axis_k(end) = [];
kmat = repmat(axis_k,nt/2+1,1);
xt_fk = xt_fk.*exp( 1i*2*pi*kmat*-xleft );

% Use loops to calculate Tau - P data for each P
np = length(p);
fp = zeros(nt/2+1,np); % fp stands for 1D fft of tp in vertical direction
min_k = min(axis_k);
max_k = max(axis_k);
for n = 1:np
    % find searching range
    f_index = 1:size(fp,1);
    k_search = -p(n) * axis_f;
    n_rm = find(k_search < min_k | k_search > max_k);
    k_search(n_rm) = [];
    f_index(n_rm) = [];
    
    k_search = cqwrap2range(k_search, 0, 2*knq);
    
    for m = 1:length(k_search)
        this_ind = f_index(m);
        fp(this_ind,n) = ...
            cqcomplex_interp( xt_fk(this_ind,:), axis_k, ...
            k_search(m), term );
    end
end

% fft back to tau - p domain
fp = [real(fp);flipud(real(fp(2:end-1,:)))] + ...
    1i * [imag(fp);-flipud(imag(fp(2:end-1,:)))];
tp = real(ifft(fp,[],1));



end

