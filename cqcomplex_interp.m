function [ coeffq ] = cqcomplex_interp( coeff, dt, fq, term )
% complex number interpolation using Fourier series
% method referenced to Wade's UH Thesis, 1989
%
% input
% -----
% coeff = Fourier coefficient with length 2 * N, N must be integer and
%         base of 2. The coeff is arranged in [0 ~ 2*pi) using matlab
%         fft method. (DON'T use fftshift)
% dt = sampling rate before fourier transform to get coeff
% fq = query frequency (df = 1/dt)
%      (fq must be lying within -1/2/dt to 1/2/dt)
% term = 6 = how many terms you want to use to interpolate
% 
% output
% ------
% coeffq = complex value at the query frequency, it is a column vector
%
% Papare coeff
% ------------
% 1. time_series is a vector with length N
% 2. coeff = fft(time_series,2*N)
% 3. The principle is to make sure the coeff's length must be two
%    times larger than the useful data length in time_series

if ~exist('term','var')||isempty(term)
    term = 6;
end

% check input
N_2 = length(coeff);
N = N_2/2;
if mod(N_2,2)~=0
    error('Invalud coeff, its length should be even!');
end

fnyq = 1/2/dt;
if max(fq) > fnyq || min(fq) < -fnyq
    error('Invalid query frequency, it should not out-bound Nyquist!');
end

fq = fq(:);
% recompose input data from [0,2*pi) to [-pi,pi]
coeff = coeff(:); % convert to column vector
% now coeff is from [-pi ~ pi]
coeff = [coeff(end);coeff(N_2/2+2:end);coeff(1:N_2/2+1)];
sgm = fq * 2 * N * dt;

% find nearest term number of data to interpolate and do interpolation
% point-by-point
gma = -N:N;
gma = gma(:); % integer index from 0 to 2N
coeffq = zeros(length(fq),1);
for k = 1:length(fq)
    % sort the difference and select the minimum term
    [~,I] = sort(abs(gma-sgm(k)));
    I = I(1:term);
    gma_sgm = (gma(I) - sgm(k)); % gamma minus sigma times pi
    coeffq(k) = ...
        sum( coeff(I) .* exp(1i*gma_sgm/2).*...
        sinc(gma_sgm/4).^2 .* sinc(gma_sgm) );
    % Matlab sinc function: sin(pi*x)/(pi*x)
end



end

