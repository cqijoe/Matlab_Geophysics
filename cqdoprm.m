function [ d2, d1, d0 ] = cqdoprm( d, dt, x, p, v, f0, f1, tswp, prng,trcut )
% The function remove doppler effects in x-f domain.
%
% input
% -----
% d = 2D matrix of data
% dt = sampling interval in seconds
% x = offset vector
% p = slope vector
% v = boat speed
% f0 = start frequency
% f1 = end frequency
% tswp = sweeping duration in seconds
% prng = 1:length(p)
%        indexes of p to filter
% trcut = <false>
%         whether to truncate Tau - P Data by hand after forward Tau -
%         P transform
% output
% ------
% d2 = d after doppler effect removal
% d1 = Optional output for tau - p data after filtering
% d0 = Optional output for tau - p data before filtering

% check input
if length(x)~=size(d,2)
    error('x size invalid!');
end
if ~exist('prng','var')||isempty(prng)
    prng = 1:length(p);
end
if ~exist('trcut','var')||isempty(trcut)
    trcut = false;
end
ns = size(d,1);
% x-t to x-f transform
N = 2^nextpow2(size(d,1));
d = fft(d,N,1);

% only take 0 ~ fnyq
d = d(1:N/2+1,:);
f = linspace(0,1/2/dt,N/2+1);
f = f(:);

% initialize output
d0 = zeros(size(d,1), length(p));
d1 = d0; % p-f data
d2 = zeros(size(d)); % x-f data

% shape x
x = x(:)';

% forward tau-p and phase correction filtering
nf = find(f>=f0 & f<=f1);
f01 = f1 - f0; % freq interval
% loop p
for k = 1:length(p)
    % tp transform
    dtx = p(k) * x;
    pshift = exp(1i*2*pi*f*dtx);
    trf = sum(d.*pshift,2); 
    d0(:,k) = trf;
end
d0 = d0/length(p);
if trcut
    % inverse transform p-f to p-tau domain
    rld0 = real(d0); imd0 = imag(d0);
    rld0 = [rld0;flipud(rld0(2:end-1,:))];
    imd0 = [imd0;-flipud(imd0(2:end-1,:))];
    d0 = rld0 + 1i*imd0;
    d0 = real(ifft(d0,[],1));

    % plot to let user truncate artifacts
    d0 = cqtrcut(d0);
    
    % fft d0
    d0 = fft(d0,[],1);
    d0 = d0(1:N/2+1,:);
end
for k = 1:length(p)
    trf = d0(:,k);
    if ismember(k, prng)
        % phase-correction filtering
        dpfct = -v * p(k);
        phc = zeros(size(trf));
        phc(nf) = -2*pi*dpfct*tswp*f(nf).^2./f01;
        trf = trf.*exp(1i*phc);
    end
    d1(:,k) = trf;
end

% inverse tau-p
for k = 1:length(x)
    dtx = x(k) * p;
    pshift = exp(-1i*2*pi*f*dtx);
    trf = f.*sum(d1.*pshift,2);
    d2(:,k) = trf;
end

% ifft
rld2 = real(d2); imd2 = imag(d2);
rld2 = [rld2;flipud(rld2(2:end-1,:))];
imd2 = [imd2;-flipud(imd2(2:end-1,:))];
d2 = rld2 + 1i*imd2;
d2 = real(ifft(d2,[],1));
d2 = d2(1:ns,:);
d2 = pi/2*d2/length(x);

% optional output
if nargout > 1
    rld1 = real(d1); imd1 = imag(d1);
    rld1 = [rld1;flipud(rld1(2:end-1,:))];
    imd1 = [imd1;-flipud(imd1(2:end-1,:))];
    d1 = rld1 + 1i*imd1;
    d1 = real(ifft(d1,[],1));
    d1 = pi/2*d1/length(x);
    if nargout > 2
        rld0 = real(d0); imd0 = imag(d0);
        rld0 = [rld0;flipud(rld0(2:end-1,:))];
        imd0 = [imd0;-flipud(imd0(2:end-1,:))];
        d0 = rld0 + 1i*imd0;
        d0 = real(ifft(d0,[],1));
        d0 = pi/2*d0/length(x);
    end
end

end

