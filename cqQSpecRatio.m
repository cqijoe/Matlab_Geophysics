function [ Q, A0, A, f ] = cqQSpecRatio( d0, t0, d, t, dt, fmax, QC )
% Calculate Q from Spectral Ratio Method
% 
% input
% -----
% d0 = trace at reference time t0
% t0 = reference time for d0
% d = traces in matrix at referece time vector t
% t = reference time vectorfor d
% dt = sampling rate for d0 and d
% fmax = maximum fitting frequency || [fmin, fmax]
%        in the first case, program finds max amplitude f and use the
%        searched f as fmin
%        
% QC = false
%      open figures to check fitting result and Q extraction
%
% output
% ------
% Q = estimated Q for traces in d at t
% A0 = amplitude spectrum at t0
% A = amplitude spectrum at t
% f = reference frequency vector

% input check
if size(d,2) ~= length(t)
    error('t and d are not matching!');
end

if ~exist('QC','var')||isempty(QC)
    QC = false;
end

% fft d0 and d and get reference frequency

% pad d0 and d to be the same size
d0 = d0(:);
if length(d0)~=size(d,1)
    dn = length(d0) - size(d,1);
    if dn > 0
        d = [d;zeros(dn,size(d,2))];
    else
        d0 = [d0;zeros(-dn,1)];
    end
end

% fft 
N = 2^nextpow2(length(d0));
d0 = fft(d0,N); d0 = d0(1:N/2+1); % take 0 ~ fnyq only
d = fft(d,N,1); d = d(1:N/2+1,:);
f = linspace(0,1/2/dt,N/2+1); % reference frequency 
f = f(:);

if length(fmax)==2
    fmin = fmax(1);
    fmax0 = fmax(2);
    [~,nmax] = min(abs(f-fmin)); % nmax actually is nfmin
    nmax = ceil(nmax);
else
    fmax0 = fmax;
end
[~,nfmax] = min(abs(f-fmax0));
nfmax = floor(nfmax);

A0 = abs(d0);
A = abs(d); % take amplitude spectrum

for k = 1:size(d,2)
    w = 2 * pi * f;
    c = log(A(:,k)./A0);
    if length(fmax)==1
        [~,nmax] = max(c); % find dominant frequency
    end
    k0 = polyfit(w(nmax:nfmax),c(nmax:nfmax),1);
    Q(k) = - (t(k)-t0)/2/k0(1);
    
    if QC
        figure;
        plot(f,c,'k','linew',2);
        hold on
        plot(f(nmax:nfmax),polyval(k0,w(nmax:nfmax)),'r','linew',2);
        xlabel('Frequency(Hz)'); ylabel('ln{A0/A}');
        title(['Trace = ',num2str(k),', Q = ',num2str(Q(k))]);
    end
end





end

