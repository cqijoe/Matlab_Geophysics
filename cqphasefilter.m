function [ pf, t, PF, F ] = ...
    cqphasefilter( p, dt, N,...
                   f_start, f_end, t_sweep,...
                   v_boat, f_slope)
% the function calculates phase filter in frequency domain and in time
% domain for Doppler effects
%
% input:
% ------------ Requirement ---------------------------
% p ... [1 x m] dt/dx of the trace
% dt ... sampling rate of seconds for the desired phase filter
% N ... number of samples in the output filter, must be power of 2 or
% program will change N to the next power of 2
% ------------ Source Description ---------------------
% f_start ... start sweeping frequency in Hz
% f_end ... end sweeping frequency in Hz
% t_sweep ... sweeping time from f_start to f_end in seconds
% ------------- Other Parameters -----------------------
% v_boat ... velocity of the boat
% f_slope ... slope of frequency to avoid Gibs effects (hz)
%
% output:
% -------------- Time Domain Filter --------------------
% pf ... [N x m]time domain phase filter. it has its time zero in the middle
% t ... [N x 1] corresponding time vector for pf filter
% -------------- Frequency Domain Filter ----------------
% PF ... [N x m] frequency domain phase filter. N = length(F)
% F ... [N x 1] corresponding frequency vector for PF filter.
% 
% Algorithm
% 1. Get f vector from N and dt
% 2. Calculate Doppler Factor for each p
% 3. Calculate phase spectrum for each f from f_start to f_end
% 4. Fill PF
% 5. Inverse fft to get pf

% checking N and make necesary change to the next power of 2
if N ~= 2^nextpow2(N)
    N = 2^nextpow2(N);
end

% preparing parameters
p = reshape(p,1,length(p));
m = length(p);
df = 1/dt/N;
fnyq = 1/2/dt;
f = (0:df:fnyq)'; % 0 ~ Nyquist Frequency
t = (0:dt:N/2*dt)';
t = [t;-flipud(t(2:end-1))];
t = [t(N/2+2:end);t(1:N/2+1)]; % t is arranged from -t1 to t2
F = [f;-flipud(f(2:end-1))];
% F = [F(N/2+2:end);F(1:N/2+1)]; % F is arranged from -f1 to f2

% calculate Doppler Factors
dopfact = -v_boat * p;
amp = zeros(length(f),m); % amplitude spectrum of the filter for 0 ~ Nyquist Frequency
phase = amp;
nf = find(f>=f_start&f<=f_end);
finterval = f_end - f_start;
f_start_end = f_start + f_slope;
f_end_begin = f_end - f_slope;
% calculate amplitude and phase spectrum for filters for each p
for iter = 1:length(p)
    amp(nf,iter) = 1;
    if f_slope~=0
        fvec = f(f>=f_start&f<=f_start_end);
        amp(f>=f_start&f<=f_start_end,iter) = interp1([f_start;f_start_end],[0;1],fvec,'linear');
        fvec = f(f>=f_end_begin&f<=f_end);
        amp(f>=f_end_begin&f<=f_end,iter) = interp1([f_end_begin;f_end],[1;0],fvec,'linear');
    end
    phase(nf,iter) = -2*pi*dopfact(iter) * t_sweep * f(nf).^2./finterval;
end
amp = [amp;flipud(amp(2:end-1,:))];
phase = [phase;-flipud(phase(2:end-1,:))];
PF = amp.*exp(1i*phase);

pf = ifft(PF,[],1);
pf = [pf(N/2+2:end,:);pf(1:N/2+1,:)];

end

