function [ d0 ] = cqslantshift( d, u, dt, dx )
% Slant shift seismic data in wavenumber domain
%
% input
% -----
% d = common-shot or common-receiver-gather
% u = receiver (for common-shot-gather) or shot (for common-receiver-gather)
%     speed
% dt = sampling rate in second
% dx = sampling rate in x
%
% output
% ------
% d0 = d after slant shifting

t = 0:dt:(size(d,1)-1)*dt; t = t(:);
shift = u*t;
npad = round(max(shift)/dx*2);
ntr0 = size(d,2);
if u>0
   d = [d,zeros(size(d,1),npad)]; 
else
    d = [zeros(size(d,1),npad),d];
end

ntr = size(d,2);
N = 2^nextpow2(ntr);
d = fft(d,N,2);
d = d(:,1:N/2+1);
k = linspace(0,1/2/dx,N/2+1); % sampling wavenumber
k = repmat(k,size(d,1),1);


shift = repmat(shift,1,N/2+1);
p = exp(-1i*2*pi.*k.*shift);

d0 = d.*p;
d0 = [real(d0),fliplr(real(d0(:,2:end-1)))] + ...
    1i* [imag(d0),-fliplr(imag(d0(:,2:end-1)))];

d0 = real(ifft(d0,[],2));
d0 = d0(:,1:ntr0);

end

