function [ crp ] = cqinvlinradon( taup,dt,x,p )
% [ crp ] = cqinvlinradon( taup,dt,x,p )
% Inverse linear radon transform. It simply apply phase shift filter
% to each w trace in taup data, i.e. d = Lu for each w
% 
% input
% taup ... taup domain data in column (each column is for a single p)
% dt ... sampling rate in time in seconds
% x ... offset vector for the output crp
% p ... ray parameter vector corresponding to the input taup matrix
% note: length(p) should equal to size(taup,2)
%
% output
% crp ... common receiver gather in column (this could be common midpoint gather)
%         or common shot gather

% check input
if length(p)~=size(taup,2)
    error('Invalid ray parameter vector!');
end
h = x/2; % convert to half offset which is used by Yilmaz in his book

% make sure that we have power of 2 sample number
m0 = size(taup,1);
n = size(taup,2);
m = 2^nextpow2(m0);
taup = [taup;zeros(m-m0,n)]; % pad zeros to make taup have power of 2 samples

% calculate column-wise fft for taup to compose taup_w
taup_w = fft(taup,[],1);
df = 1/dt/m;
fnyq = 1/2/dt;
f = 0:df:fnyq; % frequency vector
w = 2*pi*f; 
taup_w = taup_w(1:m/2+1,:); % take only 0 ~ fnyq

% initialize crp_w which is column-wise 0 ~ fnyq x - t domain data
crp_w = zeros(length(f),length(h));
crp_w = crp_w'; % so each column of crp_w represents one single w for all h

% loop into w to fill crp_w's each row for w
% transpose taup_w temporarily so each column is for one single w
taup_w = taup_w'; 
nh = length(h); % number of offset traces
np = length(p); % number of ray parameters
for k = 1:length(w)
    % compose L matrix
    L = zeros(nh,np);
    for iter_h = 1:nh
        for iter_p = 1:np
            L(iter_h,iter_p) = exp(1i*2*w(k)*p(iter_p)*h(iter_h));
        end
    end
    % apply d = Lu for crp_w
    crp_w(:,k) = L * taup_w(:,k);
end

% inverse taup_w back to p - t domain, now it is p - w domain
crp_w = crp_w'; % transpose back so each column is for a single p
% recompose two-sided freq domain taup data
real_part = real(crp_w);
img_part = imag(crp_w);
crp_w = [real_part;flipud(real_part(2:end-1,:))] + ...
    1i * [img_part;-flipud(img_part(2:end-1,:))];
crp = real(ifft(crp_w,[],1));

end

