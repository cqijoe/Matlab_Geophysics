function [ taup0,taup1 ] = cqlinradon( crp,dt,x,p,beta,v_boat,vib_info,dp,pfil,zo  )
% [ taup ] = cqlinradon( crp,dt,x,p,beta  )
% Slant stack using Linear Radon Transform. 
% Reference: Yilmaz's Text Book P988 (Seismic Data Analysis)
%
% input:
% crp ... [m x n] matrix where [sample number] = m and [trace number]=n
% dt ... sampling interval in second
% x ... offset vector. Length of x should be equal to column number of
%       crp
% p ... ray parameter vector. 
% beta ... damping factor (shoule be very small)
% v_boat ... (only for marine vibrator case) if this is given,
%                then program will do phase-correction using method by
%                Dragoset's geophysics paper
% vib_info ... [f_start,f_end,sweep_t] is 
%               starting sweep frequency
%               ending sweep frequency
%               sweeping time
% dp ... LMO for p in the input crp. So event with p1 is actually 
%        having slope of p1 + dp.
% pfil ... indexes indicating which p should be filtered and which
%          should not
% zo ... zeroing-out tau-p data out-side of pfil (true or false)
%
% If user don't want to do phase correction, then don't give
% ship_speed argument
% 
% output
% taup0 ... [m x k] matrix where [ray parameter number] = k
% taup1 ... same with taup0. only assigned when v_boat is given. it is
%           corresponded to the taup data after phase-correction
%           filtering

if nargin < 9
    pfil = 1:length(p);
end
if nargin < 10
    zo = true;
end

% check input
if length(x)~=size(crp,2)
    error('Invalid offset vector!');
end
h = x/2; % convert to half offset which is used by Yilmaz in his book

% make sure that we have power of 2 sample number
m0 = size(crp,1);
n = size(crp,2);
m = 2^nextpow2(m0);
crp = [crp;zeros(m-m0,n)]; % pad zeros to make crp have power of 2 samples

% calculate column-wise fft for crp to compose crp_w
crp_w = fft(crp,[],1);
df = 1/dt/m;
fnyq = 1/2/dt;
f = 0:df:fnyq; % frequency vector
w = 2*pi*f; 
crp_w = crp_w(1:m/2+1,:); % take only 0 ~ fnyq

% initialize taup_w which is column-wise 0 ~ fnyq tau-p domain data
taup_w = zeros(length(f),length(p));
taup_w = taup_w'; % so each column of taup_w represents one single w for all p

% loop into w to fill taup_w's each row for w
% transpose crp_w temporarily so each column is for one single w
crp_w = crp_w'; 
nh = length(h); % number of offset traces
np = length(p); % number of ray parameters
umat = diag(ones(np,1)); % unit matrix that would be used by eq F-23
nw = length(w);
for k = 1:length(w)
    [ s ] = cqnotify( k,nw, 25 );
    if ~isempty(s)
        disp(['cqlinradon:',s])
    end
    % compose L matrix
    L = zeros(nh,np);
    for iter_h = 1:nh
        for iter_p = 1:np
            % Yilmaz use - in front of 1i and I try not using this
            % rule, plz review the same L definition in cqinvlinradon
            % function!
            L(iter_h,iter_p) = exp(1i*w(k)*2*p(iter_p)*h(iter_h));
        end
    end
    % svd the L
    [U,S,V] = svd(L);
    % using equation F-23 in Yilmaz's book (p985 vol 1)
    taup_w(:,k) = V*((S'*S + beta*umat)\S')*U'*crp_w(:,k);
end


% inverse taup_w back to p - t domain, now it is p - w domain
taup_w = taup_w'; % transpose back so each column is for a single p

% recompose two-sided freq domain taup data
real_part = real(taup_w);
img_part = imag(taup_w);
taup_w0 = [real_part;flipud(real_part(2:end-1,:))] + ...
    1i * [img_part;-flipud(img_part(2:end-1,:))];
taup0 = real(ifft(taup_w0,[],1));

% if ship_speed is given, then apply Dragoset's method for phase
% correction filtering
if exist('v_boat','var')
    p_doppler = p - dp;
    dopfactor = -v_boat * p_doppler;
    f0 = vib_info(1); % starting freq
    f1 = vib_info(2); % ending freq
    finterval = f1 - f0;
    tswp = vib_info(3); % sweeping time
    nf = find(f>=f0 & f<=f1); % sample number in the taup_w needs for correction
    phase = zeros(length(f),np); % phase correction
    
    % for every ray paramter, calculate a column of phase
    for k = 1:length(p)
        % decide whether to filer this p
        if ismember(k,pfil)
            phase(nf,k) = -2*pi*dopfactor(k)*tswp*f(nf).^2./finterval;
        else
            if zo
            % test: set taup data 0 outside the selection of p
            taup_w(:,k) = 0;
            end
        end
    end
    
    % filter out everything outside the frequency range
       taup_w([1:nf(1),nf(end):end],:) = 0;
    
    % update taup_w by phase correction filter
        taup_w = taup_w.*exp(1i*phase);
    
    % recompose two-sided freq domain taup data
    real_part = real(taup_w);
    img_part = imag(taup_w);
    taup_w = [real_part;flipud(real_part(2:end-1,:))] + ...
        1i * [img_part;-flipud(img_part(2:end-1,:))];
    taup1 = real(ifft(taup_w,[],1));
end




end

