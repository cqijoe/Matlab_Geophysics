function [ gather ] = cqdpflt( xs0, xr0, us, ur, z, ...
    v, f0, f1, tswp, taper, tuncor, tcor, dt )
% Model Doppler effects for a flat reflector
% 
% input
% -----
% xs0 = vector for source locations
% xr0 = vector for receiver locations
% *****
% program determine gather type based on xs0 and xr0 size
% If ns = 1, program assumes to generate common-shot-gather
% If nr = 1, program assumes to generate common-receiver-gather
% If ns = nr, program assumes to generate common-midpoint-gather
% *****
% us = scalor for source speed (m/s)
% ur = scalor for receiver speed. Both speeds are constants. (m/s)
% z = vector for reflector depth 
% v = scalor for medium (usually) water speed
% f0 = starting sweeping frequency in Hz
% f1 = ending sweeping frequency in Hz
% tswp = sweeping time length in seconds
% taper = taper at both ends in seconds
% tuncor = maximum time in uncorrelated gather in seconds
% tcor = maximum time in correlated gather in seconds
% dt = sampling rate in seconds
%
% output
% ------
% gather = structure with output information
%    dcor = correlated gather
%    duncor = uncorrelated gather

% check input
ns = length(xs0);
nr = length(xr0);
if ns~=nr && nr~=1 && ns~=1
    error('Invalid xs0 and xr0 length!');
end

t = 0:dt:tuncor; % get recording time 
t = t(:);
nsamp = length(t); % total samples per trace
tuncor = t(end); % in case tuncor is not an integer times dt

% Create pilot trace
t_pilot = 0:dt:tswp; t_pilot = t_pilot(:);
tr_plt = cqchirp(t_pilot,f0,tswp,f1,0,taper);

t_align = cell(length(z),1);

if ns == nr
    dcor = zeros(nsamp,ns);
    gather.type = 'common-midpoint-gather';
elseif ns == 1
    dcor = zeros(nsamp,nr);
    gather.type = 'common-shot-gather';
else
    dcor = zeros(nsamp,ns);
    gather.type = 'common-receiver-gather';
end
duncor = dcor;

for nz = 1:length(z)
if ns == nr
    % create common-midpoint-gather
    for k = 1:ns
        ts = search(xs0(k),xr0(k),z(nz));
        % create uncorrelated records
        duncor(:,k) = duncor(:,k) + cqchirp(ts,f0,tswp,f1,0,taper);
        t_align{nz}(k) = interp1(ts, t, 0, 'linear', 'extrap');
    end
    
elseif ns == 1
    % create common-shot-gather
    for k = 1:nr
        ts = search(xs0,xr0(k),z(nz));
        % create uncorrelated records
        duncor(:,k) = duncor(:,k) + cqchirp(ts,f0,tswp,f1,0,taper);
        t_align{nz}(k) = interp1(ts, t, 0, 'linear', 'extrap');
    end
    
else
    for k = 1:ns
        ts = search(xs0(k),xr0,z(nz));
        % create uncorrelated records
        duncor(:,k) = duncor(:,k) + cqchirp(ts,f0,tswp,f1,0,taper);
        t_align{nz}(k) = interp1(ts, t, 0, 'linear', 'extrap');
    end
end
end

for ncorr = 1:size(duncor,2)
    % create correlated records
    [s,lag] = xcorr(duncor(:,ncorr),tr_plt);
    nt0 = find(lag==0);
    dcor(:,ncorr) = s(nt0:nt0+nsamp-1);
end

dcor = dcor(1:floor(tcor/dt),:);


gather.data_correlate = dcor;
gather.data_uncorrelate = duncor;
gather.pilot_trace = tr_plt;
gather.dt = dt;
gather.source_speed = us;
gather.receiver_speed = ur;
gather.source_x = xs0;
gather.receiver_x = xr0;
gather.freq_low = f0;
gather.freq_high = f1;
gather.taper = taper;
gather.sweeping_time = tswp;
gather.reflector_depth = z;
gather.nmo_time = t_align;
gather.tuncor = tuncor;
gather.tcor = tcor;


% =================
    function ts = search(xs0, xr0,z)
        % search for source time
        % xs0 and xr0 are only scalar
        D = xr0 + t*ur - xs0 - t*us; % D is broadcast by tuncor
        A = us^2 - v^2;
        B = 2 * D * us; % broadcast
        C = D.^2 + 4*z^2; % broadcast
        ts = t - (-B - sqrt(B.^2-4*A*C))/(2*A);
    end


end

