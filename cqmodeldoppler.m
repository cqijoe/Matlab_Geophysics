function [ s, sori, spilot,tr,t_pilot,t0 ] = ...
    cqmodeldoppler( xs0,xr0,zsrc,zrev,...
                    speed_boat, speed_rev,...
                    f_start,f_end,t_sweep,dt,t_total,taper,phase,...
                    v_medium,z_dfr,...
                    ifgeospread)
% A simple program to model doppler effects
% assume homogeneous lower half
% space. This assume a diffractor.
%
% input:
% ----------------Geometry part ---------------------------
% xs0 ... starting source x position
% xr0 ... starting reciever x position
% zsrc ... source z position, it will not change
% zrev ... reciever z position, it will not change when boat is moving
%   GEOMETRY DEFINITION ILLUSTRATION DIAGRAM
%   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~boat>>~~~~~~~~~~~~~~~~ <-- water surf
%   |<--xr------>|<---xs---->|
%   * receiver   |           |
%                |           * source
%      negative  |  positive
%                |
%                | <-- x = 0
%   -------------# <-- diffractor ---------------------
%
% ---------------Boat description -------------------------
% speed_boat ... speed of the boat
% speed_rev ... speed of the reviever.
%   (for OBC case set it to 0 and for streamer case with no relative 
%    moving set it to speed_boat)
% ---------------Source description ------------------------
% NOTE: WE ASSUME LINEAR SWEEP
% f_start ... starting sweep frequency
% f_end ... ending sweep frequency
% t_sweep ... total sweep time in seconds
% dt ... sampling rate.
%   NOTE: s WILL BE SAMPLED BY dt SAMPLING RATE TOO
% taper ... in seconds. The taper will be applied at both beginning
% and ending of the sweep. taper = 0 means no tapering
% phase ... in radiance illustrating the constant phase rotation
% --------------------Medium description -----------------------
% v_medium ... vp of the lower half space
% z_dfr ... depth of the refractor
% --------------------Phenomenon---------------------------------
% ifgeospread ... <true> to consider geometrical spreading
%
% Output
% s ... seismic trace at dt sampling rate recorded by receiver (correlated)
% sori ... original seismic trace recorded by receiver without
%    correlating the pilot sweep
% t0 ... first arrival time for v_boat = 0 case
% =============STEPS OF MODELLING============================
% 1. Determine maximum reciver recording time
% 2. Search for the source time ts for each of the reciever time tr:
%   a. calculate tu - upward raypath time by the equation
%       tu = sqrt( xr_tr.^2 + zr^2 )/v_medium
%       zr = z_dfr - zrev
%       xr_tr = xr0 + speed_rev*tr
%   b. calculate ts - source time corresponding to the reciever time
%      tr by the equation
%       (speed_boat^2+v_medium^2)*ts.^2 + ...
%       (2*xs0*speed_boat+2*(tr-tu)*v_medium^2)*ts + ...
%           xs0^2+zs^2-(tr-tu).^2*speed_boat^2;
%       zs = z_dfr - zsrc
%      A = (speed_boat^2+v_medium^2);
%      B = (2*xs0*speed_boat+2*(tr-tu)*speed_boat^2);
%      C = xs0^2+zs^2-(tr-tu).^2*speed_boat^2;
%      ts = [-B + sqrt(B^2-4*A*C)] ./ (2*A);
%      In the solution, only tr and tu are vectors with the same size.
%      So before calculation, make sure tr and tu are the same size.
%   c. calculate sample at ts using analytical expression for a linear
%      sweep. This can be vectorized using matlab built-in equation "chirp"
%   d. sori is derive using ts by chirp function
%   e. pillot trace spilot is calculated using tr by chirp function
%   f. crosscorrelate spilot with sori to get s

if ~exist('ifgeospread','var')||isempty(ifgeospread)
    ifgeospread = true;
end

% Define recording time
tr = 0:dt:t_total; tr = tr';

% Search for ts vector from tr
xr_tr = xr0 + speed_rev*tr; % receiver-side x coordinate change with time
zr = z_dfr - zrev; % refractor depth to receiver depth
Ru = sqrt(xr_tr.^2 + zr^2);
tu = Ru./v_medium; % upward raypath time

% Search for ts vector from tu
zs = z_dfr - zsrc;
A = (speed_boat^2 - v_medium^2);
B = (2*xs0*speed_boat + 2*(tr-tu)*v_medium^2);
C = xs0^2 + zs^2 - (tr-tu).^2*v_medium^2;
ts = real((sqrt(B.^2-4*A*C) - B)/(2*A));
xs_ts = xs0 + speed_boat*ts;
Rd = sqrt(xs_ts.^2 + zs.^2);

% Fill sori by cqchirp function
sori = cqchirp(ts,f_start,t_sweep,f_end,phase,taper);
if ifgeospread
    if ~any(Ru.*Rd==0)
        sori = sori.*(1./(Ru.*Rd));
    end
end

% Create pilot trace
t_pilot = 0:dt:t_sweep; t_pilot = t_pilot';
spilot = cqchirp(t_pilot,f_start,t_sweep,f_end,phase,taper);

% Cross correlate spilot with sori to get s
[s,lag] = xcorr(sori,spilot);
nt0 = find(lag==0);
s = s(nt0:nt0+length(sori)-1);

if nargout > 5
    % calculate t0
    t0 = interp1(ts,tr,0,'linear');
end

end

