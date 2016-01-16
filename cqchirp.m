function x = cqchirp(t,f0,t1,f1,phase,taper)
%x = cqchirp(t,f0,t1,f1) generates samples of a linear swept-frequency
%   signal at the time instances defined in timebase array t.  The instantaneous
%   frequency at time 0 is f0 Hertz.  The instantaneous frequency f1
%   is achieved at time t1.
%   The argument 'phase' is optional. It defines the initial phase of the 
%   signal degined in radians. By default phase=0 radian
% 
% NOTE: WHEN T IS OUTSIDE [0,T1], X = 0
% 
% input
% t ... time vector
% f0 ... frequency at time 0
% t1 ... sweep time from time 0 to time t1
% f1 ... frequency at time t1
% phase ... phase of the sweep in radians
% taper ... time of taper. if taper is scalar then it is applied to
% start and end of the sweep; if taper is [a,b] then a is applied for
% the beginning of the sweep and b is applied for the ending of the
% sweep.

 
    
if nargin==4
    phase=0;
end
if ~exist('taper','var')||isempty(taper)
    taper = [0,0];
end
if isscalar(taper)
    taper = [taper,taper];
end
t = t(:); % making t a column vector
k=(f1-f0)/t1;
x = zeros(length(t),1);
nval = find(t>=0&t<=t1);
x(nval) = cos(2*pi*(k/2*t(nval)+f0).*t(nval)+phase);

% apply taper to the sweep start
if taper(1)~=0
    ngain = find(t>=0&t<=taper(1)); % sample for applying gain filter
    valgain = interp1([0;taper(1)],[0;1],t(ngain),'linear'); % interpolate to get gain value at sampled time
    x(ngain) = x(ngain).*valgain;
end

% apply taper to the sweep end
if taper(2)~=0
    ngain = find(t>=(t1-taper(2))&t<=t1); % sample for applying gain filter
    valgain = interp1([t1-taper(2);t1],[1;0],t(ngain),'linear'); % interpolate to get gain value at sampled time
    x(ngain) = x(ngain).*valgain;
end

end