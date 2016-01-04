function [ d2 ] = cqlinshift( d, dt, t0, intpmth )
% shift traces in d by t0 vector
% 
% input
% -----
% d = 2D input matrix traces are columns
% dt = sampling interval in seconds
% t0 = vector of time shift for each trace
%      positve for down shifting and negative for up shifting
% intpmth = 'linear' 
%            other interpolation methods are:
%            'linear'   - (default) linear interpolation
%            'nearest'  - nearest neighbor interpolation
%            'next'     - next neighbor interpolation
%            'previous' - previous neighbor interpolation
%            'spline'   - piecewise cubic spline interpolation (SPLINE)
%            'pchip'    - shape-preserving piecewise cubic interpolation
%            'cubic'    - same as 'pchip'
%            'v5cubic'  - the cubic interpolation from MATLAB 5, which does not
%                         extrapolate and uses 'spline' if X is not equally
%                         spaced.
%
% output
% ------
% d2 = 2D output matrix after shifting

if nargin < 4
    intpmth = 'linear';
end

% check input
if isscalar(t0)
    t0 = repmat(t0,1,size(d,2));
else
    if length(t0)~=size(d,2)
        error('Invalud t0 vector size!');
    end
end

d2 = zeros(size(d));
t = 0:dt:(size(d,1)-1)*dt; 
t = t(:); % reshape t to column so output will be column
for k = 1:size(d,2)
    tnew = t - t0(k);
    d2(:,k) = interp1(t,d(:,k),tnew,intpmth,0); % extrapolate with 0
end








end

