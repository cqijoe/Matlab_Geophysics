function [ data2,data_p ] = cqslantstack( data,x,t, p,tau,ifrho,str_interp1 )
% The function will do the slant stacking for the input data matrix
% Input:
% data ... 2D matrix, x is distance and y is time
% x ... distance coordinates
% t ... time coordinates for the data, x and t must match data size in
%       the first and second axis direction
% p ... the vector for the p values used in tau-p transform
% note: data2 will be in the matrix with horizontal axis p and
%       vertical axis t
% tau ... by default, tau equals to t but it can change
% ifrho ... true to apply rho filter before stacking <true>
% str_interp1 ... the method passed to the interp1 equation
%       'linear'   - (default) linear interpolation
%       'nearest'  - nearest neighbor interpolation
%       'next'     - next neighbor interpolation
%       'previous' - previous neighbor interpolation
%       'spline'   - piecewise cubic spline interpolation (SPLINE)
%       'pchip'    - shape-preserving piecewise cubic interpolation
%       'cubic'    - same as 'pchip'
%       'v5cubic'  - the cubic interpolation from MATLAB 5, which does not
%                    extrapolate and uses 'spline' if X is not equally
%                    spaced.

if ~exist('ifrho','var')||isempty(ifrho)
    ifrho = true;
end
if ~exist('str_interp1','var')||isempty(str_interp1)
    str_interp1 = 'linear';
end
% 1. initialize data2
data2 = zeros(length(tau),length(p));
% 2. loop first in x then in tau
if ~exist('tau','var')||isempty(tau)
    tau = t(:);
end
tau = tau(:);
t = t(:);
x = reshape(x,1,length(x));
p = p(:);

taumat = repmat(tau,1,length(x));
xmat = repmat(x,length(tau),1);

% >>>>>>>>>>>>>>>>>>>>>>>>>..
for iter1 = 1:length(p)
    s = cqnotify(iter1,length(p),25);
    if ~isempty(s)
        disp(['cqslantstack:',s]);
    end
    data_taup = zeros(length(tau),size(data,2));
    tsearch = taumat + p(iter1)*xmat;
    for iter2 = 1:length(x)
        data_taup(:,iter2) = ...
            interp1(t,data(:,iter2),tsearch(:,iter2),str_interp1,0);
    end
    data2(:,iter1) = sum(data_taup,2);
    if nargout>1
        data_p(:,:,iter1) = data_taup;
    end
end
data2 = data2/length(p);
% do rho filter
if ifrho
    dtau = tau(2) - tau(1); 
    Nrho = 2001; % Nrho points tau-domain filter
    % design rho filter
    rhofilter = zeros(length(Nrho),1);
    nrho = 1;
    for iterj = -(Nrho-1)/2:(Nrho-1)/2
        if iterj~=0&&( mod(iterj,2)==0 )
            rhofilter(nrho) = -4/dtau^2/iterj^2;
        end
        nrho = nrho + 1;
    end
    rhofilter((Nrho+1)/2) = -sum(rhofilter);
%     rhofilter((Nrho+1)/2) = (2*pi^2)/(dtau^2);
    % >>>>>>>>>>>>>>>>>>>>>>
    % apply rho filter do tau-p domain date trace by trace
    for iterp = 1:length(p)
        data2(:,iterp) = cqi_conv(data2(:,iterp),rhofilter,(Nrho+1)/2);
    end
end


end

