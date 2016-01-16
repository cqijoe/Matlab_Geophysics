function [ data,data_x ] = cqinvslantstack( data2,x,t,p,tau,ifrho,str_interp1 )
% inverse tau-p transform
%
% Input:
% data2 ... 2D matrix, horizontal is p and vertical is tau
% p ... the vector for the p values used in tau-p transform
% tau ... tau should be the same length as size(data2,1)
% x ... distance coordinates
% t ... time coordinates for the data, x and t must match data size in
%       the first and second axis direction (by default, t=tau)
% ifrho ... <true> to apply rho filter before stacking
% Output:
% data ... 2D matrix, horizontal is distance (x) and vertical is time (t)
% data_x ... if provided, program outputs LMO tau-p data to data_x(:,:,nx)
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

if ~exist('str_interp1','var')||isempty(str_interp1)
    str_interp1 = 'linear';
end
% data condition
p = reshape(p,1,length(p));
tau = tau(:);
x = reshape(x,1,length(x));
if ~exist('t','var')
    t = tau;
end
t = t(:);
if ~exist('ifrho','var')
    ifrho = true;
end
% 1. initialize data
data = zeros(length(t),length(x));
% 3. do inverse tau-p transform
tmat = repmat(t,1,length(p));
pmat = repmat(p,length(t),1);

for iter1 = 1:length(x)
    s = cqnotify(iter1,length(x),25);
    if ~isempty(s)
        disp(['cqinvslantstack:',s]);
    end
    data_taup = zeros(length(t),length(p));
    tausearch = tmat - x(iter1)*pmat;
    for iter2 = 1:length(p)
        data_taup(:,iter2) = ...
            interp1(tau,data2(:,iter2),tausearch(:,iter2),str_interp1,0);
    end
    data(:,iter1) = sum(data_taup,2)/size(data_taup,2);
    if nargout>1
        data_x(:,:,iter1) = data_taup;
    end
end
data = data/length(x);
% do rho filter
if ifrho
    dt = t(2) - t(1); 
    Nrho = 2001; % Nrho points tau-domain filter
    % design rho filter
    rhofilter = zeros(length(Nrho),1);
    nrho = 1;
    for iterj = -(Nrho-1)/2:(Nrho-1)/2
        if iterj~=0&&( mod(iterj,2)==0 )
            rhofilter(nrho) = -4/dt^2/iterj^2;
        end
        nrho = nrho + 1;
    end
    rhofilter((Nrho+1)/2) = -sum(rhofilter);
%     rhofilter((Nrho+1)/2) = (2*pi^2)/(dtau^2);
    % >>>>>>>>>>>>>>>>>>>>>>
    % apply rho filter do tau-p domain date trace by trace
    for iterp = 1:length(x)
        data(:,iterp) = cqi_conv(data(:,iterp),rhofilter,(Nrho+1)/2);
    end
end


end

