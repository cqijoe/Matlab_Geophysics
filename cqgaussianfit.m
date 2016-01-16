function [ fitpar,fitobj, gof, opts ] = ...
    cqgaussianfit( coldata,coltick,maxiter, maxfunevals )
% Fit each column of the input matrix with Gaussian curve defined as
% following:
%   * w = a*exp(-((x-b)/c)^2) ---- equation 1
% Program takes care of upper and lower boundaries and starting point
% Input: 
% * coldata ... data matrix where each column is a trace
% * coltick ... column tick. If each trace is wavelet, then coltick is
%               the time vector
% * maxiter ... maximum iteration time allowed (default=100)
% * maxfunevals ... maximum function evaluation times (default=100). 
%                  Total evaluation times for one trace:
%                           maxiter * maxfunevals
% Output:
% * fitpar ... fitting result structure
%   * a
%   * b
%   * c (refer to the equation 1)
% * fitobj ... fitting object store in a cell array
% * gof ... goodness of fitting cell array for each column in coldata
% * opts ... option structure used in the fitting

% convert coltick to column vector and check its length
coltick = coltick(:);
try
    if length(coltick)~=size(coldata,1)
        error(...
            'coltick length(%d) should be equal to coldata row size(%d)',...
            length(coltick),size(coldata,1));
    end
catch ME
    disp(ME.message);
    disp('Program stopped!')
    return;
end


% if maxiter is not provided then use 100
if ~exist('maxiter','var')
    maxiter = 100;
end

% if maxfuneval is not provided then use 100
if ~exist('maxfunevals','var')
    maxfunevals = 100;
end

% Set up fittype and options.
ft = fittype( 'gauss1' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Robust = 'LAR';
opts.MaxIter = maxiter;
opts.MaxFunEvals = maxfunevals;

% get some parameters
ncol = size(coldata,2);

% initialize output
fitpar = zeros(ncol,3);
fitobj = cell(ncol,1);
gof = fitobj;

% report progress every 1/4 of the process
cycrpt = round(ncol/4);
% fit 
disp('Start fitting ************')
% get time range for c_start
tmin = min(coltick);
tmax = max(coltick);
trg = tmax - tmin;
tstart = tic;
for iter = 1:ncol
    if mod(iter,cycrpt)==0
        fprintf('%d out of %d finished...elapse time %.1f seconds\n',iter,ncol,toc(tstart));
    end
    trace = coldata(:,iter);
    [a0,b0] = max(abs(trace));
    b0 = coltick(b0);
    c0 = 2.35/sqrt(2)*trg;
    opts.StartPoint = [a0,b0,c0];
    opts.Lower = [-2*a0,tmin,0];
    opts.Upper = [2*a0,tmax,inf];
    [fitobj{iter}, gof{iter}] = fit( coltick,trace, ft, opts );
    fitpar(iter,:) = ...
        [fitobj{iter}.a1, fitobj{iter}.b1, fitobj{iter}.c1];
end
fprintf('\n');
disp('Finish fitting *************')

end

