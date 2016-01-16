%% Read Las File
% Picking Coal and Shale
% Monte Carlo Simulation
cqi.coalmodel;
%% Load data from previous section
load('C:\Users\Qi\Documents\MATLAB\cqi\cqsFreqAndDelay\wlog_mc.mat');
%% Test On ones Well for Kernel Fitting: Decide to use Gaussion Kernel
%% Generate Synthetic and Wavefield Dictionary
imp = wlog_mc.vp_model_time(:,end).*wlog_mc.rhob_model_time(:,end);
rc = cqi.imp2rc(imp);
pm = cqi.pm(rc);
[~,wd] = cqi.reflectivity(rc,0,200);
%% 
coltick = 0:0.001:(size(wd,1)-1)*0.001;
% use Gaussian kernel to fit
[ fitpar,fitobj, gof,opt_fit ] = cqgaussianfit( ...
    wd,coltick, 20,20);
%% Check goodness of fitting
rmse = zeros(length(gof),1);
for iter = 1:length(gof)
    rmse(iter) = gof{iter}.rmse;
end
hold on;
plot(rmse,'r');
%% Visualize the difference
wd_fit = zeros(size(wd));
for iter = 1:size(wd,2)
    wd_fit(:,iter) = fitobj{iter}(coltick(:));
end
%% Plot to check fitting result
% generate sampled data curves
wd_fit = zeros(size(wd));
for iter = 1:size(wd,2)
    wd_fit(:,iter) = fitobj{iter}(coltick);
end

% generate cross-plot
figure;
scatter(fitpar(:,2),fitpar(:,3))
xlabel('b')
ylabel('c')
title('b-c plot')
%% The result looks good however the searching takes time
% In the loop over each synthetic well, the wavefield is sparsed
imp_model = wlog_mc.vp_model_time.*wlog_mc.rhob_model_time;
rc_model = cqi.imp2rc(imp_model);
tolmod = size(imp_model,2);
notify_interval = tolmod/10;
t_start = tic;
nskip = 10;
nwf = 100; % length of wavefield
coltick = 0:0.001:(nwf-1)*0.001;
fitpar_model = cell(tolmod,1);
gof_model = cell(tolmod,1);
%% For loop takes time, do execute it if not necessary
for iter = 1:size(imp_model,2);
    clc;
    % notify
    if mod(iter,notify_interval)==0
        disp('********************************');
        fprintf('%d outof %d model is inverted...elapse time %f\n',...
            iter,tolmod,toc(t_start));
    end
    pm = cqi.pm(rc_model(:,iter));
    [~,wd] = cqi.reflectivity(rc_model(:,iter),0,nwf);
    wd = wd(:,1:nskip:end);
    % use Gaussian kernel to fit
    [ fitpar_model{iter},~,gof_model{iter} ] = cqgaussianfit( ...
        wd,coltick, 50,50);
end
%% Cross-plot all fitpar_model par2 and par3 (b and c) (time shift and broadness)
figure; hold on;
for iter = 1:tolmod
    x = fitpar_model{iter}(:,3);
    x(x>20)=0;
    fitpar_model{iter}(:,3)=x;
    scatter(fitpar_model{iter}(:,2),fitpar_model{iter}(:,3));
end
%% Convert the cell array into matrix
wfpermod = size(fitpar_model{1},1); % wavefield per model
fitpar_matrix = zeros(tolmod*wfpermod,3);
for iter = 1:tolmod
    indx = (iter-1)*wfpermod+1:iter*wfpermod;
    for itercol = 1:3
        fitpar_matrix(indx,itercol) = fitpar_model{iter}(:,itercol);
    end
end
%%
xsmp = 50;
ysmp = 50;
[ f_xy ] = cqi_bidistr( fitpar_matrix(:,2),fitpar_matrix(:,3),...
    linspace(min(fitpar_matrix(:,2)),...
    max(fitpar_matrix(:,2)),xsmp),...
    linspace(min(fitpar_matrix(:,3)),...
    max(fitpar_matrix(:,3)),ysmp),...
    true,'x_name','Time delay(second)',...
    'y_name','Broadness factor');
%% Quality control
rmserr = zeros(wfpermod,tolmod);
indx = 1:nskip:size(rc_model,1); % index of fitted wf per wd
wdshort = zeros(nwf,length(indx)); % for fitted curve
fgauss = @(x,abc) abc(1).*exp(-((x-abc(2))./abc(3)).^2);
for iter = 1:tolmod
    [~,wd] = cqi.reflectivity(rc_model(:,iter),0,100);
    % get fitted truncated wavefield dictionary and cal rms error
    for iter2 = 1:length(indx)
        wdshort(:,iter2) = fgauss(coltick,fitpar_model{iter}(iter2,:));
        rmserr(iter2,iter) = rms(wd(:,indx(iter2))-wdshort(:,iter2));
    end
end
%% Next use linear fitting to fit data and get the distribution of k
b_fit = fitpar_matrix(end-100*wfpermod:end,2);
c_fit = fitpar_matrix(end-100*wfpermod:end,3);
%% Give a range of c, what is the pdf for b and what is the hitogam for b




