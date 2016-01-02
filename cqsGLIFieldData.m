%% The script uses GLI from low freq to high freq to invert field data.
cd('C:\Users\Qi\Documents\MATLAB\cqi\cqsGLIFieldData');
%% 1. Data loading
sgy = read_segy_file('0093_pstm_CD_0000_1000_stk_gab_sigtr_IL680.sgy');
sgy = s_resample(sgy,1);
traces = sgy.traces;
trace = traces(:,407);

% load las file for well tie to determine source wavelet and scaling
% factor
wlog = read_las_file('Aratna_vp.las');
wlog = cqi.l_depth2time(wlog,{'vp','rhob'},0.001,'vp');
imp_time = wlog.vp_time.*wlog.rhob_time;
rc = cqi.imp2rc(imp_time);
pm = cqi.pm(rc,0);
par = cqi_quick_tie;
par.welltrnum = 107;
par.wlog = wlog;
trace_around_well = traces(:,300:440);
cqi_quick_tie(trace_around_well,rc,par);
% After I have run well tie, I found the best tying parameter for
% source wavelet is:
% ormsby filter | phase angle = 130 | timing = -121 sample points in
% 1ms sampling rate
%% Then find the scaling factor
% 1. generate seismogram in a window and synthetic data in a window
% 2. using the best tying parameters, calculate energy in the both
% traces
% 3. search for the best scaling factor with minimum residual energy
imp_coal = imp_time(1700+121:2400+121);
sori = trace(1700:2400);
s = sori;
d = pm(121:end);
d = d(1700:2400);
rc_coal = rc(121:end); rc_coal = rc_coal(1700:2400);
d_impulse = d;
[wl,dori] = cqi.genwav('ormsby',0.001,[2,10,60,80],130,d,'dc_co',true,...
    'norm','max','nsample',2001);

% both normalize to 1 according to the maximum peak
% Filter seismogram into very low freuqnecy
[~,s] = cqi.genwav('ormsby',0.001,[1,5,10,15],0,s,'dc_co',true,...
    'norm','max','nsample',2001);
[~,dori] = cqi.genwav('ormsby',0.001,[1,5,10,15],0,dori,'dc_co',true,...
    'norm','max','nsample',2001);
[~,wl] = cqi.genwav('ormsby',0.001,[1,5,10,15],0,wl,'dc_co',true,...
    'norm','max','nsample',2001);
s = s./max(s);
dmod = dori./max(dori);
wl_scale = wl./max(dori);
figure; plot(s); hold on; plot(dmod,'r');
%% give frequency
f = [1,5,10,15];
f3 = 10:5:80;
fmat = zeros(length(f3),4);
for iter = 1:length(f3)
    fmat(iter,:) = [1,5,0,0] + [0,0,f3(iter),f3(iter)+5];
end
result_cell = cell(length(f3),1);
%% 2. Fill option structure
opt= glis.opt;
opt.InitModel = ones(length(s)+1,1)*mean(imp_coal);
[wl,dori] = cqi.genwav('ormsby',0.001,[2,10,60,80],130,d,'dc_co',true,...
    'norm','max','nsample',2001);
for iter = 1:length(f3)
    opt.Figure = true;
    opt.MuNumber = 20;
    opt.MaxIter = 5;
    % filter seismic data and get estimated source wavelet
    [~,stemp] = cqi.genwav('ormsby',0.001,fmat(iter,:),0,sori,'dc_co',true,...
        'norm','max','nsample',2001);
    [~,wltemp] = cqi.genwav('ormsby',0.001,fmat(iter,:),0,wl,'dc_co',true,...
        'norm','max','nsample',2001);
    [~,dtemp] = cqi.genwav('ormsby',0.001,fmat(iter,:),0,dori,'dc_co',true,...
    'norm','max','nsample',2001);
    stemp = stemp./max(stemp);
    wl_scale = wltemp./max(dtemp);
    % -----------------------------------
   opt.Source = wl_scale;
   opt.TrueModel = imp_coal;
   result_cell{iter} = glis.dogli(stemp,opt);
   opt.InitModel = result_cell{iter}.m;
end
%% 2. Fill option structure
opt= glis.opt;
opt.InitModel = ones(length(s)+1,1)*mean(imp_coal);
[wl,dori] = cqi.genwav('ormsby',0.001,[2,10,60,80],130,d,'dc_co',true,...
    'norm','max','nsample',2001);
for iter = 1:length(f3)
    opt.Figure = true;
    opt.MuNumber = 20;
    opt.MaxIter = 40;
    % filter seismic data and get estimated source wavelet
    [~,stemp] = cqi.genwav('ormsby',0.001,fmat(iter,:),0,sori,'dc_co',true,...
        'norm','max','nsample',2001);
    [~,wltemp] = cqi.genwav('ormsby',0.001,fmat(iter,:),0,wl,'dc_co',true,...
        'norm','max','nsample',2001);
    [~,dtemp] = cqi.genwav('ormsby',0.001,fmat(iter,:),0,dori,'dc_co',true,...
    'norm','max','nsample',2001);
    stemp = stemp./max(stemp);
    wl_scale = wltemp./max(dtemp);
    % -----------------------------------
   opt.Source = wl_scale;
   opt.TrueModel = imp_coal;
   result_cell{iter} = glis.dogli(stemp,opt);
   opt.InitModel = result_cell{iter}.m;
end
% 2
opt= glis.opt;
opt.InitModel = ones(length(s)+1,1)*mean(imp_coal);
[wl,dori] = cqi.genwav('ormsby',0.001,[2,10,60,80],130,d,'dc_co',true,...
    'norm','max','nsample',2001);
for iter = 1:length(f3)
    opt.Figure = true;
    opt.MuNumber = 20;
    opt.MaxIter = 40;
    % filter seismic data and get estimated source wavelet
    [~,stemp] = cqi.genwav('ormsby',0.001,fmat(iter,:),0,sori,'dc_co',true,...
        'norm','max','nsample',2001);
    [~,wltemp] = cqi.genwav('ormsby',0.001,fmat(iter,:),0,wl,'dc_co',true,...
        'norm','max','nsample',2001);
    [~,dtemp] = cqi.genwav('ormsby',0.001,fmat(iter,:),0,dori,'dc_co',true,...
    'norm','max','nsample',2001);
    stemp = stemp./max(stemp);
    wl_scale = wltemp./max(dtemp);
    % -----------------------------------
   opt.Source = wl_scale;
   opt.TrueModel = imp_coal;
   result_cell2{iter} = glis.dogli(dtemp,opt);
   opt.InitModel = result_cell2{iter}.m;
end
save('matlab.mat');
%%
opt= glis.opt;
opt.InitModel = ones(length(s)+1,1)*mean(imp_coal);
for iter = 1:1
    opt.Figure = true;
    opt.MuNumber = 20;
    opt.MaxIter = 5;
    % filter seismic data and get estimated source wavelet
    stemp = sori;
    wltemp = wl;
    stemp = stemp./max(stemp);
    wl_scale = wltemp./max(dori);
    % -----------------------------------
   opt.Source = wl_scale;
   opt.TrueModel = imp_coal;
   result_cell{iter} = glis.dogli(stemp,opt);
   opt.InitModel = result_cell{iter}.m;
end
%% TV-Compression
[~,wd_coal] = cqi.reflectivity(rc_coal,0,100);
% Apply time-varying compression by extracting wd
[~,td] = max(wd_coal); td = td(:); td = [0;td];
td = td - 1;
tori = 1:length(result_cell{end}.m); tori = tori(:);
% do time compression on the inverted result
imp_inv_tvcmpress = interp1(tori,result_cell2{end}.m,tori+td);
%% Create freq filtered impedance profile
imp_wls = [];
demean_imp = imp_coal - smooth(imp_coal,100);
for iter = 1:length(f3)
    [~,temp] = cqi.genwav('ormsby',0.001,fmat(iter,:),0,demean_imp,'dc_co',true,...
    'norm','max','nsample',2001);
    imp_wls(:,iter) = cqi.norm2range1(temp);
end
%% Convert cell data to matrix data
inv_imp_mat = []; inv_imp_mat2 = [];
for iter = 1:length(f3)
    inv_imp_mat(:,iter) = cqi.norm2range1(result_cell{iter}.m);
    inv_imp_mat2(:,iter) = cqi.norm2range1(result_cell2{iter}.m);
end
%% Display
figure; cqi_plotmatrix(inv_imp_mat-0.5,'fillc','none');
hold on
cqi_plotmatrix(imp_wls-0.5,'fillc','none','linecolor','r','linewidth',2);
%% TV_Compression for each inversion result
inv_imp_mat_compress = [];
inv_imp_mat_compress2 = [];
for iter = 1:length(f3)
    inv_imp_mat_compress(:,iter) = ...
        interp1(tori,inv_imp_mat(:,iter),tori+td,'linear');
    inv_imp_mat_compress2(:,iter) = ...
        interp1(tori,inv_imp_mat2(:,iter),tori+td,'linear');
end
%% Display: tv compressed imp from seismogram
figure; cqi_plotmatrix(inv_imp_mat_compress2-0.5,'fillc','none');
hold on
%%
cqi_plotmatrix(imp_wls-0.5,'fillc','none','linecolor','r','linewidth',2);
%% 
[~,imp_coal_wl] = cqi.genwav('ormsby',0.001,[2,10,60,80],0,demean_imp,'dc_co',true,...
    'norm','max','nsample',2001);
imp_coal_wl = cqi.norm2range1(imp_coal_wl);
%% Below is the part for Surat , the previous part is stored in matlab.mat file
% in the folder
%+++++++++++++++++++++++++++++++++++++++++++++++
%+++++++++++++++++++++++++++++++++++++++++++++++
ccc
cd('C:\Users\Qi\Documents\MATLAB\cqi\cqsGLIFieldData');
% sgy = read_segy_file('gums1_cdplbls_1351.sgy');
sgy = read_segy_file('gums_mid_cdplbls_1351.sgy');
sgy = s_resample(sgy,1);
wlog = read_las_file('GUM_1.las'); 
wlog = l_unit_conversion(wlog,{'m','ft'});
%%
wlog = cqi.l_depth2time(wlog,{'vp','rhob'},0.001,'vp');
imp_ori = wlog.vp_time.*wlog.rhob_time;
rc_ori = cqi.imp2rc(imp);
cdp_gums = 333778;
% >>>>>>>>>>>>> search for trace number <<<<<<<<<<<<<<<<<<<<<<<
ntrace = find(sgy.headers(4,:) == cdp_gums);
if isempty(ntrace)
    disp('trace not found!')
else
    disp(['trace numer = ',num2str(ntrace)])
end
% >>>>>>>>>>>>> search for trace number <<<<<<<<<<<<<<<<<<<<<<<
% the trace number in the loaded sgy file is 244
trace = sgy.traces(:,ntrace);
traces = sgy.traces(1:2000,200:300);
par = cqi_quick_tie;
par.welltrnum = ntrace - 200;
%% Do Well Tie
cqi_quick_tie(traces,rc_ori,par);
%% Shift = 214, Phase = -100, freq = [1,5,23,55]
rc_tie = [zeros(214,1);rc_ori];
imp_tie = [imp(1)*ones(214,1);imp];
% give frequency
f = [1,5,10,15];
f3 = [10,15,23];
f4 = [20,35,55];
fmat = zeros(length(f3),4);
for iter = 1:length(f3)
    fmat(iter,:) = [1,5,0,0] + [0,0,f3(iter),f4(iter)];
end
result_cell = cell(length(f3),1);
s = trace(1:1500);
imp_coal = [imp_tie;ones(1500-length(imp_tie),1)*imp_tie(end)];
pm = cqi.pm(rc); 
%% 2. Fill option structure
opt= glis.opt;
opt.InitModel = ones(length(s)+1,1)*mean(imp);
[wl,dori] = cqi.genwav('ormsby',0.001,[1,5,23,55],-131,pm,'dc_co',true,...
    'norm','max','nsample',2001);
%%
for iter = 1:length(f3)
    opt.Figure = true;
    opt.MuNumber = 20;
    opt.MaxIter = 40;
    % filter seismic data and get estimated source wavelet
    [~,stemp] = cqi.genwav('ormsby',0.001,fmat(iter,:),0,s,'dc_co',true,...
        'norm','max','nsample',2001);
    [~,wltemp] = cqi.genwav('ormsby',0.001,fmat(iter,:),0,wl,'dc_co',true,...
        'norm','max','nsample',2001);
    [~,dtemp] = cqi.genwav('ormsby',0.001,fmat(iter,:),0,dori,'dc_co',true,...
    'norm','max','nsample',2001);
    stemp = stemp./max(stemp);
    wl_scale = wltemp./max(dtemp);
    % -----------------------------------
   opt.Source = wl_scale;
   opt.TrueModel = imp_coal;
   result_cell{iter} = glis.dogli(stemp,opt);
   opt.InitModel = result_cell{iter}.m;
end
save('matlab2.mat');
% system('shutdown -s');
%% Convert cell structure to matrix
inv_imp_mat = [];
for iter = 1:length(f3)
    inv_imp_mat(:,iter) = cqi.norm2range1(result_cell{iter}.m);
end
%% Do time-varying compression
[~,wd] = cqi.reflectivity(rc,0,100);
[~,td] = max(wd); td = td(:); td = [0;td];
td = td - 1;
tori = 1:length(result_cell{end}.m); tori = tori(:);
inv_imp_mat_compress = [];
for iter = 1:length(f3)
    inv_imp_mat_compress(:,iter) = ...
        interp1(tori,inv_imp_mat(:,iter),tori+td,'linear');
end
%% Create freq filtered impedance profile
imp_wls = [];
demean_imp = imp_coal - smooth(imp_coal,100);
for iter = 1:length(f3)
    [~,temp] = cqi.genwav('ormsby',0.001,fmat(iter,:),0,demean_imp,'dc_co',true,...
    'norm','max','nsample',2001);
    imp_wls(:,iter) = cqi.norm2range1(temp);
end
%% Plot
cqi_plotmatrix(inv_imp_mat,'fillco','none');
hold on
cqi_plotmatrix(imp_wls,'fillco','none','linecolor','r','linew',2);
xlim([1,4])
%% Plot Compression
cqi_plotmatrix(inv_imp_mat_compress,'fillco','none');
hold on
cqi_plotmatrix(imp_wls,'fillco','none','linecolor','r','linew',2);
xlim([1,4])
%% Generate Filter Dictionary
fitermat = [];
for iter = 1:size(fmat,1)
    filtermat(:,iter) = cqi.genwav('ormsby',0.001,fmat(iter,:),0,'dc_co',true,...
    'norm','max','nsample',2001);
end
%% Generate amplitude spectrum for each filter
amp_filter = [];
for iter = 1:length(f3)
    [X,F] = cqi_fft(filtermat(:,iter),0.001);
    amp_filter(:,iter) = abs(X);
end
%% Plot
figure;
for iter = 1:length(f3)
    subplot(length(f3)+1,1,iter);
    plot(imp_wls(:,iter),'r'); hold on;
    plot(inv_imp_mat(:,iter),'k')
%     plot(inv_imp_mat_compress(:,iter),'k')
    set(gca,'ytick',[]);
    set(gca,'xticklabel','');
end
subplot(length(f3)+1,1,length(f3)+1);
plot(cqi.norm2range1(imp_coal),'b');
set(gca,'ytick',[]);
linkaxes(findobj('type','axes'));
% xlim([1,1500])
xlim([1,700])
% set(gcf,'pos',[358         398        1274         580]);
set(gcf,'pos',[358          71        1428         907]);