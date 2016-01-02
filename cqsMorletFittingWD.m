% Data loading
load('C:\Users\Qi\Documents\MATLAB\cqi\cqsFreqAndDelay\wlog_mc.mat');
%% MF decomp for a certain model
nmod = 100;
rc_model = cqi.imp2rc(wlog_mc.impmodel_time);
rc_this = rc_model(:,nmod);
pm_this = cqi.pm(rc_this);
[~,wd] = cqi.reflectivity(rc_this,0,120);
opt = smp.opt;
%% MF Option
opt.fmax = 100;
opt.nmax = 1;
opt.rphi = 0;
%% Generate wd_wl
bandf = [1,5,80,100];
[wl,wd_wl] = cqi.genwav('ormsby',0.001,bandf,0,wd,'nsample',1001);
%% Do MF on wavefield dictionary
rcpos = 300;
result = cell(rcpos,1);
for n = 1:rcpos
    [ result{n} ] = smp.main( wd_wl(:,n), opt );
end
%% Unparse wavepar to a large matrix
wavepar = zeros(rcpos,5);
for n = 1:rcpos
    wavepar(n,:) = result{n}.wavepar;
end
%% Scatter plot u (1) versus omega (3)
scatter(wavepar(:,1),wavepar(:,3)/2/pi);
%% do decomp on seismogram
[wl,pmwl] = cqi.genwav('ormsby',0.001,bandf,0,pm_this,'nsample',1001);
optpm = opt;
optpm.nmax = 300;
result_pm = smp.main(pmwl,optpm);
%% get wavepar
wavepar_pm = result_pm.wavepar;
scatter(wavepar_pm(:,1),wavepar_pm(:,3)/2/pi)
hold on
scatter(wavepar(:,1)+0.001*(0:299)',wavepar(:,3)/2/pi,'r')
%% Quanlity Control
mwave = zeros(size(result{1}.wlmat,1),length(result));
for n = 1:length(result)
    mwave(:,n) = result{n}.wlmat;
end
%% 
cqi_plotmatrix(wd_wl,'skip',10);
%%
hold on
cqi_plotmatrix(mwave,'skip',10,'fill_co','none','linecolor','r','linew',2);
%% filtering the seismogram with morlet wavelet 
freq = (2:0.5:150)';
omega = 2*pi*freq; % freq is from 1 to 20 Hz with 1Hz increment
t = (-500:500)'*0.001; % time vector
nzero = find(t==0);
sigma = 0.5; % from the model 100 fitting result
u = 0;
morlet_wd = zeros(length(t),length(omega));
pm = [pm_this;zeros(300,1)];
pm_filtered = zeros(length(pm),length(freq));

for n = 1:length(omega)
    morlet_wd(:,n) = smp.morlet_wavelet([u,sigma,omega(n),0],t);
    pm_filtered(:,n) = ...
        cqi.convsig(pm,real(morlet_wd(:,n)),nzero);
end
%%
time = 0:0.001:(size(pm_filtered,1)-1)*0.001;
imagesc2(freq,time,pm_filtered); xlabel('Frequency(HZ)'); ylabel('Time(s)')
colormap jet;
result = cqpick('window_max',pm_filtered,freq,time);
%%
hold on
plot(freq,result.peak_z,'k','linewidth',2);
figure;
plot(freq,result.peak); xlabel('Frequency(HZ)'),ylabel('Amplitude at Peak');
figure
plot(wavepar(:,3)/pi/2);







