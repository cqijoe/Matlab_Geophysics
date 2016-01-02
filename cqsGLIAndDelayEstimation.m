%% Data Loading and Modification
cd('C:\Users\Qi\Documents\MATLAB\cqi\Data\Data_Cooper_Basin\');
wlog = read_las_file('Aratna_vp.las');
load('C:\Users\Qi\Documents\MATLAB\cqi\cqsVelAnylasisAndTotalThickness\AratnaCoalPicking.mat');
%% GLI on the aratna well
vp = l_gc(wlog,'vp');
rhob = l_gc(wlog,'rhob');
index = 1:length(vp); index = index(:);
cInd = find(index>7960&vp<10000&index<10000);
sInd = find(index>7960&vp>=10000&index<10000);
% Uncomment this afterwards
% vp(cInd) = 7700;
% rhob(cInd) = 1.35;
% count coal layer number
n = 1; ctop = []; cbot = []; cthick = [];
ctop(1) = cInd(1);
for m = 1:length(cInd)-1
    if (cInd(m+1)-cInd(m))~=1
        n = n + 1;
        ctop(n) = cInd(m+1);
        cbot(n-1) = cInd(m);
    end
end
cbot = [cbot,cInd(end)];
cthick = ctop - cbot + 1;

% vp(index>7960&vp>=10000) = 14000;
% rhob(index>7960&vp>=10000) = 2.5;

% depth to time conversion
dt = 0.001;
wlog = cqi.l_depth2time(wlog,[vp,rhob],dt,vp,{'vpt','rhobt'});

% get coal beds part only
rangec = (1900:length(wlog.vpt))';
vpc = wlog.vpt(rangec); 
rhobc = wlog.rhobt(rangec);

% Uncomment this afterwards
% vpc = [ones(100,1)*vpc(1);vpc];
% rhobc = [ones(100,1)*rhobc(1);rhobc];

% get imp and rc and pm
imp = vpc.*rhobc;
rc = cqi.imp2rc(imp);
pm = cqi.pm(rc,0);
[~,wd] = cqi.reflectivity(rc,0,150);

% convolve rc and pm with a broadband wavelet with DC value removed
[wl0,rcwl] = ...
    cqi.genwav('ormsby',dt,[10,20,30,40],0,rc,'dc_co',true,...
    'norm','max','nsample',2001);
[wl,pmwl] = ...
    cqi.genwav('ormsby',dt,[1,5,10,15],0,pm,'dc_co',true,...
    'norm','max','nsample',2001);
%% Plot
figure;wtva(rcwl); hold on
wtva(pmwl+1,1:length(pmwl),'k',1); plot(cqi.norm2range1(imp));
set(gcf,'pos',[339         419        1387         567]);
xlim([1,600]);
set(gca,'fontsize',20);
grid off
%% GLI Inversion
opt = glis.opt;
opt.MaxIter = 10;
opt.InitModel = imp(1)*ones(length(imp),1);
% opt.InitModel = smooth(imp,50);
opt.MuRange = [1e-20,1e-10];
opt.Source = wl;
opt.TrueModel = imp;
resultwl = glis.dogli(pmwl, opt);
%% Plot
% filter to low freq for imp
[wl,imp_low] = ...
    cqi.genwav('ormsby',dt,[1,5,10,15],0,imp-mean(imp),'dc_co',true,...
    'norm','max','nsample',2001);
imp_low = imp_low + mean(imp);
figure; plot(imp_low); hold on; plot(resultwl.m,'r'); plot(opt.InitModel,'k');
xlabel('Time(ms)'); ylabel('Impedance');
set(gcf,'pos',[508         563        1356         257]);
set(gca,'fontsize',15); xlim([1,600])
%% GLI on impulse
opt = glis.opt;
opt.MaxIter = 10;
opt.InitModel = imp(1)*ones(length(imp),1);
opt.MuRange = [1e-20,1e-10];
opt.TrueModel = imp;
result = glis.dogli(pm, opt);
%% Plot
figure; plot(imp); hold on; plot(result.m,'r'); plot(opt.InitModel,'k');
xlabel('Time(ms)'); ylabel('Impedance');
set(gcf,'pos',[508         563        1356         257]);
set(gca,'fontsize',15); xlim([1,600])
%%
% Apply time-varying compression by extracting wd
[~,td] = max(wd); td = td(:); td = [0;td];
td = td - 1;
tori = 1:length(resultwl.m); tori = tori(:);
% do time compression on the inverted result
imp_inv_tvcmpress = interp1(tori,resultwl.m,tori+td);

%% Filtering the imp to seismic bandwidth
[~,impwl] = cqi.genwav('ormsby',0.001,[1,5,10,15],0,imp-mean(imp),...
    'dc_co',true,'nsample',2001,'norm','max');
impwl = cqi.norm2range1(impwl);
%% Test: impedance matrix and weighted summing
f3 = [10:10:100];
f4 = [15:10:115];
f = zeros(length(f3),4);
wlmat = [];
pmmat = []; % matrix for pm filtering by corresponding wavelet
for n = 1:length(f3)
    % frequency for filtering
    f(n,:) = [1,5,f3(n),f4(n)];
    % calculate wavelet and convolve with impulse syn
    [wlpm,pmmat(:,n)] = cqi.genwav('ormsby',0.001,f(n,:),0,...
        pm,'dc_co',true,'nsample',2001,'norm','max');
    [wlmat(:,n),~] = cqi.genwav('ormsby',0.001,f(n,:),0,...
        pm,'dc_co',true,'nsample',2001,'norm','max');
end

invcell = cell(length(f3),1);
imp0 = imp(1)*ones(length(imp),1);
% do GLI
for n = 1:length(f3)
    opt = glis.opt;
    opt.MaxIter = 5;
    opt.InitModel = imp0;
    opt.MuRange = [1e-15,1e-10];
    opt.Source = wlmat(:,n);
    opt.TrueModel = imp;
    invcell{n} = glis.dogli(pmmat(:,n),opt);
    imp0 = invcell{n}.m;
end
beep on;
for n = 1:3
    beep;
end
disp('All frequency inverted!!!!')
%% Display Inversion Result
% parse the resulted impedance to impmat
impmat = [];
impmat_norm = [];
for n = 1:length(invcell)
    impmat(:,n) = invcell{n}.m;
    impmat_norm(:,n) = cqi.norm2range1(impmat(:,n));
end

cqi_plotmatrix(impmat_norm,'fillc','none','lineco','r','linew',2);
%% Calculate corresponding true model filtered by wavelets
imptrue = [];
imptrue_norm = [];
for n = 1:length(f3)
    [~,imptrue(:,n)] = cqi.genwav('ormsby',0.001,f(n,:),0,...
        imp-mean(imp),'dc_co',true,'nsample',2001,'norm','max');
    imptrue_norm(:,n) = cqi.norm2range1(imptrue(:,n));
end
hold on
cqi_plotmatrix(imptrue_norm,'fillc','none','linecolor','k',...
    'linewidth',2)
pos2 = [898   328   467   570];
set(gcf,'pos',pos2);
xlim([1,10])
%% Apply tv-compression based on delay curve from wavefield dictionary
impmat_tvcompress = [];
impmat_tvcompress_norm = [];
for n = 1:length(f3)
    impmat_tvcompress(:,n) = ...
        interp1(tori,impmat(:,n),tori+td,'linear',imp(end));
    impmat_tvcompress_norm(:,n) = ...
        cqi.norm2range1(impmat_tvcompress(:,n));
end
cqi_plotmatrix(imptrue_norm,'fillc','none','linew',2);
hold on
cqi_plotmatrix(impmat_tvcompress_norm,'fillc','none','linecolor','r',...
    'linewidth',2);
pos2 = [898   328   467   570];
set(gcf,'pos',pos2);
xlim([1,10])
%% MORE PLOT
figure; plot(imptrue_norm(:,end)); hold on;xlim([1,500])
plot(impmat_tvcompress_norm(:,end),'r')
set(gcf,'pos',pos);
%% Calculate correlation coefficient
corr_coef = [];
for n = 1:length(f3)
    temp = corrcoef(imptrue_norm(1:500,n),impmat_tvcompress_norm(1:500,n));
    corr_coef(n) = temp(1,2);
end
%% Load modeling result
load('wlog_mc.mat');
%% Or regenerate one with more sample logs
cqi.coalmodel; % output to wlog_mc
%% Calculate delay curve data base
% unparse some matrix
impmodel = wlog_mc.impmodel_time;
rcmodel = cqi.imp2rc(impmodel);

% using a for loop, calculate time delay curve and store them, only
% calculate wd not storing wd for 1000 models
tdmodel = [];
parfor n = 1:size(rcmodel,2)
    fprintf('n = %d\n',n);
    [~,wdmodel] = cqi.reflectivity(rcmodel(:,n),0,100);
    [~,temp] = max(wdmodel);
    temp = temp - 1;
    temp = temp(:);
    tdmodel(:,n) = temp;
end
%% build workspace data for curve fitting app
index = 1:size(tdmodel,1); index = index(:);
fitpar = cell(size(tdmodel,2),1);
parfor n = 1:length(fitpar)
    fprintf('n = %d\n',n);
    [~,bb] = max(tdmodel(:,n));
    fitpar{n} = ...
        cqsGLIAndDelayEstimation_createFit(index,tdmodel(:,n),[50,bb,30],false);
end
%% Unparse fitpar structure to extract relationship between d-a and b
b = [];
a = [];
c = []; % max delay curve
h = []; % thickness of the coal
for n = 1:length(fitpar)
    temp = coeffvalues(fitpar{n});
    b(n) = temp(2);
    a(n) = temp(1);
    c(n) = temp(3);
end
h = b - a;
b(h<100&h>500) = [];
c(h<100&h>500) = [];
a(h<100&h>500) = [];
h(h<100&h>500) = [];

hap  = c + h; % apparent thickness
k = c./(b-a);  % slope
% scatter plot
figure
scatter(hap,c)
ylim([0,100])
xlabel('Apparent Thickness(ms)')
ylabel('Maximum Time Delay(ms)')
set(gca,'fontsize',20); ylim([0,60])
set(gcf,'pos',[ 680         558        1089         420]);

figure
scatter(hap,k)
xlabel('Apparent Thickness(ms)')
ylabel('Slope of Delay Curve')
set(gca,'fontsize',20);ylim([0,0.3])
set(gcf,'pos',[ 680         558        1089         420]);
%% Get gaussian binomial distribution
xd = 200:450;
yd = linspace(0,0.2,100);
[ f_xy ] = cqi_bidistr( hap,k,xd,yd,false );
figure;scatter(hap,k,'k','k','filled');hold on;
contour(xd,yd,f_xy,5,'linew',2); set(gca,'fontsize',20);
axis xy;colormap jet; colorbar;xlabel('Apparent Thickness(ms)');
ylabel('Slope'); set(gcf,'pos',[781         242        1072         376]);

yd2 = linspace(0,50,100);
[ f_xy2 ] = cqi_bidistr( hap,c,xd,yd2,false );
figure;scatter(hap,c,'k','k','filled');hold on;
contour(xd,yd2,f_xy2,5,'linew',2); set(gca,'fontsize',20);
axis xy;colormap jet; colorbar;xlabel('Apparent Thickness(ms)');
ylabel('Max Delay(ms)'); set(gcf,'pos',[781         242        1072         376]);

figure
plot(yd2,f_xy2(:,80))
xlabel('Max Delay(ms)')
ylabel('PDF')
set(gcf,'pos',[680   779   560   199]);
set(gca,'fontsize',20,'xtick',[10:5:40]);
set(gca,'ytick',[])

figure
plot(yd,f_xy(:,80))
xlabel('Slope')
ylabel('PDF')
set(gcf,'pos',[680   779   560   199]);
set(gca,'fontsize',20);
set(gca,'ytick',[])

%%
ind_need = 80;
delayfunc = @(x,a,k,h) 0*(x<=a) + (k*x-k*a).*(x>a&x<(a+h))+(k*h)*(x>=(a+h));
index = 1:length(td);
esttd = delayfunc(index,50,0.1232,250);
esttd = [zeros(1,50-round(a(ind_need))),esttd];
figure; plot(esttd,'b'); hold on
plot(td,'r')
esttd = delayfunc(index,50,0.1455,250);
esttd = [zeros(1,50-round(a(ind_need))),esttd];
plot(esttd,'b--');
esttd = delayfunc(index,50,0.0929,250);
esttd = [zeros(1,50-round(a(ind_need))),esttd];
plot(esttd,'b--');
set(gca,'fontsize',17);
xlabel('Time(ms)');ylabel('Delay(ms)'); ylim([0,50])
%% 
esttd = delayfunc(index,50,0.1232,250);
esttd = [zeros(1,50-round(a(ind_need))),esttd];
% Apply time-varying compression by extracting wd
esttd = esttd(:);
esttd = esttd(1:length(imp));
tori = 1:length(resultwl.m); tori = tori(:);
% do time compression on the inverted result
imp_inv_tvcmpress_est = interp1(tori,resultwl.m,tori+esttd);
%% plot
figure;
plot(cqi.norm2range1(imp_low),'b'); hold on
plot(cqi.norm2range1(imp_inv_tvcmpress),'r');
plot(cqi.norm2range1(imp_inv_tvcmpress_est),'k');
xlim([1,500]); set(gcf,'pos',[680         585        1024         393]);
xlabel('Time(ms)');
ylabel('Normalized Scale');
set(gca,'fontsize',20);
