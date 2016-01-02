%% I am going to use Aratna Well to Stress the Problem
cd('C:\Users\Qi\Documents\MATLAB\cqi\cqsVelAnylasisAndTotalThickness\');
%% Load Aratna and pick coal and shale to make synthetic Aratna well
wlog = read_las_file('Aratna_vp.las');
vp = l_gc(wlog,'vp');
rhob = l_gc(wlog,'rhob');
figure; plot(vp);
set(gcf,'pos',[ 101         558        1746         420]);
%%
xflag = find(rhob<=1.2);
xflag = xflag(:);
%% Time domain model
% Modify original well log curves
vp = l_gc(wlog,'vp');
rhob = l_gc(wlog,'rhob');
vp(xflag==1) = 7700;
rhob(xflag==1) = 1.3;

% depth to time conversion
dt = 0.001;
wlog = cqi.l_depth2time(wlog,[vp,rhob],dt,vp,{'vpt','rhobt'});

% get coal beds part only
rangec = [1900:2400];
vpc = wlog.vpt(rangec);
rhobc = wlog.rhobt(rangec);

% get imp and rc and pm
imp = vpc.*rhobc;
rc = cqi.imp2rc(imp);
pm = cqi.pm(rc,0);
[~,wd] = cqi.reflectivity(rc,0,150);

% convolve rc and pm with a broadband wavelet with DC value removed
[~,rcwl] = ...
    cqi.genwav('ormsby',dt,[1,5,10,20],0,rc,'dc_co',true,...
    'norm','max','nsample',1001);
[wl,pmwl] = ...
    cqi.genwav('ormsby',dt,[1,5,10,20],0,pm,'dc_co',true,...
    'norm','max','nsample',1001);
%%
opt = glis.opt;
opt.MaxIter = 20;
opt.InitModel = mean(imp)*ones(length(imp),1);
opt.MuRange = [1e-20,1e-7];
opt.Source = 1;
opt.TrueModel = imp;
result = glis.dogli(pm, opt);
%% Test for recursive inversion
imp_ri = gli.ri(pm,[1,imp(1)]);


