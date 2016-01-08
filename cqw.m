function [ data_condition,hax ] = ...
    cqw( data,x,y,...
    skip,lvl,...
    hax,lvar,xtickskip,...
    iftruex)
% [ data_condition,hax ] = cqw( data,x,y,skip,lvl,hax,lvar,xtickskip )
% plot wiggle form seismogram
% 
% input:
% data ... 2D matrix of data in columns
% x ... xlabel
% y ... ylabel
% skip ... skip traces to display (default = 0)
% lvl ... gain level (0 to not normalize data)
% hax ... axes to plot
% lvar ... line plotting variables to pass to line function
%          default lvar = {'color','k','linewidth',1}
% xtickskip ... skip of the xtick label
% iftruex ... <false> not use true x as horizontal axis label
% 
% output:
% hax ... return axis
% data_condition ... data is normalized and conditioned to the display parameter
% 
% In order to skip some parameter to use default, give is empty ([])
if ~exist('iftruex','var')||isempty(iftruex)
    iftruex = false;
end
if ~exist('xtickskip','var')||isempty(xtickskip)
    xtickskip = 10;
end
if ~exist('hax','var')||isempty(hax)
    hax = gca;
end
if ~exist('skip','var')||isempty(skip)
    skip = 0;
end
if ~exist('lvar','var')||isempty(lvar)
    lvar = {'color','k','linewidth',1};
end
if ~exist('lvl','var')||isempty(lvl)
    lvl = 1;
end
if ~exist('y','var')||isempty(y)
    y = 1:size(data,1);
end
if ~exist('x','var')||isempty(x)
    x = 1:size(data,2);
end
% check xlabel and ylabel
if length(x)~=size(data,2)
    error('Invalid x length!')
end
if length(y)~=size(data,1)
    error('Invalid y length!')
end
% get dx
 
skip = skip + 1;
if skip~=0
    if iftruex
        xstd = x(1:skip:end);
    else 
        xstd = 1:size(data,2);
        xstd = xstd(1:skip:end);
    end
    x = x(1:skip:end);
end
% normalize data
if lvl~=0
    maxdata = max(max(data));
    data = data/maxdata*lvl;
end
data_condition = data;
data = data(:,1:skip:end);
data = data + repmat(xstd,length(y),1);
% plot 
axes(hax);
y = y(:); y = repmat([y;nan],length(x),1);
data = [data;ones(1,size(data,2))*nan];
ndata = size(data,1)*size(data,2);
data = reshape(data,ndata,1);
line(data,y,lvar{:});
% post-processing
xlim([min(xstd),max(xstd)]);
ylim([min(y),max(y)]);
axis ij;
xtickskip = xtickskip + 1;
set(hax,'xaxislocation','top','xticklabelrotation',60,...
    'tickdir','out');
if ~iftruex
    set(hax,'xtick',xstd(1:xtickskip:end),'xticklabel',x(1:xtickskip:end));
end


end

