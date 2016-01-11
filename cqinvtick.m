function cqinvtick(invtick,rotate,xory,ax)
% inverse tick label to be its inverse, i.e. 0.001 will be as 1000 in
% displaying
%
% input
% -----
% invtick = vector of ticks (after inverse) to be labeled
% xory = 'x' = which axis you want to inverse the tick label
%      = 'y'
% rotate = 45 = degree of rotation
% ax = gca = which axes you want to modify

if nargin < 2
    rotate = 45;
end
if nargin < 3
    xory = 'x';
end
if nargin < 4
    ax = gca;
end
% ----------------
invtick_no_0 = invtick(invtick~=0);
invinvtick = 1./invtick_no_0;
[invinvtick,I] = sort(invinvtick);
invtick_no_0 = invtick_no_0(I);
if strcmp(xory,'x')
    set(ax,'xticklabelmode','manual',...
        'xtick',invinvtick,...
        'xticklabel',invtick_no_0,...
        'xticklabelrotation',rotate);
end
if strcmp(xory,'y')
    set(ax,'yticklabelmode','manual',...
        'ytick',invinvtick,...
        'yticklabel',invtick_no_0,...
        'yticklabelrotation',rotate);
end



end