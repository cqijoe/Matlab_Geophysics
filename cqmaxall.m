function cqmaxall
% maximize all figures
fig = findobj('type','figure');
window_size = get(0,'screensize');
set(fig,'position',window_size);

end

