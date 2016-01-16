function [ PF ] = cqphasefilter_f( p, f,...
                   f_start, f_end, t_sweep,...
                   v_boat)
% frequency domain phase filter
p = reshape(p,1,length(p));
doppler = -p * v_boat;
amp = zeros(length(f),length(p));
phase = amp;


end

