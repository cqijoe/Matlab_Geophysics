function [ p ] = cqd2p( d,dt,dx )
% d is  degree and the function will find the p which is slope
% according to dt and dx.

d = d/180*pi;
p = dt/dx*tan(d);


end

