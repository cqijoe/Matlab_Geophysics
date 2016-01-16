function [ d ] = cqp2d( p, dt, dx )
% given p I return d in degree
d = atan(dx/dt*p);
d = d/pi*180;


end

