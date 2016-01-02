function [ d2 ] = cqgeogain( d, dt, t0, a )
% The function apply t^a geometric spreading gain to each column of
% input data d

if nargin < 3
    t0 = 0;
end
if nargin < 4
    a = 2;
end

t = 0:dt:(size(d,1)-1)*dt + t0;
t = t';
t = t.^a;
d2 = d.*repmat(t,1,size(d,2));


end

