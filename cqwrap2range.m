function [ d ] = cqwrap2range( d, min_range, max_range )
% Analogy to wrapTo2Pi, this function wrap data to a given range
% assuming data is periodic
% 
% input
% -----
% d = vector (row or column) to be wrapped
% min_range = minimum range 
% max_range = maximum range. d will be wrapped to [min_range, max_range]
%
% output
% ------
% d = after wrapping

% map date to range for [0,2*pi]
d_rad = (d - min_range)/(max_range - min_range)*2*pi;
% ^-- so when d = min_range, d_rad = 0; when d = max_range, d_rad = 2*pi
d_rad = wrapTo2Pi(d_rad);

% inverse d_rad to d
d = d_rad * (max_range - min_range) / 2 / pi + min_range;


end

