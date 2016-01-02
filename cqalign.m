function [ x2 ] = cqalign( x, dt, t_first, t0 )
% align a matrix in columns based on time vector and corresponding
% alignment time
% 
% input
% x ... 2D matrix in columns
% dt ... sample rate for x in vertical
% t_first ... first arrival time for each column in x
% t0 ... corresponding alignment time. all traces in x would be
%       shifted to t0 based on t_first
% output
% x2 ... x after alignment

if length(t_first)~=size(x,2)
    error('Invalid first arrival time t_first');
end
x2 = x;
nshift = round((t_first - t0)/dt);

trace_up = find(nshift>0);
trace_down = find(nshift<0);

for k = 1:length(trace_up)
    shift = nshift(trace_up(k));
    x2(1:end-shift,trace_up(k)) = x(shift+1:end,trace_up(k));
end

for k = 1:length(trace_down)
    shift = -nshift(trace_down(k));
    x2(:,trace_down(k)) = [zeros(shift,1);x(1:end-shift,trace_down(k))];
end


end

