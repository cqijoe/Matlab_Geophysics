function [ dshift ] = cqshiftcor( d0 ,d, win )
% Shift traces with correlation
% 
% input
% -----
% d0 = data to be correlated
% d = data to be shifted and will be correlated with d0
% win = [top,bottom] range to be correlated
% 
% output
% ------
% dshift = shifted d

dshift = d;
if nargin < 3
    win = [1,size(d0,1)];
end
for k = 1:size(d0,2)
    [cor,lag] = xcorr(d0(win(1):win(2),k),d(win(1):win(2),k));
    [~,nmax] = max(cor);
    nshift = lag(nmax); % negative to shift up and positive to shift down
    if nshift < 0
        dshift(:,k) = [dshift(-nshift+1:end,k);zeros(-nshift,1)];
    elseif nshift > 0
        dshift(:,k) = [zeros(nshift,1);dshift(1:end-nshift,k)];
    end
end


end

