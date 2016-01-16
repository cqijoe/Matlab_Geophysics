function [ s ] = cqnotify( iter,total, dpercent )
% Output a string with percentage norification in a loop
% For example, in a loop with 100 iterations, iteration is 20 and
% dpercent is 20% then the cqnotify will out put a string:
%   ' 20% finished'
% If dpercent is 1% and there are 100 total iterations then after each
% iteration cqnotify will output information

n_notify = round(total*dpercent/100);
if iter==total
    s = sprintf('%3.0d%% finished',100);
else
    if mod(iter,n_notify)==0
        s = sprintf('%3.0d%% finished',floor(100*iter/total));
    else
        s = '';
    end
end

end

