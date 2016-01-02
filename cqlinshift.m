function [ d2 ] = cqlinshift( d,dn )
% The function applies linear shift trace-by-trace
% d ... input data
% dn ... shifting step
% left-most trace will not be shift

d2 = zeros(size(d));
n = abs((0:size(d,2)-1) * dn);
d2(:,1) = d(:,1);

for iter = 2:size(d,2)
    if dn >= 0
        d2(n(iter)+1:end,iter) = d(1:end-n(iter),iter);
    else
        d2(1:end-n(iter),iter) = d(n(iter)+1:end,iter);
    end
end









end

