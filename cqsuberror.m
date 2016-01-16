function [ a_b ] = cqsuberror( a,b )
% get subtraction error between two matrixes using maximum value in
% each matrix
% a and b should be in the same size

a_b = a./max(max(a)) - b./max(max(b));

end

