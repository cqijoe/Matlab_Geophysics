function [ coeffq ] = cqcomplex_interp( coeff, x, qval, term )
% complex number interpolation using Fourier series
% method referenced to Wade's UH Thesis, 1989
%
% input
% -----
% coeff = Complex(or real) Numbers with length N
% x = vector of x axis of coeff, x must increasing and equally sampled
%     if not, program will sort your x and test if it is equally
%     sampled
% qval = query values of coeff
% term = 6 = how many terms you want to use to interpolate
% 
% output
% ------
% coeffq = complex value at the query frequency, it is a column vector


if ~exist('term','var')||isempty(term)
    term = 6;
end

% check input
dx = x(2) - x(1);
if length(x) ~= length(coeff)
    error('Invalid x axis! Must equal to coeff length!');
end
[x,I] = sort(x);
coeff = coeff(I);

if any( (diff(x) - dx) > 1e-5)
    error('Invalid x axis! Must be equally sampled!');
end

% make sure the query x is within x limit
if max(qval) > max(x) || min(qval) < min(x)
    error('Invalid query position! Must within range of x!');
end

coeff = coeff(:);
qval = qval(:);



% find nearest term number of data to interpolate and do interpolation
% point-by-point
coeffq = zeros(length(qval),1);
for k = 1:length(qval)
    % sort the difference and select the minimum term
    float_index_distance = (x - qval(k))/dx;
    [~,I] = sort(abs(float_index_distance));
    I = I(1:term);
    float_index_distance = float_index_distance(I);
    float_index_distance = float_index_distance(:);
    coeffq(k) = ...
        sum( coeff(I) .* exp(1i*float_index_distance*pi/2) .* ...
        sinc(float_index_distance/4).^2 .* ...
        sinc(float_index_distance) );
    % Matlab sinc function: sin(pi*x)/(pi*x)
end



end

