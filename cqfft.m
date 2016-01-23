function [ y, f ] = cqfft( x, dt, t0, n, dim )
% fft transform with time axis. t need to be equally sampled
%
% input
% -----
% x = input data
% dt = time sampling rate
% t0 = the arrival time of the event (if t0 is a vector, it broadcast for matrix)
% n = number of samples you want to use in fft
% dim = along which dimension do you want to do fft
%
% output
% ------
% y = output data (0 ~ fnyq)
% f = frequency axis (0 ~ fnyq)
% 
% cqfft uses fft as the function to calculate fourier transform.
% then it uses time-shift theorem to modify phase spectrum

if ~exist('dim','var') || isempty(dim)
    if isrow(x)
        dim = 2;
    elseif iscolumn(x)
        dim = 1;
    else
        dim = 1;
    end
end
if ~exist('n','var') || isempty(n)
    if isrow(x) || iscolumn(x)
        n = length(x);
    else
        n = size(x,dim);
    end
end


y = fft(x, n, dim);

fnyq = 1/2/dt;
f = linspace(0,fnyq,n/2+1);
if dim == 1
     f0 = f(:);
     f0 = repmat(f0,1,size(x,2));
     y = y(1:n/2+1,:);
elseif dim == 2
    f0 = repmat(f,size(x,1),1);
    y = y(:,1:n/2+1);
end

if length(t0) == 1
    y = y .* exp(1i * 2 * pi * f0 * t0 );
else
    if dim == 1
        for k = 1:size(x,2);
            y(:,k) = y(:,k) .* exp(1i * 2 * pi * f0(:,k) * t0(k));
        end
    elseif dim == 2
        for k = 1:size(x,1);
            y(k,:) = y(k,:) .* exp(1i * 2 * pi * f0(k,:) * t0(k));
        end
    end
end






end

