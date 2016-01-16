function [ dfk, f, k ] = cqfk( d, dt, dx )
% The function do the 2d fft on data, and output f-k data where f is
% from 0 ~ fnyq and k is from -knyq ~ knyq
%
% input
% d ... input matrix and each column correspondes to one x
% dt ... time sampling rate
% dx ... spatial sampling rate
%
% output
% dfk ... output matrixwith each column corresponds to one wavenumber
% f ... frequency sampling vector
% k ... wave number sampling vector
% -------------------------
%
% note: wavenumber is defined as 1/L where L is the wavelength. The
%       number of wavelength per unit distance

% ----------------------------

nf = 2^nextpow2(size(d,1));
nk = 2^nextpow2(size(d,2));
dfk = fft2(d,nf,nk);

dfk = dfk(1:nf/2+1,:);
dfk = [dfk(:,nk/2+2:end),dfk(:,1:nk/2+1)];

f = linspace(0,1/2/dt,nf/2+1);
dk = 1/dx/nk;
k = (-1/2/dx+dk):dk:1/2/dx;







end

