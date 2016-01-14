function [ xt, xt_bg, x ] = cqfktpinv(tp, dt, p, nx, dx, fk_bg, term)
% Inverse Tau-P transform using fk domain method
% 
% input
% -----
% tp = tau - p data matrix
% dt = sampling rate in seconds
% p = vector of the p in tau - p matrix
% nx = number of traces do you want
% dx = x sampling rate in meters
% xleft = what is the offset of the left-most trace in your
%         reconstructed trace?
% fk_bg = (optional) background fk data generated by cqfktp method
%
% output
% ------
% xt = x - t data matrix without background fk filling
% xt_bg = x - t data matrix with background fk filling

% Check input
% ------------
if length(p)~=size(tp,2)
    error('Invalid p vector. Length error!')
end


% Transform tp data to fp data
% ---------------------------
nt0 = size(tp,1);
nt = 2^nextpow2(nt0);
fp = fft(tp,nt,1);
fp = fp(1:nt/2+1,:);
fnyq = 1/2/dt;
f = linspace(0,fnyq,nt/2+1);
f_index = 1:(nt/2+1); 

% Output initialization
% ---------------------
fk = zeros(length(f),nx);
knyq = 1/2/dx;
k = linspace(0,2*knyq,nx+1);
k(end) = [];
f_reserve = cell(nx,1); % each cell stores indexes to be filled by background

% Fill the fk data trace-by-trace via k
% -------------------------------------
for trace = 1:length(k)
    this_f_index = f_index;
    this_p = [0,-k(trace)./f(2:end)]; % zero in f can not be included
    this_f = f;
    n_rm = [1,find(this_p>max(p) | this_p<min(p))];
    f_reserve{trace} = n_rm;
    this_f(n_rm) = []; % remove outsiders
    this_f_index(n_rm) = [];
    this_p(n_rm) = [];
    % loop into every f that's within p range to get its fp value
    for m = 1:length(this_f)
        fk(this_f_index(m),trace) = ...
            cqcomplex_interp(fp(this_f_index(m),:),p,this_p(m),term);
    end
end

% Inverse fk before filling
% --------------------------
fk0 = [real(fk);flipud(real(fk(2:end-1,:)))] + ...
    1i * [imag(fk);-flipud(imag(fk(2:end-1,:)))];
xt = real(ifft2(fk0));
xt = xt(1:nt0,:);
% ifft shift
xt = [xt(:,nx/2+2:end),xt(:,1:nx/2+1)];
x = (-(nx/2-1)*dx):dx:(nx/2)*dx;


% If fk_bg is provided, fill the hole in fk domain by fk_bg
% ----------------------------------------------------------
if exist('fk_bg','var')
    if size(fk_bg,2)==nx && size(fk_bg,1)==(nt/2+1)
        for trace = 1:length(k)
            fk(f_reserve{trace},trace) = ...
                fk_bg(f_reserve{trace},trace);
        end
    else
        warning('Background fk data can''t be used due to size issue!')
    end
end

% Inverse fk after filling
% --------------------------
fk1 = [real(fk);flipud(real(fk(2:end-1,:)))] + ...
    1i * [imag(fk);-flipud(imag(fk(2:end-1,:)))];
xt_bg = real(ifft2(fk1));
xt_bg = xt_bg(1:nt0,:);
% ifft shift
xt_bg = [xt_bg(:,nx/2+2:end),xt_bg(:,1:nx/2+1)];

end
