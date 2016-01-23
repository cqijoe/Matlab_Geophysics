function [data]= cqmirror(data,angRef,pointRef)
% cqmirror:
% Mirror data about a reference line passing through a reference point
% at a reference angle, using reflection matrix
% data: [x y] vector format
% angRef: radian
% pointRef: [x y] format
% Example:
% data = [-1 2;-3 1;-2 5;-1 2];
% angRef = 45*pi/180;
% pointRef = [0 1];
%
% [dataMirrored]= transf_Reflect(data,angRef,pointRef);
% x = -4:4;
% y = x+1;
% figure
%
% plot(data(:,1),data(:,2),x,y,dataMirrored(:,1),dataMirrored(:,2))
% axis equal
% grid

% Prior to mirroing, make zero point of x-y coordinate about the reference point
[data]= transf_Translate(data,-pointRef);
data = data'; % transpose for matrix multiplication


matrixCoeff = [cos(angRef*2) sin(angRef*2); sin(angRef*2) -cos(angRef*2)]; % for reflection
data = matrixCoeff*data;

data = data'; % back to original vector format
[data]= transf_Translate(data,pointRef); % transform data about original x-y coordinate

function [data]= transf_Translate(data,shift)
data(:,1) = data(:,1) + shift(1);
data(:,2) = data(:,2) + shift(2);