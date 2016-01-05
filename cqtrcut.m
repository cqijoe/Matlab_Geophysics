function [ d ] = cqtrcut( d )
% Truncate data interactively using polygons
%
% input
% -----
% d = input data
% 
% output
% ------
% d = input data after truncation

figure;
imagesc2(d);

input('Press Enter to Draw Polygons...');
s = 'n';
xq = 1:size(d,2);
yq = 1:size(d,1);
[xq, yq] = meshgrid(xq, yq);
while ~strcmp(s,'y')
    h = impoly;
    pos = wait(h);
    xv = pos(:,1);
    yv = pos(:,2);
    [in] = inpolygon(xq,yq,xv,yv);
    d = d.*(~in);
    s = input('Press any key to continue, press y to quit -> ','s');
end

end

