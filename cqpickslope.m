function [ result ] = cqpickslope( ax )
% The function starts a while loop for continuing picking slopes
% press q to stop, left click to select two points for slope
% calculation

if ~exist('ax','var')
    ax = gca;
end
button = 1;
disp('press q to stop left click to pick...');
disp('please pick continued two points for slope calculation');
result = [];
n = 1;
axes(ax);
h_line = line(...
    'color','r',...
    'marker','o',...
    'markersize',10,...
    'linestyle','--');
h_text = text(...
    'fontsize',20,...
    'color','b');
h_line.XData = []; h_line.YData = [];
while button~=double('q')
    % continue picking untill user press q
    for k = 1:2
        [x(k),y(k),button] = ginput(1);
        if button==double('q')
            try
                delete(h_line);
                delete(h_text);
            catch
            end
            return
        else
            
            if k==1
                h_text.String = [];
                h_line.XData = x(k);
                h_line.YData = y(k);
            else
                h_line.XData(k) = x(k);
                h_line.YData(k) = y(k);
            end
        end
    end
    if diff(x)~=0
        result(n) = diff(y)/diff(x);
    else
        result(n) = nan;
    end
    % give text
    h_text.Position = [mean(x),mean(y)];
    h_text.String = num2str(result(n));
    n = n + 1;
end
end


