function [ result ] = cqpick( mode, par, x, z, ax )
% pick module for different types of picking using ginput
% mode .. {<'window_max'>,'well','seismic')
% par ... different mode require different parameter
%   window_max ... matrix data on the screen
%   well ... no need to provide anything, just []
%   seismic ... no need to provide anything, just []
% x ... horizontal axis
% z ... vertical axis
% ax ... axis to pick, by default it is current axis, namely gca()
% result ... different mode gives differnet result structure

if ~exist('ax','var')
    ax = gca;
end
axes(ax);
switch mode
    case 'window_max'
        result.wintop = [];
        result.winbot = [];
        result.peak_z = []; % store z of the peak
        result.peak = []; % store peak value
        result.peak_zind = []; % store z index of the peak
        % ----auxillary lines------------
        axes(ax);
        h_line(1) = line(...
            'Color','m',...
            'LineStyle','--',...
            'Marker','>',...
            'MarkerFaceColor','r',...
            'MarkerSize', 5); % window top
        h_line(2) = line(...
            'Color','k',...
            'LineStyle','--',...
            'Marker','<',...
            'MarkerFaceColor','y',...
            'MarkerSize', 5); % window bottom
        h_line(3) = line(...
            'Color','r',...
            'LineStyle','-'); % picked peak
        %% 1. pick window boundary
        % top window *********************************
        input('Pick window top (Left pick, Right finish)');
        button = 1;
        n = 1;
        while(button==1)
            [result.wintop(n,1),result.wintop(n,2),button] = ginput(1);
            if n == 1
                h_line(1).XData = result.wintop(n,1);
                h_line(1).YData = result.wintop(n,2);
            else
                h_line(1).XData = [h_line(1).XData, result.wintop(n,1)];
                h_line(1).YData = [h_line(1).YData, result.wintop(n,2)];
            end
            n = n + 1;
        end
        result.wintop = ...
            interp1(result.wintop(:,1),result.wintop(:,2),...
            x, 'linear', 'extrap');
        h_line(1).XData = x;
        h_line(1).YData = result.wintop;
        h_line(1).Marker = 'none';
        % bottom window ******************************
        input('Pick window bottom (Left pick, Right finish)');
        button = 1;
        n = 1;
        while(button==1)
            [result.winbot(n,1),result.winbot(n,2),button] = ginput(1);
            if n == 1
                h_line(2).XData = result.winbot(n,1);
                h_line(2).YData = result.winbot(n,2);
            else
                h_line(2).XData = [h_line(2).XData, result.winbot(n,1)];
                h_line(2).YData = [h_line(2).YData, result.winbot(n,2)];
            end
            n = n + 1;
        end
        result.winbot = ...
            interp1(result.winbot(:,1),result.winbot(:,2),...
            x, 'linear', 'extrap');
        h_line(2).XData = x;
        h_line(2).YData = result.winbot;
        h_line(2).Marker = 'none';
        %% 2. program auto pick the peak values within window
        for n = 1:length(x)
            top = min([result.winbot(n),result.wintop(n)]);
            bot = max([result.winbot(n),result.wintop(n)]);
            winindex = find(z<=bot&z>=top);
            [result.peak(n),ind] = max(par(winindex,n));
            result.peak_zind(n) = winindex(ind);
            result.peak_z(n) = z(result.peak_zind(n));
        end
        h_line(3).XData = x;
        h_line(3).YData = result.peak_z;
        %% 3. delete?
        ifdelete = input('Delete lines?(y/n)','s');
        if strcmp(ifdelete,'y')
            delete(h_line);
        end
    case 'well'
        % in this mode, program suppose you have a well log curve
        % plotted in the way that: x indicates position and z(y)
        % indicates logging value
        % result will be a cell array with structures in the following
        % manner: 
        %   n is the number of lithology
        %   result{n}.name ... lithology name
        %   result{n}.top ... lithology top positions
        %   result{n}.bottom ... lithology bottom position
        %                   note: top and bottom will have same length
        %   result{n}.xflag ... all zero copy of input x with the zones
        %                 with such lithologies being assigned value 1
        %   result{n}.thickness ... thickness in each block
        n = 1; % lithology number
        result = cell(n,1);
        name = input('Lithology Name: ','s');
        axes(ax);
        valinf = 1e5;
        xrange0 = get(gca,'xlim');
        set(gca,'ylimmode','manual');
        while(~isempty(name))
            disp('Start picking, left pick right delete q finish...');
            xrange = xrange0;
            top = [];
            bottom = [];
            h_line(1) = line('Color','r'); % top line
            h_line(1).XData = [];
            h_line(1).YData = [];
            h_line(2) = line('Color','b'); % bottom line
            h_line(2).XData = [];
            h_line(2).YData = [];
            m = 1; % 1 for top 2 for bottom other for ascii key
            % loop for this lithology
            while true
                [xpick,~,button] = ginput(1);
                if button~=double('q')
                    if button<4 && button>0
                        % correspond to the nearest x
                        [~,temp] = min(abs(xpick-x));
                        xpick = x(temp);
                        if button==1
                            % left click
                            switch m
                                % pick
                                case 1
                                    top = [top,xpick];
                                    h_line(1).XData = [h_line(1).XData,xpick,xpick,nan];
                                    h_line(1).YData = [h_line(1).YData,-valinf,valinf,nan];
                                    m = 2;
                                case 2
                                    bottom = [bottom,xpick];
                                    h_line(2).XData = [h_line(2).XData,xpick,xpick,nan];
                                    h_line(2).YData = [h_line(2).YData,-valinf,valinf,nan];
                                    m = 1;
                            end
                        elseif button==3
                            % right click delete, not completed yet
                            % seach for tops
                           temp = xpick - top;
                           temp(temp<0) = inf;
                           [~,temp] = min(temp);
                           try
                               bottom(temp) = [];
                               top(temp) = [];
                               for iter = 1:2
                                   h_line(iter).XData(3*(temp-1)+1:3*(temp-1)+3) = [];
                                   h_line(iter).YData(3*(temp-1)+1:3*(temp-1)+3) = [];
                               end
                           catch ME
                               if strcmp(ME.identifier,'MATLAB:subsdeldimmismatch')
                                   warning('Cannot delete: not in a block!');
                               end
                           end 
                        end
                    else
                        switch button
                            case double('d')
                                % right move
                                set(gca,'xlim',xrange+diff(xrange)/3);
                            case double('s')
                                % left move
                                set(gca,'xlim',xrange-diff(xrange)/3);
                            case double('w')
                                % shrink
                                set(gca,'xlim',xrange+[diff(xrange)/5,-diff(xrange)/5]);
                            case double('e')
                                % expand
                                set(gca,'xlim',xrange+[-diff(xrange)/5,diff(xrange)/5]);
                        end
                    end
                    xrange = get(gca,'xlim');
                else
                    % when q is pressed, do not do any thing when
                    % bottom is not picked yet
                    if m==1
                        % indicate bottom is picked
                        break;
                    end
                end
            end
            result{n}.top = top;
            result{n}.bottom = bottom;
            xflag = zeros(size(x));
            thickness = zeros(length(top),1);
            % calculate flag and thickness
            for k = 1:length(top)
                xflag(x>=top(k)&x<=bottom(k)) = 1;
                thickness(k) = -top(k) + bottom(k);
            end
            result{n}.xflag = xflag;
            result{n}.thickness = thickness;
            name = input('Lithology Name: ','s');
            % delelet the auxillary lines
            delete(h_line);
            n = n + 1;
        end
        disp('quit')       
    case 'seismic'
        button = 1;
        h_line = line(...
            'color','r',...
            'marker','o');
        h_line.XData = [];
        h_line.YData = [];
        n = 1;
        result.pick = [];
        result.pick_interp = [];
        xpick = [];
        ypick = [];
        while button~=3
            [xpick(n),ypick(n),button] = ginput(1);
            h_line.XData = [h_line.XData, xpick(n)];
            h_line.YData = [h_line.YData, ypick(n)];
            n = n + 1;
        end
        % interpolate
        result.pick_interp = interp1(xpick,ypick,x,'linear','extrap');
        h_line.XData = x;
        h_line.YData = result.pick_interp;
        h_line.Marker = 'none';
end


end

