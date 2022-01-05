function PlotSetup(position, x_label, y_label, x_lim, legend_show)
    % Set up plots in my personal way
    % position : int
    %   Number from 1 to 4 for the 4 quadrants of the screen
    % x_label : str
    %   String with the label of the x axis
    % y_label : str
    %   String with the label of the y axis
    % x_lim : double or array with 2 entries
    %   Min and max x axis values (if one number, used for min and max)
    % legend_show : bool
    %   Show or not the legend (default: true)
    
    if nargin < 5
        legend_show = true;
    end
    
    if position == 1 % top left
        pos = [100, 550, 700, 450];
    end
    if position == 2 % bottom left
        pos = [100, 0, 700, 450];
    end
    if position == 3 % top right
        pos = [1000, 550, 700, 450];
    end
    if position == 4 % bottom right
        pos = [1000, 0, 700, 450];
    end
    
    if length(x_lim) == 1
        x_max = x_lim;
        x_min = -x_lim;
    else
        x_max = x_lim(2);
        x_min = x_lim(1);
    end
    
    figure('Position', pos);
    hold on;
    xlabel(x_label);
    ylabel(y_label);
    grid on;
    if legend_show
        legend show;
    end    
    legend('Location', 'southeast');
    ax = gca;
    set(gca,'GridLineStyle','--');
    xlim([x_min x_max]);
%     ax.XTick = x_min:2:x_max;
    
    
    % xlim([-9 9]);ax.XTick = -9:2:9;
    % ylim([-1600 1600]);ax.YTick = -1600:500:1600;
    % set(gca,'FontSize',12);set(gca,'FontName','Times New Roman');
end



