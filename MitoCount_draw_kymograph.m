function [h,y_axis] = MitoCount_draw_kymograph(filename,m)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%% Plot kymograph within timeframe specified. 
%%% Here, 'kymograph' is the x-cordinates of the centre of each object
%%% throughout the time-course of the imaging experiment (i.e. no
%%% information on object size, shape etc.)
%%%
%%% Filename is the .mat file name of experiment for which you
%%% wish to plot the kymograph (created from .csv file using MitoCount)
%%% m is the time interval you want to plot (1 = baseline). 
%%%
%%% This function is called from within MitoCount_stationary.m but can also
%%% be called independently. 
%%%
%%% Sample function call: 
%%% MitoCount_draw_kymograph('sample_cellprofiler_output1', 1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load(filename)

% Time interval to plot (1 = baseline)
plot_time_int = m;      
start = time_interval(1,plot_time_int) - 1;

kymo_size = time_interval(2,plot_time_int) - time_interval(1,plot_time_int) +1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Define y-axis. Row size = number of timepoints you want to plot. Column
%%% size = size of x_center (i.e. number of objects) x_axis when plotting 
%%% (x_center) should be the same size
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
y_axis = ones(kymo_size,size(x_center,2));  

for i = 1:kymo_size
    % y-axis goes from kymo_size (450) to 1
    y_axis(i,:) = y_axis(i,:) * kymo_size-i;     
end

% Get screensize so can draw the kymograph the width of the screen...
scrsz = get(0,'ScreenSize');            
h = figure('Position',[10 scrsz(4)/3 (scrsz(3)/1.05)*(mito_most_dist/1230) ...
    scrsz(4)/2],'Name',strcat(filename,'_',num2str(plot_time_int)));
hold on

i = 1:kymo_size;
% Plot rows of x_center within selected time period (start -> start + kymo_size)
temp_x = x_center(i + start,:);     
% Plot corresponding values from y-axis
temp_y = y_axis(i,:);                    

% Scatter points onto kymograph (all points the same size)
scatter(temp_x(:),temp_y(:),'.k')       

% Set limits of x axis
xlim([0,mito_most_dist]);    

ylabel('Inverted image number (i.e. time goes from top to bottom)')
xlabel('X coordinates')

end
