%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Colours mitochondria based on stationary decision
%%% Requires kymograph to be previously plotted
%%% 
%%% Called from within MitoCount_stationary.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Mito_kymocolour_stationary_decision(filename,plot_time_int,x_range,y_axis)

%%% plot_time_int: time interval that is plotted (1 = baseline)
%%% x_range: range within which a mitochondria can move while still
%%% considered stationery. x_range must be used to classify
%%% the mitochondria as otherwise you're not plotting the mitochondria that
%%% are classified as stationary!
%%% kymo_size: size (in images) of the plotted kymograph (generally 450)
%%% y_axis: y-axis used to plot the kymograph (450 - 1)

load(filename)

% Row in stationary_decision variable that contains stationary decision!
temp_row = plot_time_int * 2;   

start_row = time_interval(1,plot_time_int);
finish_row = time_interval(2,plot_time_int);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Obtain coordinates of stationary mitochondria
%%% Coordinates are first row of stationary_decision where 2nd row = 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
stationary_coords = stationary_decision(temp_row - 1, stationary_decision(temp_row,:) == 1);  
stationary_coords(stationary_coords == 0) = NaN;

% For every stationary object
for i = 1:size(stationary_coords,2) 
   min_x(i) = stationary_coords(1,i) - x_range;       
   max_x(i) = stationary_coords(1,i) + x_range;
       
   %%% For each timepoint (start_row to finish_row), find coordinates 
   %%% of objects that lie within 'stationary' range (location only) 
   %%% and colour them red
   
   j = start_row:finish_row;
   % Indices (row,col) of values that are within range
   [row col] = find(x_center(j,:) < max_x(i) & x_center(j,:) > min_x(i)); 
      
   % Need to correct for whichever time interval you're looking at
   row2 = row + (start_row-1); 
   % Obtains linear index for the matrix (from row,col get one value)
   temp_index = sub2ind(size(x_center),row2,col);    
   
   % Values of x_center that are within range
   temp_x = x_center(temp_index);  
   temp_y = y_axis(row);   
   %(:) turns matrix into vector (all cols under eachother). scatter can only plot vectors
   scatter(temp_x(:),temp_y(:),'.r')      
   
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% At the first timepoint, plot stationary mitochondria as red 
% and moving mitochondria as blue
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
stationary = stationary_decision(temp_row,:) == 1;
moving = stationary_decision(temp_row,:) == 0;

scatter(stationary_decision(temp_row-1,moving),y_axis(1,moving),100,'db','filled','MarkerEdgeColor','k')
scatter(stationary_decision(temp_row-1,stationary),y_axis(1,stationary),100,'dr','filled','MarkerEdgeColor','k')

end
