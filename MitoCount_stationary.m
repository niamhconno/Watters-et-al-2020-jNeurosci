%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% This function
%%%     --Plots kymograph from .mat file. 
%%%     --Identifies potentially stationary (identified according to 
%%%         conditions defined within function!)
%%%     --Asks user to verify all objects identified as stationary [red]
%%%         (or close to stationary [cyan])  
%%%     --Plots all confirmed stationary objects as red on kymograph
%%% NB. If you enter a value outside those recognised by program (1 or 2),
%%% the algorithm will keep its original decision (red objects
%%% classified as stationary, cyan objects classified as non-stationary)
%%%
%%% Example function call: 
%%% MitoCount_stationary('sample_cellprofiler_output1')
%%% Press Ctrl+C at any time to quit
%%%
%%% If you use this code, please cite:
%%% Watters, Connolly et al., (2020) J Neurosci
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function MitoCount_stationary(filename)

% load .mat file
load(filename)  
% Clear any variables related to previous stationary analysis
clear stationary*           
disp('===================================================')
disp('===================================================')
fprintf('Filename: %s\n',filename)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Define range (distance) and % of images within which the objects 
%%% (mitochondria) can move and still be defined as 'stationary'. 
%%% If mitochondrium is present within the defined range in greater than x% 
%%% of images then it is labelled as stationary
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Range = 1% of total image width 
% Uses most distal mitochondria (over entire expt) to define image width.
% This was necessary, as although the um width of each image was constant, 
% the pixel-width varied between images/experiments 
% Nevertheless, it should be noted that this will affect experiments with
% e.g. very sparse mitochondria that do not move (so no mitochondria may 
% reach image edge throughout experiment)
% For Watters et al each image ~=354 um, so range ~= 3.5 um 
x_range = mito_most_dist * (1 / 100) ;   
fprintf('Image width: %i\n', mito_most_dist)

% Object must exist in 90% of timepoints
stationary_percent = 0.9;  
disp('To be considered stationary, object must remain within')
fprintf('+/- %0.2i of starting coordinate for %0.2i percent of the image\n', x_range, stationary_percent)
disp('Objects just outside this range will be coloured cyan and warnings will be generated')
input('Is this ok? (Enter or Ctrl+C to quit)')
disp('===================================================')

% Generate warnings if:
% Object exists in 70-90% of timepoints
stationary_percent_outside = 0.7;    
% > 1 object is detected within range for 50% of the timepoints
multiple_percent = 0.5;     
% Object is not found for 30 consecutive images
consec_images = 30;             

% Number of timepoints in each time interval of kymograph
kymo_size = time_interval(2,1) - time_interval(1,1);        

% For each time interval...
for m = 1:size(time_interval,2)  
    
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%% Draw kymograph
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   disp('===================================================')
   fprintf('\nTime Interval %i (Images %i - %i):\n',m,...
       time_interval(1,m),time_interval(2,m));

   start_row = time_interval(1,m);
   finish_row = time_interval(2,m);
    
   temp_name = strcat(filename,'_Stationary_TimeInt',num2str(m),...
       '_Images',num2str(start_row),'-',num2str(finish_row),'.fig');
   temp = 1;
    
   % If the kymograph already exists, check if user wants to overwrite
   if exist(temp_name,'file') == 2
       fprintf('The kymograph for this time interval already exists.\n')
       temp = input('Do you want to redo analysis and overwrite? (yes = 1 / no = 2) ');
   end
     
   if temp == 1
       disp('============================================================')
       disp('Please note:')
       disp('If stationary analysis was performed previously, and you are now redoing it,')
       disp('stationary variables may need to be deleted. This should be done manually.')
       disp('Please check protocol!')
       input('Press enter to continue')
       disp('============================================================')
       fprintf('Drawing...\n\n')
       
       % Draw kymograph, return y-axis
       [h,y_axis] = MitoCount_draw_kymograph(filename,m);  
       
       % # cols in x_center minus those that are NaN 
       temp_size = size(x_center(start_row,:),2) - sum(isnan(x_center(start_row,:))); 
       
       % Define (and zero) matrices for each time interval
       stationary_check = zeros(size(time_interval,2),temp_size);        
       multiple_check = zeros(1,temp_size);
       consec_empty = zeros(1,temp_size);
       temp_consec_empty = 0;
        
       % First row of 'check' = coordinates in FIRST image of time interval
       stationary_check(1,:) = x_center(start_row,1:temp_size);  

       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       %%% Populate stationary_check variable with 1 if no object exists
       %%% at that coordinate, or 0 if an object exists
       %%% This determines if that object will be checked to see if it
       %%% remains stationary for the time interval
       %%% Note that this means that an object that appears after the first
       %%% time interval will not be considered
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       % For each coordinate in FIRST image of time interval
       for i = 1:temp_size  
           % Define range around object coordinates within which object can 
           % 'move', based on coordinates of object in first image only
           min_x(i) = x_center(start_row,i) - x_range;        
           max_x(i) = x_center(start_row,i) + x_range;        
            
           % for each time-point within time interval
           for j = 1: finish_row - start_row  
            
            % If coordinate = NaN (i.e. no object at that coordinate, just 
            % sparse matrix value)...
            if isnan(x_center(start_row,i))
                % ...assign 1 (1 => mitochondria is not there!)
               stationary_check(j+1,i) = 1;    
            else
                %%% Identify stationary mito based on x coordinates.
                %%% Return index of first object within defined range
                %%% at that timepoint. Return empty matrix if no
                %%% object exists at that timepoint. Therefore
                %%% stationary_check is populated with 1 if no object
                %%% exists, or 0 if an object exists
                temp = find(x_center(start_row+j,:) < max_x(i) & x_center(start_row+j,:) > min_x(i)); 
                stationary_check(j+1,i) = isempty(temp);    
                
                % Count occasions where >1 object exists within allowed range
                if size(temp,2) > 1       
                    multiple_check(i) = multiple_check(i) + 1;
                end
                
                % Count consecutive images where no object is found
                if isempty(temp)            
                    temp_consec_empty = temp_consec_empty + 1; 
                else temp_consec_empty = 0;
                end
                
                % If no object is found for 30 images, generate warning below
                if temp_consec_empty > consec_images
                    consec_empty(i) = 1;     
                end                        
            end
           end
       end
    
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       %%% Invert stationary_check so that '1's correspond to the presence of 
       %%% an object in that image at that coordinate (easier to comprehend)
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       stationary_check(stationary_check ==0) = 2; %Convert 0s to 2
       stationary_check(stationary_check ==1) = 0; %Convert 1s to 0
       stationary_check(stationary_check ==2) = 1; %Convert 2s back to 1
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       %%% Determine if objects are conditionally classified as stationary, 
       %%% and verify with user. 
       %%% --Confirm all objects classified as stationary (red) 
       %%% --Check objects 'almost' classified as stationary (cyan). 
       %%% --Generate warnings when appropriate
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
       % First row contains object coordinates
       stationary_decision(m + (m-1),1:size(stationary_check,2)) = stationary_check(1,:);       
       % any remaining cols = NaN
       stationary_decision(m + (m-1),(size(stationary_check,2)+1):end) = NaN;        
       
       % for each object
       for i = 1:size(stationary_check,2)  
          % If the object is NaN - ignore!
          if isnan(stationary_check(1,i))    
          else
             % Leftmost point of rectangles to be drawn around objects
             x_rect = stationary_check(1,i) - x_range; 
             width_rect = x_range * 2; % width of rectangle 
       
             temp_string = '';
       
             %%% Identify indices and values of objects that fall within 
             %%% range (for subsequent plotting)
             j = start_row:finish_row;
             % Indices (row,col) of values that are within range
             [row,col] = find(x_center(j,:) < max_x(i) & x_center(j,:) > min_x(i)); 
             % Need to correct for whichever time interval you're looking at
             row2 = row + (start_row-1);   
             % Obtains linear index for the matrix (from row,col get one value)
             temp_index = sub2ind(size(x_center),row2,col);    
             % Values of x_center that are within range
             temp_x = x_center(temp_index);  
             temp_y = y_axis(row); 
           
             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
             %%% If object will be highlighted and investigated by user, 
             %%% print appropriate warnings to screen
             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
             if sum(stationary_check(2:end,i)) > kymo_size * stationary_percent_outside
                
                % If adjacent object is within range of the current object, 
                % generate warning (don't do it for the first or last objects, 
                % as you're checking both sides of the objects)
                if i > 1 && i < size(stationary_check,2)    
                    if (stationary_check(1,i+1) - stationary_check(1,i)) < width_rect ...
                            || (stationary_check(1,i) - stationary_check(1,i-1)) < width_rect
                        fprintf('Warning 1!! Two objects in close proximity - this may affect classification.\n')
                        temp_string = ('1');   % Text to be added on scatter plot
                    end
                end
            
                % If there are a high number of instances where > 1 object 
                % is detected within allowed range
                if multiple_check(i) > kymo_size * multiple_percent
                    fprintf('Warning 2!! Many instances where >1 object is detected within range for this object.\n')
                    temp_string = char(temp_string,'2');
                end
            
                % If object was not detected for a period of consecutive
                % images, display warning
                if consec_empty(i) == 1
                     fprintf('Warning 3!! There is a period (> %i images) during which this object is not detected.\n',consec_images)
                     temp_string = char(temp_string,'3');
                end
             end
       
             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
             %%% If object is identified within range in > 90% of timepoints
             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
             % (Ignore the first row as this includes object coordinates!)
             if sum(stationary_check(2:end,i)) > kymo_size * stationary_percent   
                % Surround current object with red diamond
                scatter(stationary_check(1,i),y_axis(1,i),100,'dr')  
                                     
                %%% Colour all objects that fall within this range 
                %(:) turns matrix into vector (all cols under eachother). scatter can only plot vectors
                temp_h = scatter(temp_x(:),temp_y(:),'.g');      
                    
                % Write warnings (numbers) onto plot above (+20) top most point
                text(stationary_check(1,i)-5, y_axis(1,i)+25, temp_string);       
           
                fprintf('This object (highlighted green) was automatically classified as stationary.\n')
                temp = input('Is this correct? (yes = 1 / no = 2)','s');
                
                switch temp
                   case '2'
                      %%% If user enters '2', stationary_decision = 0 (moving), &
                      %%% colour object blue
                      fprintf('Object is now classified as moving. \n\n')
                      stationary_decision(m*2,i) = 0;
                      scatter(stationary_check(1,i),y_axis(1,i),100,'db','filled','MarkerEdgeColor','k')
                      % Turn track back to black (just delete green points)
                      delete(temp_h);         
                      % Draw rectangle of appropriate width around object
                      rectangle('Position',[x_rect, 0, width_rect, kymo_size],'EdgeColor','r');    
                      
                    otherwise
                      %%% If user enters '1' (or anything else)
                      %%% stationary_decision = 1 (stationary) & colour 
                      %%% object red
                      fprintf('Object is now classified as stationary. \n\n')
                      stationary_decision(m*2,i) = 1;
                      % Add red diamond to top of plot
                      scatter(stationary_check(1,i),y_axis(1,i),100,'dr','filled','MarkerEdgeColor','k') 
                      delete(temp_h)
                      scatter(temp_x(:),temp_y(:),'.r')   % Colour track red
                end
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
                %%% elseif object is identified in > 80% of timepoints,
                %%% generate warning
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
             elseif sum(stationary_check(2:end,i)) > kymo_size * stationary_percent_outside 
           
                scatter(stationary_check(1,i),y_axis(1,i),100,'db')
                               
                %%% Colour all objects that fall within this range (cyan)
                temp_h = scatter(temp_x(:),temp_y(:),'.c');
                % Write warnings (numbers) onto plot above line
                text(stationary_check(1,i)-5, y_axis(1,i)+25, temp_string);       
           
                fprintf('This object (highlighted cyan) falls just outside stationary criteria (>80).\n')
                temp = input('Is this object stationary? (yes = 1 / no = 2)','s');
                
                switch temp
                   case '1'
                       %%% If user enters 1, stationary_decision = 1
                       %%% (stationary), & colour object red
                       fprintf('Object is now classified as stationary. \n\n')
                       stationary_decision(m*2,i) = 1;
                       scatter(stationary_check(1,i),y_axis(1,i),100,'dr','filled','MarkerEdgeColor','k')
                       delete(temp_h);
                       scatter(temp_x(:),temp_y(:),'.r')
                                
                    otherwise
                       %%% If user enters 2 (or anything else)
                       %%% stationary_decision = 0 (moving) & colour object blue
                       fprintf('Object is now classified as moving. \n\n')
                       stationary_decision(m*2,i) = 0;
                       scatter(stationary_check(1,i),y_axis(1,i),100,'db','filled','MarkerEdgeColor','k')
                       delete(temp_h);
                       % Draw rectangle of appropriate width around object
                       rectangle('Position',[x_rect, 0, width_rect, kymo_size],'EdgeColor','c');    
                end
            
                %%% If object is identified in < 80% of timepoints (i.e.
                %%% 'definitely' moving), stationary_decision = 0 (moving) & colour
                %%% object (in first image) blue
                else
                    stationary_decision(m*2,i) = 0;
                    scatter(stationary_check(1,i),y_axis(1,i),100,'db','filled','MarkerEdgeColor','k')
              end
          end
       end
        
       % Any remaining 0 value coordinates (except those in first column) 
       % will be changed to NaN
       for i = 1:2:(size(stationary_decision,1)-1)        
          temp = stationary_decision(i,:)==0;
          temp(:,1) = 0;
          stationary_decision(i,temp) = NaN;  % First
       end
            
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       %%% Count the total number of stationary mitochondria 
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       stationary_count(m,1) = sum(stationary_decision(m*2,:));
    
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       %%% Append stationary_decision and _count to workspace and save
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       save(filename, 'stationary_decision', 'stationary_count', '-append')
            
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       %%% Colour all stationary mitochondria traces red 
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       % Colour stationary mitochondria
       Mito_kymocolour_stationary_decision(filename,m,x_range,y_axis)     
       % Allow user to exit if they're not happy
       temp = input('Final figure: OK? (yes = 1 / no = 2)','s');    
       switch temp
          case '2'
             fprintf('Press CTRL+C to exit programme (figure and counts will not be saved)\n')
             input('(or just press any other key to continue [*figure will be saved*])');
       end
       fprintf('Stationary mitochondria in time interval %0.0i: %0.0i\n', ...
           interval_range(m), stationary_count(m))
        
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       %%% Save figure with filename & time_interval
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       saveas(h,temp_name)
       close all       % close figure so next one can be drawn
   else
       fprintf('No files have been changed, moving onto next time interval\n')
   end
end

disp('===================================================')
disp('===================================================')
disp('Stationary count complete for all time intervals!')
for i = 1:size(interval_range,2)
    fprintf('Time Interval %0.0i: %0.0i\n', interval_range(i), stationary_count(i))
end
end
