%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% This function imports .csv file into matlab, defines time intervals for 
%%% kymographs, and counts the mitochondria at each timepoint.
%%% The .csv file should follow a similar format to
%%% "sample_cellprofiler_output1.csv" 
%%% Columns required: ImageNumber, AreaShape_Center_X (=x-coordinate of
%%% center pixel of object)
%%%
%%% NB. If the .mat file already exists, please delete before executing this
%%% code, as otherwise the variables assigned above the file import will be 
%%% when the file is imported
%%%
%%% Example function call: 
%%% MitoCount_stationary_combined('sample_cellprofiler_output1', 450, 450, 1, 2)
%%%
%%% If you use this code, please cite:
%%% Watters, Connolly et al., (2020) J Neurosci
%%% DOI: 10.1523/JNEUROSCI.2067-19.2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function MitoCount_stationary(filename, kymo_size, drug_add, image_num_col, center_x_col)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Function inputs specified by user 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% kymo_size: 
% Size of time intervals into which data will be split 
% 450 = 450 images (timepoints) in each time interval (default)
% With 4 images/s, this equates to 30 min time intervals

% drug_add:
% Drug addition time (kymograph before this time are considered baseline)
% Set to 1 if there was no drug addition (default)

% image_num_col, center_x_col
% column indices for ImageNumber and AreaShape_Center_X within input file
% default values 1,2

% Funtion to check if inputs are integers 
isaninteger = @(x)isfinite(x) & x==floor(x);

% Assign default values to input parameters if they have not been specified
switch nargin
    case 0
        error('A filename must be entered for this function')
    case 1
        disp('Default values have been set for all input parameters')
        kymo_size = 450;  
        drug_add = 1;
        image_num_col = 1;
        center_x_col = 2;
    case 2
        disp('Default values have been set for some input parameters')
        if isaninteger(kymo_size)
        else; error('Input parameter kymo_size must be an integer') 
        end
        drug_add = 1;
        image_num_col = 1;
        center_x_col = 2;
    case 3
        disp('Default values have been set for some input parameters')
        if isaninteger(kymo_size) && isaninteger(drug_add)
        else; error('Input parameters kymo_size, drug_add must be integers') 
        end
        image_num_col = 1;
        center_x_col = 2;
    case 4
        disp('Default values have been set for some input parameters')
        if isaninteger(kymo_size) && isaninteger(drug_add) && isaninteger(image_num_col)
        else; error('Input parameters kymo_size, drug_add, image_num_col must be integers') 
        end
        center_x_col = 2;
    case 5
        if isaninteger(kymo_size) && isaninteger(drug_add) && isaninteger(image_num_col) && isaninteger(center_x_col)
        else; error('Input parameters kymo_size, drug_add, image_num_col, center_x_col must be integers') 
        end
    otherwise
        error('A maximum of 5 input arguments can be accepted for this function')
end

disp('###############################')
disp('If you use this code, please cite:')
disp('Watters, Connolly et al., (2020) J Neurosci. DOI: 10.1523/JNEUROSCI.2067-19.2020')
disp('###############################')
disp('###############################')
disp('Before proceeding, please ensure .mat file does not already exist!')
disp('If .mat file already exists, please quit (Ctrl+C) and delete file.')
input('Otherwise press enter')
disp('###############################')

% Verify that input values are ok
fprintf('---Each time interval is %i images (''kymo_size'').\n', kymo_size)
fprintf('---Drug addition is specified at image #%i (''drug_add'').\n', drug_add)
if drug_add == 1; disp('This implies no drug addition'); end
fprintf('---Columns %i and %i will be read from imported file.\n', image_num_col, center_x_col)
input('Press Enter if these values are ok. Otherwise quit (Ctrl+C) and edit function inputs.')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Import .csv file and save as .mat file for further analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Read from 2nd row (ignore title row)
M = csvread(strcat(filename,'.csv'),1,0);  
    
% Sort M by image number column and then by x-coordinate column
M = sortrows(M,[image_num_col center_x_col]);  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Define time intervals (determined by kymo_size)
%%%% e.g. 1 -> 450; 451 ->900; 901 ->1350 etc. = 30 min time intervals
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Specify number of intervals (maximum of 5)
% = number of images/number of images per kymograph
num_intervals = min(5,max(M(:,image_num_col)/kymo_size));

% If there is no drug addition AND ONLY kymo_size IMAGES (i.e. 1 kymograph)
if drug_add == 1
   time_interval(1,1) = 1;
   time_interval(2,1) = kymo_size;
   % Number of time intervals/kymographs
   interval_range = 1;
% If drug addition time is < kymo_size (450), 1st time interval will be shorter
elseif drug_add < kymo_size
   time_interval(1,1) = 1;
   time_interval(2,1) = drug_add;
   interval_range = 1:num_intervals; 
else
   time_interval(1,1) = drug_add - (kymo_size-1);
   time_interval(2,1) = drug_add;
   interval_range = 1:num_intervals;
end
    
% For each time interval (have already assigned 1st interval)
for i = 1:num_intervals-1         
% Only define time intervals if images exist for ALL images within
% that time interval - otherwise the time interval is incomplete 
    if time_interval(2,i) + kymo_size <= max(M(:,image_num_col))    
        time_interval(1,i+1) = time_interval(2,i) +1;
        time_interval(2,i+1) = time_interval(1,i+1) + kymo_size-1;
% Generate last time interval if there are 375 or more images (supposed to be 450)
    elseif time_interval(2,i) + (kymo_size-100) <= max(M(:,image_num_col))        
        time_interval(1,i+1) = time_interval(2,i) +1;
        time_interval(2,i+1) = max(M(:,image_num_col)); 
    else
        time_interval(1,i+1) = max(M(:,image_num_col));
        disp('WARNING: ONLY 1 TIME INTERVAL')
    end        
end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Identify most proximal and distal mitochondria
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mito_most_prox = min(M(:,center_x_col));  
mito_most_dist = max(M(:,center_x_col));   
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Count mitochondria 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For all images [max(M(:,image_num_col))]
for i = 1:max(M(:,image_num_col))        
    % Assign 1 to all rows associated with this image, 
    % and 0 otherwise - to count separately for each timestep
    timestep = M(:,image_num_col) == i;      
        
    % Count the number of mitochondria at each timestep 
    mito_total(i,:) = sum(timestep);
end                        
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Combine all x center coordinates into one matrix, where
%%% each row is a timepoint and each column is a coordinate
%%% corresponding to the centre of each object at that timepoint.
%%% We will use this matrix to draw kymograph.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pre-define matrix for speed
x_center = zeros(max(M(:,image_num_col)),max(mito_total));  
    
% for all images
for i = 1 : max(M(:,image_num_col))   
    % for experiments where out of focus images were removed
    if isnan(mito_total(i))
       x_center(i,1:end) = NaN;        
    else
       % Put x-center coordinates and mito shape into one row 
       x_center(i,1:mito_total(i)) = M(M(:,image_num_col) == i,center_x_col)'; 
    end
end
      
temp = x_center==0; % temp has 1 wherever x_center =0
temp(:,1) = 0;      % remove any 1s from first column (as coordinate of first object may =0)
x_center(temp) = NaN;           % Change 0s to NaN 
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Clear redundant variables and save workspace
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear i j m timestep temp* %file_list
save(filename)
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Print variables to screen 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('###############################')
disp('###############################')

fprintf('Filename: %s\n', filename)
fprintf('Number of images in .csv file: %i\n',size(mito_total,1))
fprintf('Split into %i time intervals (%i images in each interval)\n',...
    max(interval_range), kymo_size)
fprintf('Average # mitochondria at each timepoint: %0.2f\n',nanmean(mito_total))
disp('Mitochondrial counts are stored in mito_total variable of .mat file')
fprintf('\n')
disp('Now it''s time to count the stationary mitochondria!')
input('Press Enter if you are happy to continue. (Ctrl+C will quit)')

% Function to identify and classify stationary mitochondria
stationary(filename)

end
