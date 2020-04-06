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
%%% Example function call: MitoCount('sample_cellprofile_output1')
%%% User input is required (marked as !!!)
%%%
%%% If you use this code, please cite:
%%% Watters, Connolly et al., (2020) J Neurosci
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function MitoCount(filename)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% User-specified inputs !!!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Size of time intervals into which data will be split 
% 450 = 450 images (timepoints) in each time interval
% With 4 images/s, this equates to 30 min time intervals
kymo_size = 450;  
fprintf('Each time interval is i% images.\n', kymo_size)

% Drug addition time (kymograph before this time are considered baseline)
% Set to 1 if there was no drug addition
drug_add = 450;
fprintf('Drug addition specified at image #%i.\n', drug_add)

% Identify columns for ImageNumber and AreaShape_Center_X within csv file
image_num_col = 1;
center_x_col = 2;
fprintf('Columns %i and %i will be read from imported file.\n', image_num_col, center_x_col)
input('Press Enter if these values are ok. Otherwise quit (Ctrl+C) and edit these values.')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Import .csv file as .mat file for further analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('################################################')
disp('Please ensure .mat file does not already exist!')
disp('If .mat file already exists, please quit (Ctrl+C) and delete file.')
input('Otherwise press enter')
disp('################################################')
% Read from 2nd row (ignore title row)
M = csvread(strcat(filename,'.csv'),1,0);  
    
% Sort M by first column (image number) and then by 4th column
% (x-coordinates)
M = sortrows(M,[image_num_col center_x_col]);  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Define time intervals (determined by kymo_size)
%%%% e.g. 1 -> 450; 451 ->900; 901 ->1350 etc. = 30 min time intervals
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
% If there is no drug addition AND ONLY kymo_size IMAGES
if drug_add == 1
   time_interval(1,1) = 1;
   time_interval(2,1) = kymo_size;
   % Number of time intervals
   interval_range = 1;
% If drug addition time is < kymo_size (450), 1st time interval will be shorter
elseif drug_add < kymo_size
   time_interval(1,1) = 1;
   time_interval(2,1) = drug_add;
   interval_range = 1:5; % 5 time intervals
else
   time_interval(1,1) = drug_add - (kymo_size-1);
   time_interval(2,1) = drug_add;
   interval_range = 1:5;
end
    
% Each experiment must fit into a maximum of 5 time intervals
for i = interval_range         
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
        disp('WARNING: ONLY 1 IMAGE INTERVAL')
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
%%% Combine shape areas also (for stationery calculations). Value at
%%% each location will correspond to size of mitochondria at the
%%% coordinate from the corresponding place in x_center
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x_center = zeros(max(M(:,image_num_col)),max(mito_total));  % Pre-define matrix for speed
    
% for all images
for i = 1 : max(M(:,image_num_col))   
    % for experiments where out of focus images were removed
    if isnan(mito_total(i))
       x_center(i,1:end) = NaN;        
    else
       % Put x-center coordinates and mito shape into one row 
       % x coordinates in 4th column of M
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
disp('--------------------------------------------------------------------')
disp('Filename:')
disp(filename)
fprintf('Number of images in .csv file: %i\n',size(mito_total,1))
fprintf('Split into %i time intervals\n', max(interval_range))
fprintf('Average # mitochondria at each timepoint: %0.2f\n',nanmean(mito_total))
disp('Mitochondrial counts are stored in mito_total variable of corresponding .mat file')
disp('--------------------------------------------------------------------')
        
end
