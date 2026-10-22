%%
volS_clean = volS;

%% DCSCAN Duct
[idx_x, idx_y, idx_z] = ind2sub(size(volS), find(volS == 2));
points = [idx_x, idx_y, idx_z];

epsilon = 22;
minPts = 5;
labels = dbscan(points, epsilon, minPts, 'Distance', 'euclidean');

num_clusters = max(labels);
duct_clusters = cell(1, num_clusters);
for k = 1:num_clusters
    duct_clusters{k} = points(labels == k, :);
end

%%
labeled_vol = ones(size(volS));

for k = 1:num_clusters
    coords = duct_clusters{k};
    
    if ~isempty(coords)
        linear_idx = sub2ind(size(labeled_vol), coords(:,1), coords(:,2), coords(:,3));
        
        labeled_vol(linear_idx) = k + 1;
    end
end
%% 
volumeViewer(labeled_vol,VolumeType="labels")

%%
volumeViewer(volS,VolumeType="labels")

%%
cmap=[233 244 163;... % Acinar; yellow
      72 100 174;... % Duct; Blue
      126 204 162;... % LG PanIN; green
      253 205 230;... % HG PanIN; pink
      244 109 64;... % Cancer; orange
      251 217 128;... % Vascular; Orange
      206 212 218;... % stromal; grey
      255 255 255;... % Background; white
      202 175 212;... % Neural purple
      0 0 0
      ];

%% Duct spatial visualization + annotation
max_layer_list = NaN(num_clusters, 1);

cmap=[251 217 128;... % Acinar; yellow
      72 100 174;... % Duct; blue
      126 204 162;... % LG PanIN; green
      253 205 230;... % HG PanIN; pink
      244 109 64;... % Cancer; orange
      233 244 163;... % Vascular; orange
      206 212 218;... % Stromal; grey
      255 255 255;... % Background; white
      202 175 212;... % Neural; purple
      0 0 0
      ];
cmap = cmap / 255;

% --- 2. Loop through each cluster ---
num_clusters = numel(duct_clusters);
fprintf('Starting generation and saving of cross-sectional images...\n');

for k = 1:num_clusters
    coords = duct_clusters{k};
    
    % Skip empty cells
    if isempty(coords)
        continue;
    end
    
    % A. Find the layer containing the largest number of points
    %    (assuming the third column represents the Z axis)
    z_vals = coords(:, 3);
    z_max = mode(z_vals); % mode directly returns the most frequently occurring value
    
    max_layer_list(k) = z_max;

    % B. Extract the corresponding layer (2D slice) from volS
    slice_2D = volS(:, :, z_max);
    
    % C. Find all points from the current cluster located on this layer
    % Extract the X and Y coordinates of rows where Z equals z_max
    pts_on_slice = coords(z_vals == z_max, :); 
    
    % Convert X and Y coordinates to linear indices in the 2D matrix
    idx_2D = sub2ind(size(slice_2D), pts_on_slice(:, 1), pts_on_slice(:, 2));
    
    % D. Change the values of these specific points to 10
    slice_2D(idx_2D) = 10;
    
    % --- 3. Plot and save ---
    % Create an invisible figure to prevent numerous windows from opening
    % and slowing down the computer
    fig = figure('Visible', 'off'); 
    
    % Display the image
    imagesc(slice_2D);
    
    % Apply the custom colormap and map label values to the defined colors
    colormap(cmap);
    clim([1, 11]); % For older MATLAB versions, replace clim with caxis
    
    % Adjust the display aspect ratio and hide the axes
    axis image; 
    axis off;   
    title(sprintf('Cluster %d - Max Area Layer Z=%d', k, z_max));
    
    % Save as a JPG image
    output_dir = 'F:\P1_3D\10x\classification_05_07_2025\check_Duct\';
    file_name = sprintf('Layer_%d_Cluster_%03d.jpg', z_max, k);
    filename = fullfile(output_dir, file_name);
    exportgraphics(gca, filename, 'Resolution', 300); % Save at 300 DPI
    
    % Close the current figure to release memory
    close(fig); 
end

fprintf('All images have been saved to the output folder!\n');
%%
% 1. Define the source and destination folder paths
src_dir = 'F:\P1_3D\10x\classification_05_07_2025\check_classification\';
dst_dir = 'F:\P1_3D\10x\classification_05_07_2025\check_Duct\';

% Create the destination folder if it does not exist
if ~exist(dst_dir, 'dir')
    mkdir(dst_dir);
    fprintf('Created destination folder: %s\n', dst_dir);
end

% 2. Find all "PDAC_*.jpg" files in the source directory
file_pattern = fullfile(src_dir, 'PDAC_*.jpg');
files_info = dir(file_pattern);

if isempty(files_info)
    error('No files starting with PDAC_ and ending with .jpg were found in the source folder!');
end

% Store all file names in a cell array
file_names = {files_info.name}; 

% 3. Extract the numbers from the file names and sort numerically
% Use regexp to extract consecutive digits and convert them to double-precision numbers
layer_nums = cellfun(@(x) str2double(regexp(x, '\d+', 'match', 'once')), file_names);

% Sort the extracted layer numbers in ascending order
[~, sort_idx] = sort(layer_nums);
sorted_file_names = file_names(sort_idx); % Sorted list of file names

fprintf('Successfully read and sorted %d image files.\n', length(sorted_file_names));

% 4. Clean max_layer_list (remove duplicates and exclude NaN)
valid_indices = max_layer_list(~isnan(max_layer_list));
unique_indices = unique(valid_indices);

fprintf('Preparing to copy %d selected images...\n', length(unique_indices));
copy_count = 0;

% 5. Select and copy files using the sorted indices
for i = 1:length(unique_indices)
    % target_idx is the index of the Nth image to select
    target_idx = unique_indices(i);
    
    % Check that the index is positive and does not exceed the number of available images
    if target_idx > 0 && target_idx <= length(sorted_file_names)
        
        % Get the original name of the file at target_idx in the sorted list
        src_file_name = sorted_file_names{target_idx};
        
        % Construct the full source file path
        src_file = fullfile(src_dir, src_file_name);
        
        % --- Construct the destination file path ---
        % Option A: Copy using the original file name
        % dst_file_name = src_file_name; 
        
        % Option B: Rename to Layer_XXX.jpg when copying (as previously requested)
        dst_file_name = sprintf('Layer_%d.jpg', target_idx);
        % --------------------------
        
        dst_file = fullfile(dst_dir, dst_file_name);
        
        % Copy the file
        copyfile(src_file, dst_file);
        copy_count = copy_count + 1;
        
    else
        warning('Index %d is out of range (%d images available); skipping.', target_idx, length(sorted_file_names));
    end
end

fprintf('Copying completed! Successfully processed %d files.\n', copy_count);

%%
duct_to_lg = [12 16];

%% Reclassify Duct as LG PanIN
for idx = duct_to_lg
    % 1. Get the coordinate matrix (N x 3) at the current index
    coords = duct_clusters{idx}; 
    
    if ~isempty(coords)
        % 2. Extract the x, y, z coordinates
        % Note: MATLAB matrix indexing typically uses VolS(row, col, page)
        % Assume the coordinate columns represent rows, columns, and layers, respectively:
        r = coords(:, 1);
        c = coords(:, 2);
        z = coords(:, 3);
        
        % 3. Convert 3D coordinates to linear indices
        % The first argument of sub2ind is the size of the target volume
        linIdx = sub2ind(size(volS_clean), r, c, z);
        
        % 4. Set the corresponding positions in VolS to 3
        volS_clean(linIdx) = 3;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DCSCAN LG PanIN
[idx_x, idx_y, idx_z] = ind2sub(size(volS), find(volS == 3));
points = [idx_x, idx_y, idx_z];

% Step 2: Apply DBSCAN (k-d tree-based implementation)
epsilon = 22;         % Maximum neighborhood radius (in voxels)
minPts = 5;          % Minimum number of neighbors for each core point
labels = dbscan(points, epsilon, minPts, 'Distance', 'euclidean');

% Step 3: Organize the clustering results into coordinate lists for each cluster
num_clusters = max(labels);
lg_clusters = cell(1, num_clusters);
for k = 1:num_clusters
    lg_clusters{k} = points(labels == k, :);
end

%%
labeled_vol = ones(size(volS));

% 2. Loop through all clusters and assign labels 2, 3, 4, ... (reserving 1 for the background)
for k = 1:num_clusters
    % Get the 3D coordinates of the kth cluster
    coords = lg_clusters{k};
    
    % If the cluster is not empty
    if ~isempty(coords)
        % Convert [x, y, z] coordinates to linear indices in the matrix
        linear_idx = sub2ind(size(labeled_vol), coords(:,1), coords(:,2), coords(:,3));
        
        % Assign a unique label (k + 1) to this cluster
        labeled_vol(linear_idx) = k + 1;
    end
end
%% 
volumeViewer(labeled_vol,VolumeType="labels")

%%
volumeViewer(volS,VolumeType="labels")

%%
cmap=[233 244 163;... % Acinar; yellow
      72 100 174;... % Duct; Blue
      126 204 162;... % LG PanIN; green
      253 205 230;... % HG PanIN; pink
      244 109 64;... % Cancer; orange
      251 217 128;... % Vascular; Orange
      206 212 218;... % stromal; grey
      255 255 255;... % Background; white
      202 175 212;... % Neural purple light purple
      0 0 0
      ];

%% LG PanIN spatial visualization + annotation correction
max_layer_list = NaN(num_clusters, 1);

cmap=[251 217 128;... % Acinar; yellow
      72 100 174;... % Duct; Blue
      126 204 162;... % LG PanIN; green
      253 205 230;... % HG PanIN; pink
      244 109 64;... % Cancer; orange
      233 244 163;... % Vascular; Orange
      206 212 218;... % stromal; grey
      255 255 255;... % Background; white
      202 175 212;... % Neural purple light purple
      0 0 0

      ];
cmap = cmap / 255;

% --- 2. Loop through each cluster ---
num_clusters = numel(lg_clusters);
fprintf('Starting generation and saving of cross-sectional images...\n');

for k = 1:num_clusters
    coords = lg_clusters{k};
    
    % Skip empty cells
    if isempty(coords)
        continue;
    end
    
    % A. Find the layer with the most points (assuming column 3 represents the Z-axis layer)
    z_vals = coords(:, 3);
    z_max = mode(z_vals); % mode directly returns the most frequently occurring value
    
    max_layer_list(k) = z_max;

    % B. Extract the corresponding layer (2D slice) from VolS
    slice_2D = volS(:, :, z_max);
    
    % C. Find all points from the current cluster on this layer
    % Extract the X and Y coordinates of rows where Z equals z_max
    pts_on_slice = coords(z_vals == z_max, :); 
    
    % Convert X and Y coordinates to linear indices in the 2D matrix
    idx_2D = sub2ind(size(slice_2D), pts_on_slice(:, 1), pts_on_slice(:, 2));
    
    % D. Set the values of these specific points to 10
    slice_2D(idx_2D) = 10;
    
    % --- 3. Plot and save ---
    % Create an invisible figure to avoid opening many windows and slowing down the computer
    fig = figure('Visible', 'off'); 
    
    % Display the image
    imagesc(slice_2D);
    
    % Apply the custom colormap and map values 0-10 to the 11 color rows
    colormap(cmap);
    clim([1, 10]); % For older MATLAB versions, replace clim with caxis
    
    % Adjust the display aspect ratio and hide the axes
    axis image; 
    axis off;   
    title(sprintf('Cluster %d - Max Area Layer Z=%d', k, z_max));
    
    % Save as a JPG image
    output_dir = 'F:\P1_3D\10x\classification_05_07_2025\check_LG\';
    file_name = sprintf('Layer_%d_Cluster_%03d.jpg', z_max, k);
    filename = fullfile(output_dir, file_name);
    exportgraphics(gca, filename, 'Resolution', 300); % Save at 300 DPI
    
    % Close the current figure to release memory
    close(fig); 
end

fprintf('All images have been saved to the current folder!\n');
%%
% 1. Define the source and destination folder paths
src_dir = 'F:\P1_3D\10x\classification_05_07_2025\check_classification\';
dst_dir = 'F:\P1_3D\10x\classification_05_07_2025\check_LG\';

% Create the destination folder if it does not exist
if ~exist(dst_dir, 'dir')
    mkdir(dst_dir);
    fprintf('Created destination folder: %s\n', dst_dir);
end

% 2. Find all "PDAC_*.jpg" files in the source directory
file_pattern = fullfile(src_dir, 'PDAC_*.jpg');
files_info = dir(file_pattern);

if isempty(files_info)
    error('No files starting with PDAC_ and ending with .jpg were found in the source folder!');
end

% Store all file names in a cell array
file_names = {files_info.name}; 

% 3. Extract the numbers from the file names and sort numerically
% Use regexp to extract consecutive digits and convert them to double-precision numbers
layer_nums = cellfun(@(x) str2double(regexp(x, '\d+', 'match', 'once')), file_names);

% Sort the extracted layer numbers in ascending order
[~, sort_idx] = sort(layer_nums);
sorted_file_names = file_names(sort_idx); % Sorted list of file names

fprintf('Successfully read and sorted %d image files.\n', length(sorted_file_names));

% 4. Clean max_layer_list (remove duplicates and exclude NaN)
valid_indices = max_layer_list(~isnan(max_layer_list));
unique_indices = unique(valid_indices);

fprintf('Preparing to copy %d selected images...\n', length(unique_indices));
copy_count = 0;

% 5. Select and copy files using the sorted indices
for i = 1:length(unique_indices)
    % target_idx is the index of the Nth image to select
    target_idx = unique_indices(i);
    
    % Check that the index is positive and does not exceed the number of available images
    if target_idx > 0 && target_idx <= length(sorted_file_names)
        
        % Get the original name of the file at target_idx in the sorted list
        src_file_name = sorted_file_names{target_idx};
        
        % Construct the full source file path
        src_file = fullfile(src_dir, src_file_name);
        
        % --- Construct the destination file path ---
        % Option A: Copy using the original file name
        % dst_file_name = src_file_name; 
        
        % Option B: Rename to Layer_XXX.jpg when copying (as previously requested)
        dst_file_name = sprintf('Layer_%d.jpg', target_idx);
        % --------------------------
        
        dst_file = fullfile(dst_dir, dst_file_name);
        
        % Copy the file
        copyfile(src_file, dst_file);
        copy_count = copy_count + 1;
        
    else
        warning('Index %d is out of range (%d images available); skipping.', target_idx, length(sorted_file_names));
    end
end

fprintf('Copying completed! Successfully processed %d files.\n', copy_count);

%%
lg_to_stroma = [7 9 15 18 20];
lg_to_cancer = [17, 24];

%% Reclassify LG PanIN as Cancer
for idx = lg_to_cancer
    % 1. Get the coordinate matrix (N x 3) at the current index
    coords = lg_clusters{idx}; 
    
    if ~isempty(coords)
        % 2. Extract the x, y, z coordinates
        % Note: MATLAB matrix indexing typically uses VolS(row, col, page)
        % Assume the coordinate columns represent rows, columns, and layers, respectively:
        r = coords(:, 1);
        c = coords(:, 2);
        z = coords(:, 3);
        
        % 3. Convert 3D coordinates to linear indices
        % The first argument of sub2ind is the size of the target volume
        linIdx = sub2ind(size(volS_clean), r, c, z);
        
        % 4. Set the corresponding positions in VolS to 5
        volS_clean(linIdx) = 5;
    end
end

%% Reclassify LG PanIN as Stroma
for idx = lg_to_stroma
    % 1. Get the coordinate matrix (N x 3) at the current index
    coords = lg_clusters{idx}; 
    
    if ~isempty(coords)
        % 2. Extract the x, y, z coordinates
        % Note: MATLAB matrix indexing typically uses VolS(row, col, page)
        % Assume the coordinate columns represent rows, columns, and layers, respectively:
        r = coords(:, 1);
        c = coords(:, 2);
        z = coords(:, 3);
        
        % 3. Convert 3D coordinates to linear indices
        % The first argument of sub2ind is the size of the target volume
        linIdx = sub2ind(size(volS_clean), r, c, z);
        
        % 4. Set the corresponding positions in VolS to 5
        volS_clean(linIdx) = 7;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DCSCAN HG PanIN
[idx_x, idx_y, idx_z] = ind2sub(size(volS), find(volS == 4));
points = [idx_x, idx_y, idx_z];

% Step 2: Apply DBSCAN (k-d tree-based implementation)
epsilon = 22;         % Maximum neighborhood radius (in voxels)
minPts = 5;          % Minimum number of neighbors for each core point
labels = dbscan(points, epsilon, minPts, 'Distance', 'euclidean');

% Step 3: Organize the clustering results into coordinate lists for each cluster
num_clusters = max(labels);
hg_clusters = cell(1, num_clusters);
for k = 1:num_clusters
    hg_clusters{k} = points(labels == k, :);
end

%%
labeled_vol = ones(size(volS));

% 2. Loop through all clusters and assign labels 2, 3, 4, ... (reserving 1 for the background)
for k = 1:num_clusters
    % Get the 3D coordinates of the kth cluster
    coords = hg_clusters{k};
    
    % If the cluster is not empty
    if ~isempty(coords)
        % Convert [x, y, z] coordinates to linear indices in the matrix
        linear_idx = sub2ind(size(labeled_vol), coords(:,1), coords(:,2), coords(:,3));
        
        % Assign a unique label (k + 1) to this cluster
        labeled_vol(linear_idx) = k + 1;
    end
end
%% 
volumeViewer(labeled_vol,VolumeType="labels")

%%
volumeViewer(volS,VolumeType="labels")

%%
cmap=[233 244 163;... % Acinar; yellow
      72 100 174;... % Duct; Blue
      126 204 162;... % LG PanIN; green
      253 205 230;... % HG PanIN; pink
      244 109 64;... % Cancer; orange
      251 217 128;... % Vascular; Orange
      206 212 218;... % stromal; grey
      255 255 255;... % Background; white
      202 175 212;... % Neural purple light purple
      0 0 0
      ];

%% HG PanIN spatial visualization + annotation correction
max_layer_list = NaN(num_clusters, 1);

cmap=[251 217 128;... % Acinar; yellow
      72 100 174;... % Duct; Blue
      126 204 162;... % LG PanIN; green
      253 205 230;... % HG PanIN; pink
      244 109 64;... % Cancer; orange
      233 244 163;... % Vascular; Orange
      206 212 218;... % stromal; grey
      255 255 255;... % Background; white
      202 175 212;... % Neural purple light purple
      0 0 0
      ];
cmap = cmap / 255;

% --- 2. Loop through each cluster ---
num_clusters = numel(hg_clusters);
fprintf('Starting generation and saving of cross-sectional images...\n');

for k = 1:num_clusters
    coords = hg_clusters{k};
    
    % Skip empty cells
    if isempty(coords)
        continue;
    end
    
    % A. Find the layer with the most points (assuming column 3 represents the Z-axis layer)
    z_vals = coords(:, 3);
    z_max = mode(z_vals); % mode directly returns the most frequently occurring value
    
    max_layer_list(k) = z_max;

    % B. Extract the corresponding layer (2D slice) from VolS
    slice_2D = volS(:, :, z_max);
    
    % C. Find all points from the current cluster on this layer
    % Extract the X and Y coordinates of rows where Z equals z_max
    pts_on_slice = coords(z_vals == z_max, :); 
    
    % Convert X and Y coordinates to linear indices in the 2D matrix
    idx_2D = sub2ind(size(slice_2D), pts_on_slice(:, 1), pts_on_slice(:, 2));
    
    % D. Set the values of these specific points to 10
    slice_2D(idx_2D) = 10;
    
    % --- 3. Plot and save ---
    % Create an invisible figure to avoid opening many windows and slowing down the computer
    fig = figure('Visible', 'off'); 
    
    % Display the image
    imagesc(slice_2D);
    
    % Apply the custom colormap and map values 0-10 to the 11 color rows
    colormap(cmap);
    clim([1, 10]); % For older MATLAB versions, replace clim with caxis
    
    % Adjust the display aspect ratio and hide the axes
    axis image; 
    axis off;   
    title(sprintf('Cluster %d - Max Area Layer Z=%d', k, z_max));
    
    % Save as a JPG image
    output_dir = 'F:\P1_3D\10x\classification_05_07_2025\check_HG\';
    file_name = sprintf('Layer_%d_Cluster_%03d.jpg', z_max, k);
    filename = fullfile(output_dir, file_name);
    exportgraphics(gca, filename, 'Resolution', 300); % Save at 300 DPI
    
    % Close the current figure to release memory
    close(fig); 
end

fprintf('All images have been saved to the current folder!\n');
%%
% 1. Define the source and destination folder paths
src_dir = 'F:\P1_3D\10x\classification_05_07_2025\check_classification';
dst_dir = 'F:\P1_3D\10x\classification_05_07_2025\check_HG\';

% Create the destination folder if it does not exist
if ~exist(dst_dir, 'dir')
    mkdir(dst_dir);
    fprintf('Created destination folder: %s\n', dst_dir);
end

% 2. Find all "PDAC_*.jpg" files in the source directory
file_pattern = fullfile(src_dir, 'PDAC_*.jpg');
files_info = dir(file_pattern);

if isempty(files_info)
    error('No files starting with PDAC_ and ending with .jpg were found in the source folder!');
end

% Store all file names in a cell array
file_names = {files_info.name}; 

% 3. Extract the numbers from the file names and sort numerically
% Use regexp to extract consecutive digits and convert them to double-precision numbers
layer_nums = cellfun(@(x) str2double(regexp(x, '\d+', 'match', 'once')), file_names);

% Sort the extracted layer numbers in ascending order
[~, sort_idx] = sort(layer_nums);
sorted_file_names = file_names(sort_idx); % Sorted list of file names

fprintf('Successfully read and sorted %d image files.\n', length(sorted_file_names));

% 4. Clean max_layer_list (remove duplicates and exclude NaN)
valid_indices = max_layer_list(~isnan(max_layer_list));
unique_indices = unique(valid_indices);

fprintf('Preparing to copy %d selected images...\n', length(unique_indices));
copy_count = 0;

% 5. Select and copy files using the sorted indices
for i = 1:length(unique_indices)
    % target_idx is the index of the Nth image to select
    target_idx = unique_indices(i);
    
    % Check that the index is positive and does not exceed the number of available images
    if target_idx > 0 && target_idx <= length(sorted_file_names)
        
        % Get the original name of the file at target_idx in the sorted list
        src_file_name = sorted_file_names{target_idx};
        
        % Construct the full source file path
        src_file = fullfile(src_dir, src_file_name);
        
        % --- Construct the destination file path ---
        % Option A: Copy using the original file name
        % dst_file_name = src_file_name; 
        
        % Option B: Rename to Layer_XXX.jpg when copying (as previously requested)
        dst_file_name = sprintf('Layer_%d.jpg', target_idx);
        % --------------------------
        
        dst_file = fullfile(dst_dir, dst_file_name);
        
        % Copy the file
        copyfile(src_file, dst_file);
        copy_count = copy_count + 1;
        
    else
        warning('Index %d is out of range (%d images available); skipping.', target_idx, length(sorted_file_names));
    end
end

fprintf('Copying completed! Successfully processed %d files.\n', copy_count);
%%
hg_to_cancer = [2 11 14 17 19 20 22 27 30];
hg_to_lg = [21 18 29];
hg_to_duct = [9 23];
hg_to_acinar = [3 12 13 15 16 24 28 32];
hg_to_stroma = [33 7 34 10 34];
hg_to_nerve = [25];
%% Reclassify HG PanIN as Cancer
for idx = hg_to_cancer
    % 1. Get the coordinate matrix (N x 3) at the current index
    coords = hg_clusters{idx}; 
    
    if ~isempty(coords)
        % 2. Extract the x, y, z coordinates
        % Note: MATLAB matrix indexing typically uses VolS(row, col, page)
        % Assume the coordinate columns represent rows, columns, and layers, respectively:
        r = coords(:, 1);
        c = coords(:, 2);
        z = coords(:, 3);
        
        % 3. Convert 3D coordinates to linear indices
        % The first argument of sub2ind is the size of the target volume
        linIdx = sub2ind(size(volS_clean), r, c, z);
        
        % 4. Set the corresponding positions in VolS to 5
        volS_clean(linIdx) = 5;
    end
end
%% Reclassify HG PanIN as LG PanIN
for idx = hg_to_lg
    % 1. Get the coordinate matrix (N x 3) at the current index
    coords = hg_clusters{idx}; 
    
    if ~isempty(coords)
        % 2. Extract the x, y, z coordinates
        % Note: MATLAB matrix indexing typically uses VolS(row, col, page)
        % Assume the coordinate columns represent rows, columns, and layers, respectively:
        r = coords(:, 1);
        c = coords(:, 2);
        z = coords(:, 3);
        
        % 3. Convert 3D coordinates to linear indices
        % The first argument of sub2ind is the size of the target volume
        linIdx = sub2ind(size(volS_clean), r, c, z);
        
        % 4. Set the corresponding positions in VolS to 1
        volS_clean(linIdx) = 3;
    end
end
%% Reclassify HG PanIN as Acinar
for idx = hg_to_acinar
    % 1. Get the coordinate matrix (N x 3) at the current index
    coords = hg_clusters{idx}; 
    
    if ~isempty(coords)
        % 2. Extract the x, y, z coordinates
        % Note: MATLAB matrix indexing typically uses VolS(row, col, page)
        % Assume the coordinate columns represent rows, columns, and layers, respectively:
        r = coords(:, 1);
        c = coords(:, 2);
        z = coords(:, 3);
        
        % 3. Convert 3D coordinates to linear indices
        % The first argument of sub2ind is the size of the target volume
        linIdx = sub2ind(size(volS_clean), r, c, z);
        
        % 4. Set the corresponding positions in VolS to 1
        volS_clean(linIdx) = 1;
    end
end

%% Reclassify HG PanIN as Duct
for idx = hg_to_duct
    % 1. Get the coordinate matrix (N x 3) at the current index
    coords = hg_clusters{idx}; 
    
    if ~isempty(coords)
        % 2. Extract the x, y, z coordinates
        % Note: MATLAB matrix indexing typically uses VolS(row, col, page)
        % Assume the coordinate columns represent rows, columns, and layers, respectively:
        r = coords(:, 1);
        c = coords(:, 2);
        z = coords(:, 3);
        
        % 3. Convert 3D coordinates to linear indices
        % The first argument of sub2ind is the size of the target volume
        linIdx = sub2ind(size(volS_clean), r, c, z);
        
        % 4. Set the corresponding positions in VolS to 1
        volS_clean(linIdx) = 2;
    end
end
%% Reclassify HG PanIN as Stroma
for idx = hg_to_stroma
    % 1. Get the coordinate matrix (N x 3) at the current index
    coords = hg_clusters{idx}; 
    
    if ~isempty(coords)
        % 2. Extract the x, y, z coordinates
        % Note: MATLAB matrix indexing typically uses VolS(row, col, page)
        % Assume the coordinate columns represent rows, columns, and layers, respectively:
        r = coords(:, 1);
        c = coords(:, 2);
        z = coords(:, 3);
        
        % 3. Convert 3D coordinates to linear indices
        % The first argument of sub2ind is the size of the target volume
        linIdx = sub2ind(size(volS_clean), r, c, z);
        
        % 4. Set the corresponding positions in VolS to 1
        volS_clean(linIdx) = 7;
    end
end
%% Reclassify HG PanIN as Nerve
for idx = hg_to_nerve
    % 1. Get the coordinate matrix (N x 3) at the current index
    coords = hg_clusters{idx}; 
    
    if ~isempty(coords)
        % 2. Extract the x, y, z coordinates
        % Note: MATLAB matrix indexing typically uses VolS(row, col, page)
        % Assume the coordinate columns represent rows, columns, and layers, respectively:
        r = coords(:, 1);
        c = coords(:, 2);
        z = coords(:, 3);
        
        % 3. Convert 3D coordinates to linear indices
        % The first argument of sub2ind is the size of the target volume
        linIdx = sub2ind(size(volS_clean), r, c, z);
        
        % 4. Set the corresponding positions in VolS to 1
        volS_clean(linIdx) = 9;
    end
end
%%
volumeViewer(volS,VolumeType="labels")
%%
volumeViewer(volS_clean,VolumeType="labels")

%% Save the final 3D model after manual checking
save('F:\P1_3D\10x\classification_05_07_2025\volumetric_data\P1_volume_smoothed_cleaned_manual_checked_final.mat', 'volS_clean', '-v7.3');

