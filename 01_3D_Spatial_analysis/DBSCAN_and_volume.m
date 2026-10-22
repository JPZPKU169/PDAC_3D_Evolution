%%
function filled_volumes = calculate_filled_volumes(cluster_cell, dx, dy, dz)
% CALCULATE_FILLED_VOLUMES Calculate the true physical volume of clustered
% 3D structures after filling the luminal spaces
%
% Input parameters:
%   cluster_cell  : Cell array, where each element is an Nx3 matrix
%                   containing the [X, Y, Z] coordinates of one cluster
%   dx            : Physical length represented by one voxel in the X direction
%                   (e.g., um or mm)
%   dy            : Physical length represented by one voxel in the Y direction
%                   (e.g., um or mm)
%   dz            : Physical length represented by one voxel in the Z direction
%                   (e.g., um or mm)
%
% Output parameter:
%   filled_volumes: 1xN vector containing the physical volume of each cluster
%                   after slice-by-slice filling of the luminal spaces

voxel_vol = dx * dy * dz;

num_clusters = length(cluster_cell);
filled_volumes = zeros(1, num_clusters);

fprintf('Starting calculation of filled volumes for %d clusters...\n', num_clusters);

for k = 1:num_clusters
    pts = round(cluster_cell{k}); 
    
    if size(pts, 1) < 3
        filled_volumes(k) = size(pts, 1) * voxel_vol;
        continue;
    end
    
    min_coords = min(pts, [], 1);
    max_coords = max(pts, [], 1);
    
    grid_size = max_coords - min_coords + 1;
    
    local_grid = false(grid_size(1), grid_size(2), grid_size(3));
    
    local_pts = pts - min_coords + 1;
    
    linear_idx = sub2ind(grid_size, local_pts(:,1), local_pts(:,2), local_pts(:,3));
    local_grid(linear_idx) = true;
    
    filled_grid = local_grid;
    for z = 1:grid_size(3)
        filled_grid(:,:,z) = imfill(local_grid(:,:,z), 'holes');
    end
    
    filled_voxel_count = sum(filled_grid(:));
    filled_volumes(k) = filled_voxel_count * voxel_vol;
    
    if mod(k, max(1, floor(num_clusters/10))) == 0 || k == num_clusters
        fprintf('  processed: %d / %d (%.0f%%)\n', k, num_clusters, (k/num_clusters)*100);
    end
end

fprintf('All clusters completed！\n');
end

%%%%%%%%%%%%%%%%%%%%%% P1 %%%%%%%%%%%%%%%%%
%% Load
load('F:\P1_3D\10x\classification_05_07_2025\volumetric_data\P1_volume_smoothed_cleaned_manual_checked_final.mat', 'volS_clean')

%%
len_per_pixel = 0.0232; % mm
height_per_layer = 39; % μm
radius=0.5 / len_per_pixel; % n

%% LG PanIN
[idx_x, idx_y, idx_z] = ind2sub(size(volS_clean), find(volS_clean == 3));
points = [idx_x, idx_y, idx_z];

% DBSCAN
epsilon = round(radius);
minPts = 5;
labels = dbscan(points, epsilon, minPts, 'Distance', 'euclidean');

num_clusters = max(labels);
lg_clusters = cell(1, num_clusters);
for k = 1:num_clusters
    lg_clusters{k} = points(labels == k, :);
end
%%
coords_dir = ...
    'F:\P1_3D\10x\classification_05_07_2025\volumetric_data\ROI_coordinates';
[lg_clusters, lg_keep_idx, lg_cluster_info] = ...
    filter_clusters_by_roi_overlap( ...
        lg_clusters, coords_dir, 20);
%%
lg_volumes = calculate_filled_volumes(lg_clusters, len_per_pixel*1000, len_per_pixel*1000, height_per_layer);

%%
Cluster_ID = (1:length(lg_volumes))'; 
Volume_Physical = lg_volumes';
vol_table = table(Cluster_ID, Volume_Physical);
writetable(vol_table, "E:\Projects\PDAC\results\3D\PanIN_volume\P1_LG_volume.csv");
save("E:\Projects\PDAC\results\3D\PanIN_volume\P1_LG_DBSCAN_clusters.mat", 'lg_clusters')

%% HG PanIN
[idx_x, idx_y, idx_z] = ind2sub(size(volS_clean), find(volS_clean == 4));
points = [idx_x, idx_y, idx_z];

% DBSCAN
epsilon = round(radius);
minPts = 5;
labels = dbscan(points, epsilon, minPts, 'Distance', 'euclidean');

num_clusters = max(labels);
hg_clusters = cell(1, num_clusters);
for k = 1:num_clusters
    hg_clusters{k} = points(labels == k, :);
end

%%
[hg_clusters, hg_keep_idx, hg_cluster_info] = ...
    filter_clusters_by_roi_overlap( ...
        hg_clusters, coords_dir, 20);

%%
hg_volumes = calculate_filled_volumes(hg_clusters, len_per_pixel*1000, len_per_pixel*1000, height_per_layer);

%%
Cluster_ID = (1:length(hg_volumes))'; 
Volume_Physical = hg_volumes';
vol_table = table(Cluster_ID, Volume_Physical);
writetable(vol_table, "E:\Projects\PDAC\results\3D\PanIN_volume\P1_HG_volume.csv");
save("E:\Projects\PDAC\results\3D\PanIN_volume\P1_HG_DBSCAN_clusters.mat", 'hg_clusters')
