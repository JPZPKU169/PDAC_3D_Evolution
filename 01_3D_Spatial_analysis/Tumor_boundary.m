%%
function [vol_out, smooth_tumor_bed, tumor_surface] = build_macro_tumor_boundary(vol_in, cancer_label, boundary_label, sigma_fuse, density_threshold)
% BUILD_MACRO_TUMOR_BOUNDARY Extract the largest macroscopic interface of
% a diffuse tumor in voxel space
%
% Input parameters:
%   vol_in            : Original 3D label matrix
%   cancer_label      : Label value of tumor cells (default: 5)
%   boundary_label    : Label value assigned to the tumor boundary (default: 10)
%   sigma_fuse        : Voxel radius used for density-field fusion
%                       (recommended: 10-30)
%   density_threshold : Threshold used to binarize the density field
%                       (recommended: 0.01-0.05)
%
% Output parameters:
%   vol_out           : New 3D volume with the boundary label (e.g., 10) added
%   smooth_tumor_bed  : Logical matrix of the solid tumor bed
%                       (used to determine the enclosure proportion)
%   tumor_surface     : Logical matrix containing only the thin boundary shell
%                       (used to extract the physical point cloud)

    if nargin < 5; density_threshold = 0.02; end
    if nargin < 4; sigma_fuse = 10; end
    if nargin < 3; boundary_label = 10; end
    if nargin < 2; cancer_label = 5; end

    fprintf('Building tumor surface (fuse radius: %d, cutoff: %.3f)...\n', sigma_fuse, density_threshold);

    all_cancer = (vol_in == cancer_label);

    if ~any(all_cancer(:))
        warning('No Tumor detected!');
        vol_out = vol_in; smooth_tumor_bed = false(size(vol_in)); tumor_surface = false(size(vol_in));
        return;
    end

    tumor_density = imgaussfilt3(double(all_cancer), sigma_fuse);
    fused_tumor_bed = tumor_density > density_threshold;

    CC = bwconncomp(fused_tumor_bed, 26);
    if CC.NumObjects == 0
        vol_out = vol_in; smooth_tumor_bed = false(size(vol_in)); tumor_surface = false(size(vol_in));
        return;
    end
    
    numPixels = cellfun(@numel, CC.PixelIdxList);
    [~, max_idx] = max(numPixels);
    macro_tumor_core = false(size(vol_in));
    macro_tumor_core(CC.PixelIdxList{max_idx}) = true;

    smooth_tumor_bed = imgaussfilt3(double(macro_tumor_core), 1) > 0.5;
    tumor_surface = bwperim(smooth_tumor_bed, 26);
    tumor_surface_thick = imdilate(tumor_surface, strel('sphere', 1));

    vol_out = vol_in;
    vol_out(tumor_surface_thick) = boundary_label;
    
    fprintf('=> Tumor morphological boundary construction completed! Boundary labeled as %d\n', boundary_label);
end

%%
function result_table = analyze_panin_spatial_distribution(cluster_cell, shrunk_tumor_bed, shrunk_surface, expanded_tumor_bed, dx, dy, dz, csv_save_path)
% ANALYZE_PANIN_SPATIAL_DISTRIBUTION Calculate the distance from PanIN
% clusters to the contracted tumor boundary and their distribution
% proportions across three spatial regions
%
% Input parameters:
%   cluster_cell       : Cell array containing PanIN point-cloud coordinates
%                        (e.g., lg_clusters)
%   shrunk_tumor_bed   : Logical matrix of the contracted solid tumor bed
%                        (all voxels inside the tumor bed are 1)
%   shrunk_surface     : Logical matrix of the contracted surface boundary shell
%                        (only boundary voxels are 1)
%   expanded_tumor_bed : Logical matrix of the expanded solid tumor bed
%   dx, dy, dz         : Physical length represented by one voxel in the
%                        X, Y, and Z directions (um)
%   csv_save_path      : Path for saving the results

    num_clusters = length(cluster_cell);
    grid_size = size(shrunk_tumor_bed);
    
    disp('=> Extracting the contracted boundary point cloud and constructing the physical search space...');
    [bound_Y, bound_X, bound_Z] = ind2sub(grid_size, find(shrunk_surface));
    
    bound_coords_phys = [bound_X * dx, bound_Y * dy, bound_Z * dz];
    
    Cluster_ID = (1:num_clusters)';
    Centroid_Dist_to_Shrunk_um = zeros(num_clusters, 1);
    Min_Dist_to_Shrunk_um = zeros(num_clusters, 1);
    
    Pct_Inside_Shrunk = zeros(num_clusters, 1);
    Pct_Between_Boundaries = zeros(num_clusters, 1);
    Pct_Outside_Expanded = zeros(num_clusters, 1);
    
    fprintf('=> Starting calculation of spatial distribution and physical distances for %d PanIN clusters...\n', num_clusters);
    
    for k = 1:num_clusters
        pts = round(cluster_cell{k});
        
        if isempty(pts)
            Centroid_Dist_to_Shrunk_um(k) = NaN; Min_Dist_to_Shrunk_um(k) = NaN;
            continue;
        end
        
        panin_Y = pts(:,1); panin_X = pts(:,2); panin_Z = pts(:,3);
        total_pts = length(panin_Y);
        
        linear_idx = sub2ind(grid_size, panin_Y, panin_X, panin_Z);
        
        is_inside_shrunk = shrunk_tumor_bed(linear_idx);
        is_inside_expanded = expanded_tumor_bed(linear_idx);
        
        pts_inside_shrunk = sum(is_inside_shrunk);
        pts_between = sum(is_inside_expanded & ~is_inside_shrunk);
        pts_outside_expanded = sum(~is_inside_expanded);

        Pct_Inside_Shrunk(k) = (pts_inside_shrunk / total_pts) * 100;
        Pct_Between_Boundaries(k) = (pts_between / total_pts) * 100;
        Pct_Outside_Expanded(k) = (pts_outside_expanded / total_pts) * 100;
        
        panin_coords_phys = [panin_X * dx, panin_Y * dy, panin_Z * dz];
        
        try
            [~, dists] = knnsearch(bound_coords_phys, panin_coords_phys);
        catch
            dists = zeros(size(panin_coords_phys, 1), 1);
            for i = 1:size(panin_coords_phys, 1)
                diff_sq = (bound_coords_phys(:,1) - panin_coords_phys(i,1)).^2 + ...
                          (bound_coords_phys(:,2) - panin_coords_phys(i,2)).^2 + ...
                          (bound_coords_phys(:,3) - panin_coords_phys(i,3)).^2;
                dists(i) = sqrt(min(diff_sq));
            end
        end
        
        dists(is_inside_shrunk) = 0;
        Min_Dist_to_Shrunk_um(k) = min(dists);
        
        centroid_sub = round(mean(pts, 1));
        centroid_sub = max(centroid_sub, 1);
        centroid_sub = min(centroid_sub, [grid_size(1), grid_size(2), grid_size(3)]);
        
        if shrunk_tumor_bed(centroid_sub(1), centroid_sub(2), centroid_sub(3))
            Centroid_Dist_to_Shrunk_um(k) = 0;
        else
            centroid_phys = [centroid_sub(2)*dx, centroid_sub(1)*dy, centroid_sub(3)*dz];
            try
                [~, c_dist] = knnsearch(bound_coords_phys, centroid_phys);
            catch
                c_dist = sqrt(min((bound_coords_phys(:,1)-centroid_phys(1)).^2 + ...
                                  (bound_coords_phys(:,2)-centroid_phys(2)).^2 + ...
                                  (bound_coords_phys(:,3)-centroid_phys(3)).^2));
            end
            Centroid_Dist_to_Shrunk_um(k) = c_dist;
        end
    end
    
    result_table = table(Cluster_ID, Centroid_Dist_to_Shrunk_um, Min_Dist_to_Shrunk_um, ...
                         Pct_Inside_Shrunk, Pct_Between_Boundaries, Pct_Outside_Expanded, ...
        'VariableNames', {'Cluster_ID', 'Centroid_Dist_to_Shrunk_um', 'Shortest_Dist_to_Shrunk_um', ...
                          'Pct_Inside_Shrunk', 'Pct_Between_Shrunk_and_Expanded', 'Pct_Outside_Expanded'});
    
    if nargin >= 8 && ~isempty(csv_save_path)
        [out_dir, ~, ~] = fileparts(csv_save_path);
        if ~isempty(out_dir) && ~isfolder(out_dir), mkdir(out_dir); end
        writetable(result_table, csv_save_path, 'Encoding', 'UTF-8');
        fprintf('=> Complete! Results are saved to: %s\n', csv_save_path);
    end
end


%% P1
load('F:\P1_3D\10x\classification_05_07_2025\volumetric_data\P1_volume_smoothed_cleaned_manual_checked_final.mat', 'volS_clean');

%%
cancer_label = 5;
boundary_label = 10;
sigma_fuse = 10; 
density_threshold = 0.02;

[volS_with_boundary, tumor_bed, tumor_shell] = build_macro_tumor_boundary(...
    volS_clean, cancer_label, boundary_label, sigma_fuse, density_threshold);
%%
volumeViewer(volS_with_boundary, VolumeType="labels")

%%
dx = 23.2; % um
dy = 23.2; % um
dz = 39.0; % um

[volS_final, shrunk_bed, shrunk_shell, expanded_bed, expanded_shell] = shrink_and_expand_boundary_v2(...
    volS_clean, tumor_bed, dx, dy, dz, 5, 11, 12, 3, 500);

%%
volumeViewer(volS_final, VolumeType="labels")

%%
dx = 23.2; % um
dy = 23.2; % um
dz = 39.0; % um
csv_path = 'E:\Projects\PDAC\results\3D\PanIN_distance\P1_LG_distance.csv';

zonal_metrics_table = analyze_panin_spatial_distribution(...
    lg_clusters, shrunk_bed, shrunk_shell, expanded_bed, dx, dy, dz, csv_path);
%%
load("E:\Projects\PDAC\results\3D\PanIN_volume\P1_HG_DBSCAN_clusters.mat", 'hg_clusters')
%%
dx = 23.2; % um
dy = 23.2; % Y um
dz = 39.0; % Z um
csv_path = 'E:\Projects\PDAC\results\3D\PanIN_distance\P1_HG_distance.csv';

zonal_metrics_table = analyze_panin_spatial_distribution(...
    hg_clusters, shrunk_bed, shrunk_shell, expanded_bed, dx, dy, dz, csv_path);

