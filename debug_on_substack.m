close all 
clear

sample_date = '2019-10-04' ;
rendered_folder_path = sprintf('/nrs/mouselight/SAMPLES/%s', sample_date) ;
roi_center_xyz = [73846.8, 18235.4, 36316.3] ;  % um
roi_corner_xyz = [74884.4, 19198.9, 36316.3] ; % um 
pad_depth_in_um = 50 ; % um
zoom_level = 4 ;
gfp_channel_index = 0 ;
background_channel_index = 1 ;

render_parameters_file_path = fullfile(rendered_folder_path, 'calculated_parameters.jl') ;
render_parameters = read_renderer_calculated_parameters_file(render_parameters_file_path) ;
max_zoom_level = render_parameters.level_step_count ;
chunk_shape_ijk = render_parameters.leaf_shape ;  % xyz order, same at all zoom levels, just the chunk count changes
spacing_at_max_zoom_xyz = render_parameters.spacing ;
origin_at_max_zoom_xyz = render_parameters.origin ;

spacing_at_zoom_level_0_xyz = 2^max_zoom_level * spacing_at_max_zoom_xyz ;
spacing_at_zoom_level_xyz = spacing_at_zoom_level_0_xyz ./ (2^zoom_level) ;

heckbert_origin_xyz = origin_at_max_zoom_xyz - spacing_at_max_zoom_xyz/2 ;  % this origin does not change with the zoom level
origin_at_zoom_level_xyz = heckbert_origin_xyz + spacing_at_zoom_level_xyz/2 ;

substack_desired_half_shape_xyz = abs(roi_corner_xyz - roi_center_xyz) ;
substack_desired_half_shape_xyz(3) = mean(substack_desired_half_shape_xyz(1:2)) ;  % want similar extent in z
substack_desired_shape_xyz = 2 * substack_desired_half_shape_xyz 
substack_desired_origin_xyz = roi_center_xyz - substack_desired_half_shape_xyz 
substack_shape_ijk = round(substack_desired_shape_xyz ./ spacing_at_zoom_level_xyz) ;
substack_origin_ijk1 = round((substack_desired_origin_xyz - origin_at_zoom_level_xyz) ./ spacing_at_zoom_level_xyz) + 1 ;
   % origin of the substack within the full stack at this zoom
substack_shape_xyz = spacing_at_zoom_level_xyz .* substack_shape_ijk 
substack_origin_xyz = origin_at_zoom_level_xyz + spacing_at_zoom_level_xyz .* (substack_origin_ijk1-1) 


pad_depth_ijk = ceil(pad_depth_in_um ./ spacing_at_zoom_level_xyz) ;  % 3 x 1
padded_substack_origin_ijk1 = substack_origin_ijk1 - pad_depth_ijk 
padded_substack_shape_ijk = substack_shape_ijk + 2*pad_depth_ijk 
padded_substack_shape_xyz = spacing_at_zoom_level_xyz .* padded_substack_shape_ijk 
padded_substack_origin_xyz = origin_at_zoom_level_xyz + spacing_at_zoom_level_xyz .* (padded_substack_origin_ijk1-1) 
    % coord of voxel center of lowest-indices voxel
padded_substack_far_corner_xyz = padded_substack_origin_xyz + (padded_substack_shape_ijk-1) .* spacing_at_zoom_level_xyz ;
    % coord of voxel center of highest-indices voxel



padded_substack_yxz = ...
    get_mouselight_rendered_substack(rendered_folder_path, gfp_channel_index, padded_substack_origin_ijk1, padded_substack_shape_ijk, zoom_level) ;  

substack_mip = max(padded_substack_yxz,[], 3) ;

intensity_threshold = 0.8 * 2^16 ;
minimum_volume = 500 ;  % um^3
maximum_volume = 15000 ;  % um^3
maximum_sqrt_condition_number = 5 ;
parameters = struct('intensity_threshold', {intensity_threshold}, ...
                    'minimum_volume', {minimum_volume}, ...
                    'maximum_volume', {maximum_volume}, ...
                    'maximum_sqrt_condition_number', maximum_sqrt_condition_number) ;

[voxel_count_from_label, ...
 component_centroid_xyz_from_label, ...
 sqrt_condition_number_from_label, ...
 max_intensity_from_label] = ... 
    compute_somalike_features_in_uint16_stack(...
        padded_substack_yxz, ...
        padded_substack_origin_xyz, ...
        spacing_at_zoom_level_xyz, ...
        substack_origin_xyz, ...
        substack_shape_xyz, ...
        parameters) ;
        % these are in the coordinate system of the full stack        
component_count = length(voxel_count_from_label) ;

volume_per_voxel = prod(spacing_at_zoom_level_xyz) ;
minimum_volume_in_voxels = minimum_volume / volume_per_voxel
maximum_volume_in_voxels = maximum_volume / volume_per_voxel

is_putative_soma_from_label = ...
    (minimum_volume_in_voxels < voxel_count_from_label) & ...
    (voxel_count_from_label < maximum_volume_in_voxels) & ...
    sqrt_condition_number_from_label < maximum_sqrt_condition_number ;
putative_somata_xyzs = component_centroid_xyz_from_label(is_putative_soma_from_label, :) ;    
        
% f = figure('color', 'w') ;
% a = axes(f, 'YDir', 'reverse') ;
% image(a, 'CData', substack_mip, ...
%          'XData', [padded_substack_origin_xyz(1) padded_substack_far_corner_xyz(1)], ...
%          'YData', [padded_substack_origin_xyz(2) padded_substack_far_corner_xyz(2)], ...
%          'CDataMapping', 'scaled') ;         
% colormap(gray(256)) ;
% axis image
% 
% hold on ;
% component_count = size(component_centroid_xyz_from_label,1) ;
% for i = 1 : component_count ,
%     centroid_xyz = component_centroid_xyz_from_label(i,:) ;
%     plot(centroid_xyz(1), centroid_xyz(2), '+', 'Color', [0 0.5 1]) ;    
% end
% hold off ; 
% 
% hold on ;
% somata_count = size(putative_somata_xyzs,1) ;
% for i = 1 : somata_count ,
%     soma_xyz = putative_somata_xyzs(i,:) ;
%     plot(soma_xyz(1), soma_xyz(2), 'r+') ;    
% end
% hold off ; 

% miss_1_xy = [74096.2049 18331.3971]  % location of a component centroid that is *not* classified as a soma, but should be
% distance_to_miss_1_xy = sqrt(sum((component_centroid_xyz_from_label(:,1:2) - miss_1_xy).^2,2)) ;
% [distance_from_miss_1_to_nearest_component_in_xy, label_of_miss_1] = min(distance_to_miss_1_xy)
% 
% centroid_xyz_of_miss_1 = component_centroid_xyz_from_label(label_of_miss_1, :)
% voxel_count_of_miss_1 = voxel_count_from_label(label_of_miss_1) 
% sqrt_condition_number_of_miss_1 = sqrt_condition_number_from_label(label_of_miss_1)
% is_putative_soma_for_miss_1 = is_putative_soma_from_label(label_of_miss_1)
% max_intesity_for_miss_1 = max_intensity_from_label(label_of_miss_1)
% % sqrt condition number is 2.6, voxel_count is 556, which just clears the
% % current minimum (515)

%%
% Load the ground-truth
[raw_target_xyzs, is_raw_target_a_soma] = load_targets() ;

% Filter out the ones outside the substack
bounding_box_lower_corner_xyz = substack_origin_xyz - spacing_at_zoom_level_xyz/2 ;
bounding_box_upper_corner_xyz = bounding_box_lower_corner_xyz + substack_shape_xyz ;
is_within_bounding_box = ...
    all(bounding_box_lower_corner_xyz <= raw_target_xyzs & raw_target_xyzs < bounding_box_upper_corner_xyz, 2)  ;
xyz_from_target_index = raw_target_xyzs(is_within_bounding_box, :) ;
is_soma_from_target_index = is_raw_target_a_soma(is_within_bounding_box) ;  

% is_soma_from_target_index is all true

target_count = length(is_soma_from_target_index) ;



% Calculate distances between all pairs
distance_matrix = zeros(target_count, component_count) ;
for i = 1 : target_count ,
    target_i_xyz = xyz_from_target_index(i,:) ;
    for j = 1 : component_count ,
        component_j_xyz = component_centroid_xyz_from_label(j,:) ;
        distance = sqrt(sum((component_j_xyz - target_i_xyz).^2)) ;
        distance_matrix(i,j) = distance ;        
    end
end

[distance_to_nearest_component_from_target_index, label_from_target_index] = min(distance_matrix, [], 2) 
is_there_a_nearby_component_from_target_index = (distance_to_nearest_component_from_target_index<10) ;
is_test_positive_from_target_index = is_there_a_nearby_component_from_target_index & is_putative_soma_from_label(label_from_target_index) ;
is_positive_from_target_index = is_soma_from_target_index ;
is_negative_from_target_index = ~is_soma_from_target_index ;

positive_count = sum(is_positive_from_target_index)
negative_count = sum(is_negative_from_target_index)

is_hit = is_positive_from_target_index & is_test_positive_from_target_index ;
is_miss = is_positive_from_target_index & ~is_test_positive_from_target_index ;
is_chase = is_negative_from_target_index & is_test_positive_from_target_index ;
is_ball = is_negative_from_target_index & ~is_test_positive_from_target_index ;

hit_count = sum(is_hit)
miss_count = sum(is_miss)
chase_count = sum(is_chase)
ball_count = sum(is_ball)

hit_rate = hit_count / positive_count 
miss_rate = miss_count / positive_count
chase_rate = chase_count / negative_count
ball_rate  = ball_count / negative_count


% Get the feature values for the component matching each target
voxel_count_from_target_index = voxel_count_from_label(label_from_target_index) ;
voxel_count_from_target_index(~is_there_a_nearby_component_from_target_index) = 0 ;

sqrt_condition_number_from_target_index = sqrt_condition_number_from_label(label_from_target_index) ;
sqrt_condition_number_from_target_index(~is_there_a_nearby_component_from_target_index) = 1 ;

max_intensity_from_target_index = max_intensity_from_label(label_from_target_index) ;
max_intensity_from_target_index(~is_there_a_nearby_component_from_target_index) = 0 ;


% Plot each target, and its nearest component centroid
f = figure('color', 'w', 'name', 'targets-and-their-matches') ;
a = axes(f, 'YDir', 'reverse') ;
image(a, 'CData', substack_mip, ...
         'XData', [padded_substack_origin_xyz(1) padded_substack_far_corner_xyz(1)], ...
         'YData', [padded_substack_origin_xyz(2) padded_substack_far_corner_xyz(2)], ...
         'CDataMapping', 'scaled') ;         
colormap(gray(256)) ;
axis image

% plot each target, and its matching component if there was one
hold on ;
for target_index = 1 : target_count ,
    target_xyz = xyz_from_target_index(target_index,:) ;
    marker = fif(is_soma_from_target_index(target_index), '+', 'o') ;
    marker_color = fif(is_soma_from_target_index(target_index), [0 0.5 1], [1 0 0]) ;
    plot(target_xyz(1), target_xyz(2), 'Marker', marker, 'Color', marker_color) ;            
    if is_there_a_nearby_component_from_target_index(target_index) ,
        label = label_from_target_index(target_index) ;
        component_xyz = component_centroid_xyz_from_label(label,:) ;
        marker = fif(is_putative_soma_from_label(label), '+', 'o') ;
        marker_color = fif(is_putative_soma_from_label(label), [0 1 0], [1 2/3 0]) ;    
        plot(component_xyz(1), component_xyz(2), 'Marker', marker, 'MarkerSize', 6, 'Color', marker_color) ;
        plot([target_xyz(1) component_xyz(1)], [target_xyz(2) component_xyz(2)], 'Color', 0*[1 1 1]) ;    
    end
end
hold off ; 









% Plot the targets in the (very wimpy) feature space
f = figure('Color', 'w', 'Name', 'features') ;
a = axes(f) ;
h_hit = plot3(max_intensity_from_target_index(is_hit), voxel_count_from_target_index(is_hit), sqrt_condition_number_from_target_index(is_hit), ...
              'Marker', '+', 'Color', [0 0.5 1], 'LineStyle', 'none') ;
hold on ;
h_miss = plot3(max_intensity_from_target_index(is_miss), voxel_count_from_target_index(is_miss), sqrt_condition_number_from_target_index(is_miss), ...
              'Marker', 'o', 'Color', [1 0 0], 'LineStyle', 'none') ;
h_chase = plot3(max_intensity_from_target_index(is_chase), voxel_count_from_target_index(is_chase), sqrt_condition_number_from_target_index(is_chase), ...
              'Marker', 'o', 'Color', [1 2/3 0], 'LineStyle', 'none') ;
h_ball = plot3(max_intensity_from_target_index(is_ball), voxel_count_from_target_index(is_ball), sqrt_condition_number_from_target_index(is_ball), ...
              'Marker', '+', 'Color', [0 1 0], 'LineStyle', 'none') ;
hold off ;
xlabel('Max GFP signal (counts)') ;
ylabel('Voxel count') ;
zlabel('SD ratio') ;
%xlim([0 2^16]) ;
%ylim([0 2^16]) ;
handles = zeros(1,0) ;
legend_labels = cell(1,0) ;
if ~isempty(h_hit) ,
    handles(1, end+1) = h_hit ;
    legend_labels{1, end+1} = 'hit' ;
end
if ~isempty(h_miss) ,
    handles(1, end+1) = h_miss ;
    legend_labels{1, end+1} = 'miss' ;
end
if ~isempty(h_chase) ,
    handles(1, end+1) = h_chase ;
    legend_labels{1, end+1} = 'chase' ;
end
if ~isempty(h_ball) ,
    handles(1, end+1) = h_ball ;
    legend_labels{1, end+1} = 'ball' ;
end

legend(handles, legend_labels, 'Location', 'northwest') ;
axis vis3d
grid on
camproj perspective


