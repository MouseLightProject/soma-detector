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

    
    
%%
% Load the ground-truth (these are all soma locations)
xyz_from_raw_target_index = load_soma_targets() ;

% Filter out the ones outside the substack
bounding_box_lower_corner_xyz = substack_origin_xyz - spacing_at_zoom_level_xyz/2 ;
bounding_box_upper_corner_xyz = bounding_box_lower_corner_xyz + substack_shape_xyz ;
is_within_bounding_box = ...
    all(bounding_box_lower_corner_xyz <= xyz_from_raw_target_index & xyz_from_raw_target_index < bounding_box_upper_corner_xyz, 2)  ;
xyz_from_target_index = xyz_from_raw_target_index(is_within_bounding_box, :) ;
target_count = size(xyz_from_target_index, 1) 
    

padded_substack_yxz = ...
    get_mouselight_rendered_substack(rendered_folder_path, gfp_channel_index, padded_substack_origin_ijk1, padded_substack_shape_ijk, zoom_level) ;  

substack_mip = max(padded_substack_yxz,[], 3) ;

intensity_threshold = 0.8 * 2^16 ;
minimum_volume = 500 ;  % um^3
maximum_volume = 15000 ;  % um^3
maximum_sqrt_condition_number = 10 ;
parameters = struct('intensity_threshold', {intensity_threshold}, ...
                    'minimum_volume', {minimum_volume}, ...
                    'maximum_volume', {maximum_volume}, ...
                    'maximum_sqrt_condition_number', maximum_sqrt_condition_number) ;

[voxel_count_from_candidate_index, ...
 xyz_from_candidate_index, ...
 sqrt_condition_number_from_candidate_index, ...
 max_intensity_from_candidate_index] = ... 
    compute_somalike_features_in_uint16_stack(...
        padded_substack_yxz, ...
        padded_substack_origin_xyz, ...
        spacing_at_zoom_level_xyz, ...
        substack_origin_xyz, ...
        substack_shape_xyz, ...
        parameters) ;
        % these are in the coordinate system of the full stack        
candidate_count = length(voxel_count_from_candidate_index) ;



% Calculate distance between each target and each candidate
target_candidate_distance_matrix = zeros(target_count, candidate_count) ;
for i = 1 : target_count ,
    target_i_xyz = xyz_from_target_index(i,:) ;
    for j = 1 : candidate_count ,
        candidate_j_xyz = xyz_from_candidate_index(j,:) ;
        distance = sqrt(sum((candidate_j_xyz - target_i_xyz).^2)) ;
        target_candidate_distance_matrix(i,j) = distance ;        
    end
end

is_match_distance_threshold = 10 ;
[is_target_candidate_match, target_candidate_match_distance] = match_targets_and_guesses(target_candidate_distance_matrix, is_match_distance_threshold) ;

is_there_a_matched_candidate_from_target_index = any(is_target_candidate_match, 2) 
is_there_a_matched_target_from_candidate_index = any(is_target_candidate_match, 1)'  % want a col 
[distance_to_matched_candidate_from_target_index, matching_candidate_index_from_target_index] = min(target_candidate_match_distance, [], 2) 
[distance_to_matched_target_from_candidate_index_as_row, matching_target_index_from_candidate_index_as_row] = min(target_candidate_match_distance, [], 1) ;
distance_to_matched_target_from_candidate_index = distance_to_matched_target_from_candidate_index_as_row' 
matching_target_index_from_candidate_index = matching_target_index_from_candidate_index_as_row'


candidate_hit_count = sum(is_there_a_matched_candidate_from_target_index)
assert( sum(is_there_a_matched_target_from_candidate_index) == candidate_hit_count ) ;
candidate_miss_count = sum(~is_there_a_matched_candidate_from_target_index)
candidate_chase_count = sum(~is_there_a_matched_target_from_candidate_index)

candidate_precision = candidate_hit_count / candidate_count 
candidate_recall = candidate_hit_count / target_count 


% Plot the MIP image
f = figure('color', 'w', 'name', 'targets-and-candidates') ;
a = axes(f, 'YDir', 'reverse') ;
image(a, 'CData', substack_mip, ...
         'XData', [padded_substack_origin_xyz(1) padded_substack_far_corner_xyz(1)], ...
         'YData', [padded_substack_origin_xyz(2) padded_substack_far_corner_xyz(2)], ...
         'CDataMapping', 'scaled') ;         
colormap(gray(256)) ;
axis image
xlabel('x (um)') ;
ylabel('y (um)') ;

% plot each target, and each candidate, and draw a line between matches
hold on ;
for target_index = 1 : target_count ,
    target_xyz = xyz_from_target_index(target_index,:) ;
    marker_color = fif(is_there_a_matched_candidate_from_target_index(target_index), [0 0.5 1], [1 0 0]) ;    
    plot(target_xyz(1), target_xyz(2), 'Marker', '+', 'Color', marker_color) ;            
    %text(target_xyz(1)+5, target_xyz(2)+5, sprintf('t%d', target_index), 'Color', 0.5*[1 1 1]) ;            
end
for candidate_index = 1 : candidate_count ,
    candidate_xyz = xyz_from_candidate_index(candidate_index,:) ;
    marker_color = fif(is_there_a_matched_target_from_candidate_index(candidate_index), [0 0.5 1], [1 0 0]) ;    
    plot(candidate_xyz(1), candidate_xyz(2), 'Marker', 'o', 'MarkerSize', 6, 'Color', marker_color) ;
    %text(candidate_xyz(1)-5, candidate_xyz(2)-5, sprintf('g%d', candidate_index), 'Color', 0.5*[1 1 1]) ;            
end
for target_index = 1 : target_count ,
    target_xyz = xyz_from_target_index(target_index,:) ;
    if is_there_a_matched_candidate_from_target_index(target_index) ,
        candidate_index = matching_candidate_index_from_target_index(target_index) ;
        candidate_xyz = xyz_from_candidate_index(candidate_index,:) ;
        plot([target_xyz(1) candidate_xyz(1)], [target_xyz(2) candidate_xyz(2)], 'Color', [0 0.5 1]) ;    
    end
end
hold off ; 











volume_per_voxel = prod(spacing_at_zoom_level_xyz) ;
minimum_volume_in_voxels = minimum_volume / volume_per_voxel
maximum_volume_in_voxels = maximum_volume / volume_per_voxel

is_guess_from_candidate_index = ...
    (minimum_volume_in_voxels < voxel_count_from_candidate_index) & ...
    (voxel_count_from_candidate_index < maximum_volume_in_voxels) & ...
    sqrt_condition_number_from_candidate_index < maximum_sqrt_condition_number ;
voxel_count_from_guess_index = voxel_count_from_candidate_index(is_guess_from_candidate_index) ;
xyz_from_guess_index = xyz_from_candidate_index(is_guess_from_candidate_index, :) ;    
sqrt_condition_number_from_guess_index =sqrt_condition_number_from_candidate_index(is_guess_from_candidate_index) ;
max_intensity_from_guess_index = max_intensity_from_candidate_index(is_guess_from_candidate_index) ;
guess_count = length(voxel_count_from_guess_index) ;

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
% candidate_count = size(candidate_centroid_xyz_from_label,1) ;
% for i = 1 : candidate_count ,
%     centroid_xyz = candidate_centroid_xyz_from_label(i,:) ;
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

% miss_1_xy = [74096.2049 18331.3971]  % location of a candidate centroid that is *not* classified as a soma, but should be
% distance_to_miss_1_xy = sqrt(sum((candidate_centroid_xyz_from_label(:,1:2) - miss_1_xy).^2,2)) ;
% [distance_from_miss_1_to_nearest_candidate_in_xy, label_of_miss_1] = min(distance_to_miss_1_xy)
% 
% centroid_xyz_of_miss_1 = candidate_centroid_xyz_from_label(label_of_miss_1, :)
% voxel_count_of_miss_1 = voxel_count_from_label(label_of_miss_1) 
% sqrt_condition_number_of_miss_1 = sqrt_condition_number_from_label(label_of_miss_1)
% is_putative_soma_for_miss_1 = is_putative_soma_from_label(label_of_miss_1)
% max_intesity_for_miss_1 = max_intensity_from_label(label_of_miss_1)
% % sqrt condition number is 2.6, voxel_count is 556, which just clears the
% % current minimum (515)




% Calculate distance between each target and each guess
target_guess_distance_matrix = zeros(target_count, guess_count) ;
for i = 1 : target_count ,
    target_i_xyz = xyz_from_target_index(i,:) ;
    for j = 1 : guess_count ,
        guess_j_xyz = xyz_from_guess_index(j,:) ;
        distance = sqrt(sum((guess_j_xyz - target_i_xyz).^2)) ;
        target_guess_distance_matrix(i,j) = distance ;        
    end
end

[is_target_guess_match, target_guess_match_distance] = match_targets_and_guesses(target_guess_distance_matrix, is_match_distance_threshold) ;

is_there_a_matched_guess_from_target_index = any(is_target_guess_match, 2) 
is_there_a_matched_target_from_guess_index = any(is_target_guess_match, 1) 
[distance_to_matched_guess_from_target_index, matching_guess_index_from_target_index] = min(target_guess_match_distance, [], 2) 
[distance_to_matched_target_from_guess_index, matching_target_index_from_guess_index] = min(target_guess_match_distance, [], 1) 

is_hit_from_target_index = is_there_a_matched_guess_from_target_index ;
is_hit_from_guess_index = is_there_a_matched_target_from_guess_index ;
is_miss_from_target_index = ~is_there_a_matched_guess_from_target_index ;
is_chase_from_guess_index = ~is_there_a_matched_target_from_guess_index ;

guess_hit_count = sum(is_hit_from_target_index)
assert( sum(is_hit_from_guess_index) == guess_hit_count ) ;
guess_miss_count = sum(is_miss_from_target_index)
guess_chase_count = sum(is_chase_from_guess_index)

guess_precision = guess_hit_count / guess_count 
guess_recall = guess_hit_count / target_count 


% Plot the MIP image
f = figure('color', 'w', 'name', 'targets-and-guesses') ;
a = axes(f, 'YDir', 'reverse') ;
image(a, 'CData', substack_mip, ...
         'XData', [padded_substack_origin_xyz(1) padded_substack_far_corner_xyz(1)], ...
         'YData', [padded_substack_origin_xyz(2) padded_substack_far_corner_xyz(2)], ...
         'CDataMapping', 'scaled') ;         
colormap(gray(256)) ;
axis image
xlabel('x (um)') ;
ylabel('y (um)') ;

% plot each target, and each guess, and draw a line between matches
hold on ;
for target_index = 1 : target_count ,
    target_xyz = xyz_from_target_index(target_index,:) ;
    marker_color = fif(is_there_a_matched_guess_from_target_index(target_index), [0 0.5 1], [1 0 0]) ;    
    plot(target_xyz(1), target_xyz(2), 'Marker', '+', 'Color', marker_color) ;            
    %text(target_xyz(1)+5, target_xyz(2)+5, sprintf('t%d', target_index), 'Color', 0.5*[1 1 1]) ;            
end
for guess_index = 1 : guess_count ,
    guess_xyz = xyz_from_guess_index(guess_index,:) ;
    marker_color = fif(is_there_a_matched_target_from_guess_index(guess_index), [0 0.5 1], [1 0 0]) ;    
    plot(guess_xyz(1), guess_xyz(2), 'Marker', 'o', 'MarkerSize', 6, 'Color', marker_color) ;
    %text(guess_xyz(1)-5, guess_xyz(2)-5, sprintf('g%d', guess_index), 'Color', 0.5*[1 1 1]) ;            
end
for target_index = 1 : target_count ,
    target_xyz = xyz_from_target_index(target_index,:) ;
    if is_there_a_matched_guess_from_target_index(target_index) ,
        guess_index = matching_guess_index_from_target_index(target_index) ;
        guess_xyz = xyz_from_guess_index(guess_index,:) ;
        plot([target_xyz(1) guess_xyz(1)], [target_xyz(2) guess_xyz(2)], 'Color', [0 0.5 1]) ;    
    end
end
hold off ; 





% Considering just the cadidates, plot them in feature space, and
% characterize each as hit/miss/chase/ball.  This makes sense b/c we cast a wide net: every target is
% also a candidate.

is_positive_from_candidate_index =  is_there_a_matched_target_from_candidate_index ;
is_negative_from_candidate_index = ~is_there_a_matched_target_from_candidate_index ;
positive_count_within_candidates = sum(is_positive_from_candidate_index) 
negative_count_within_candidates = sum(is_negative_from_candidate_index) 

is_test_positive_from_candidate_index =  is_guess_from_candidate_index ;
is_test_negative_from_candidate_index = ~is_guess_from_candidate_index ;
test_positive_count_within_candidates = sum(is_test_positive_from_candidate_index) 
test_negative_count_within_candidates = sum(is_test_negative_from_candidate_index) 

is_hit_from_candidate_index = is_test_positive_from_candidate_index & is_positive_from_candidate_index ;
is_miss_from_candidate_index = is_test_negative_from_candidate_index & is_positive_from_candidate_index ;
is_chase_from_candidate_index = is_test_positive_from_candidate_index & is_negative_from_candidate_index ;
is_ball_from_candidate_index = is_test_negative_from_candidate_index & is_negative_from_candidate_index ;

hit_count_within_candidates = sum(is_hit_from_candidate_index)
miss_count_within_candidates = sum(is_miss_from_candidate_index)
chase_count_within_candidates = sum(is_chase_from_candidate_index)
ball_count_within_candidates = sum(is_ball_from_candidate_index) 

precision_within_candidates = hit_count_within_candidates / test_positive_count_within_candidates 
recall_within_candidates = hit_count_within_candidates / positive_count_within_candidates 

hit_rate_within_candidates = hit_count_within_candidates / positive_count_within_candidates
ball_rate_within_candidates = ball_count_within_candidates / negative_count_within_candidates



% Plot the targets in the (very wimpy) feature space
f = figure('Color', 'w', 'Name', 'features') ;
a = axes(f) ;
h_hit = plot3(max_intensity_from_candidate_index(is_hit_from_candidate_index), ...
              voxel_count_from_candidate_index(is_hit_from_candidate_index), ...
              sqrt_condition_number_from_candidate_index(is_hit_from_candidate_index), ...
              'Marker', '.', 'Color', [0 0.3 1], 'LineStyle', 'none') ;
hold on ;
h_miss = plot3(max_intensity_from_candidate_index(is_miss_from_candidate_index), ...
               voxel_count_from_candidate_index(is_miss_from_candidate_index), ...
               sqrt_condition_number_from_candidate_index(is_miss_from_candidate_index), ...
              'Marker', 'd', 'Color', [1 0 0], 'LineStyle', 'none') ;
h_chase = plot3(max_intensity_from_candidate_index(is_chase_from_candidate_index), ...
                voxel_count_from_candidate_index(is_chase_from_candidate_index), ...
                sqrt_condition_number_from_candidate_index(is_chase_from_candidate_index), ...
              'Marker', 'd', 'Color', [138 43 226]/255, 'LineStyle', 'none') ;
h_ball = plot3(max_intensity_from_candidate_index(is_ball_from_candidate_index), ...
               voxel_count_from_candidate_index(is_ball_from_candidate_index), ...
               sqrt_condition_number_from_candidate_index(is_ball_from_candidate_index), ...
              'Marker', '.', 'Color', [0 0.8 0], 'LineStyle', 'none') ;
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


