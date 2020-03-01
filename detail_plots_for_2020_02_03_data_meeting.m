sample_date = '2019-10-04' ;
rendered_folder_path = sprintf('/nrs/mouselight/SAMPLES/%s', sample_date) ;
roi_center_xyz = [73846.8, 18235.4, 36316.3] ;  % um
roi_corner_xyz = [74884.4, 19198.9, 36316.3] ; % um 
pad_depth_in_um = 50 ; % um
zoom_level = 4 ;
gfp_channel_index = 0 ;
background_channel_index = 1 ;

intensity_threshold = 0.75 * 2^16 ;
minimum_volume = 500 ;  % um^3
%maximum_volume = 15000 ;  % um^3
maximum_volume = 25000 ;  % um^3
maximum_sqrt_condition_number = 20 ;

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

stack_shape_ijk = chunk_shape_ijk * 2^zoom_level ;  % the shape of the full-brain stack, at zoom level zoom_level
stack_shape_xyz = stack_shape_ijk .* spacing_at_zoom_level_xyz  % the shape of the full-brain stack in um, does not change with zoom level

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
padded_substack_heckbert_origin_xyz = padded_substack_origin_xyz - spacing_at_zoom_level_xyz/2 
    % coord of exterior corner of lowest-indices voxel
padded_substack_far_corner_xyz = padded_substack_origin_xyz + (padded_substack_shape_ijk-1) .* spacing_at_zoom_level_xyz ;
    % coord of voxel center of highest-indices voxel

    
    
%%
% Load the ground-truth (these are all soma locations)
%xyz_from_raw_target_index = load_soma_targets() ;
[xyz_from_raw_target_index, name_from_raw_target_index] = load_traceable_soma_targets_from_tracers() ;

% Filter out the ones outside the substack
bounding_box_lower_corner_xyz = substack_origin_xyz - spacing_at_zoom_level_xyz/2 ;
bounding_box_upper_corner_xyz = bounding_box_lower_corner_xyz + substack_shape_xyz ;
is_within_bounding_box = ...
    all(bounding_box_lower_corner_xyz <= xyz_from_raw_target_index & xyz_from_raw_target_index < bounding_box_upper_corner_xyz, 2)  ;
xyz_from_target_index = xyz_from_raw_target_index(is_within_bounding_box, :) ;
name_from_target_index = name_from_raw_target_index(is_within_bounding_box) ;
target_count = size(xyz_from_target_index, 1) 
    

padded_substack_yxz = ...
    get_mouselight_rendered_substack(rendered_folder_path, gfp_channel_index, padded_substack_origin_ijk1, padded_substack_shape_ijk, zoom_level) ;  
padded_background_substack_yxz = ...
    get_mouselight_rendered_substack(rendered_folder_path, background_channel_index, padded_substack_origin_ijk1, padded_substack_shape_ijk, zoom_level) ;  

substack_mip = max(padded_substack_yxz,[], 3) ;

parameters = struct('intensity_threshold', {intensity_threshold}, ...
                    'minimum_volume', {minimum_volume}, ...
                    'maximum_volume', {maximum_volume}, ...
                    'maximum_sqrt_condition_number', maximum_sqrt_condition_number) ;

feature_struct_from_candidate_index = ... 
    find_candidate_somata_in_uint16_stack(...
        padded_substack_yxz, ...
        padded_background_substack_yxz, ...
        padded_substack_origin_xyz, ...
        origin_at_zoom_level_xyz, ...
        spacing_at_zoom_level_xyz, ...
        substack_origin_xyz, ...
        substack_shape_xyz, ...
        parameters) ;
        % these are in the coordinate system of the full stack
candidate_count = length(feature_struct_from_candidate_index) ;
xyz_from_candidate_index = reshape([feature_struct_from_candidate_index.centroidoid_xyz], [3 candidate_count])' ;
voxel_count_from_candidate_index = [feature_struct_from_candidate_index.voxel_count]' ;
sqrt_condition_number_from_candidate_index = [feature_struct_from_candidate_index.sqrt_condition_number]' ;
max_intensity_from_candidate_index = [feature_struct_from_candidate_index.max_intensity]' ;
max_background_intensity_from_candidate_index = [feature_struct_from_candidate_index.max_background_intensity]' ;

% % Report the perf of the guesses
% print_performace_statistics_and_plot('candidates', ...
%                                      xyz_from_target_index, ...
%                                      feature_struct_from_candidate_index, ...
%                                      heckbert_origin_xyz, ...
%                                      stack_shape_xyz, ...
%                                      spacing_at_zoom_level_xyz, ...
%                                      substack_mip, ...
%                                      padded_substack_origin_xyz(1:2), ...
%                                      spacing_at_zoom_level_xyz(1:2)) ;
%                                      
% 
% 
% 




volume_per_voxel = prod(spacing_at_zoom_level_xyz) ;
minimum_volume_in_voxels = minimum_volume / volume_per_voxel  %#ok<NOPTS>
maximum_volume_in_voxels = maximum_volume / volume_per_voxel  %#ok<NOPTS>

is_guess_from_candidate_index = ...
    (minimum_volume_in_voxels < voxel_count_from_candidate_index) & ...
    (voxel_count_from_candidate_index < maximum_volume_in_voxels) & ...
    sqrt_condition_number_from_candidate_index < maximum_sqrt_condition_number & ...
    max_intensity_from_candidate_index > max_background_intensity_from_candidate_index ;
feature_struct_from_guess_index = feature_struct_from_candidate_index(is_guess_from_candidate_index) ;
guess_count = length(feature_struct_from_guess_index) ;





% Report the perf of the guesses
% print_performace_statistics_and_plot('guesses', ...
%                                      xyz_from_target_index, ...
%                                      feature_struct_from_guess_index, ...
%                                      heckbert_origin_xyz, ...
%                                      stack_shape_xyz, ...
%                                      spacing_at_zoom_level_xyz, ...
%                                      substack_mip, ...
%                                      padded_substack_origin_xyz(1:2), ...
%                                      spacing_at_zoom_level_xyz(1:2) ) ;

is_there_a_mip = true ;
mip_origin_xy = padded_substack_origin_xyz(1:2) ;
mip_spacing_xy = spacing_at_zoom_level_xyz(1:2) ;
do_plot_guesses = true ;
mip_clim = [min(substack_mip(:)) max(substack_mip(:))] ;

% Figure out what targets are in what guesses
%target_count = size(xyz_from_target_index, 1) ;
%guess_count = length(feature_struct_from_guess_index) ;
origin_at_zoom_level_xyz = heckbert_origin_xyz + spacing_at_zoom_level_xyz/2 ;
matching_guess_index_from_target_index = ...
    match_targets_and_components(xyz_from_target_index, feature_struct_from_guess_index, origin_at_zoom_level_xyz, spacing_at_zoom_level_xyz) ;
matching_target_index_from_guess_index = invert_partial_map_array(matching_guess_index_from_target_index, guess_count) ;
is_there_a_matched_guess_from_target_index = isfinite(matching_guess_index_from_target_index) ;
is_there_a_matched_target_from_guess_index = isfinite(matching_target_index_from_guess_index) ;

%
% End of new way
%

% Compute hit counts, precision, recall for the guesses
fprintf('\n\n%s:\n', 'guesses') ;
target_count  %#ok<NOPTS>
guess_count  %#ok<NOPTS>
hit_count = sum(is_there_a_matched_guess_from_target_index)  %#ok<NOPTS>
assert( sum(is_there_a_matched_target_from_guess_index) == hit_count ) ;
miss_count = sum(~is_there_a_matched_guess_from_target_index)  %#ok<NOPTS>
chase_count = sum(~is_there_a_matched_target_from_guess_index)  %#ok<NOPTS>

recall = hit_count / target_count   %#ok<NOPTS>
precision = hit_count / guess_count   %#ok<NOPTS>


% Plot the MIP image by itself
f = figure('color', 'k', 'name', 'detail-mip') ;
a = axes(f, 'YDir', 'reverse', 'XColor', 'w', 'YColor', 'w', 'Layer', 'top', 'Box', 'on', 'DataAspectRatio', [1 1 1]) ;
mip_shape_ji = size(substack_mip) ;
mip_shape_ij = mip_shape_ji([2 1]) ;
mip_shape_xy = mip_spacing_xy .* mip_shape_ij ;
mip_near_voxel_center_xy = mip_spacing_xy/2 ;
mip_far_voxel_center_xy = mip_shape_xy - mip_spacing_xy/2 ;
%mip_far_corner_xy = mip_origin_xy + mip_spacing_xy .* (mip_shape_ij-1) ;
image(a, 'CData', substack_mip, ...
         'XData', [mip_near_voxel_center_xy(1) mip_far_voxel_center_xy(1)], ...
         'YData', [mip_near_voxel_center_xy(2) mip_far_voxel_center_xy(2)], ...
         'CDataMapping', 'scaled') ;         
xlim([0 mip_shape_xy(1)]) ;
ylim([0 mip_shape_xy(2)]) ;
a.CLim = mip_clim ;
colormap(gray(256)) ;
%set(a, 'DataAspectRatio', [1 1 1]) ;
%axis image    
xlabel('x (um)', 'Color', 'w') ;
ylabel('y (um)', 'Color', 'w') ;
set_figure_to_wysiwyg_printing(f) ;
set_figure_size([13+1/3 7.5]) ;


% Plot the MIP image
f = figure('color', 'k', 'name', 'detail-targets-and-guesses-mip') ;
a = axes(f, 'YDir', 'reverse', 'XColor', 'w', 'YColor', 'w', 'Layer', 'top', 'Box', 'on', 'DataAspectRatio', [1 1 1]) ;
mip_shape_ji = size(substack_mip) ;
mip_shape_ij = mip_shape_ji([2 1]) ;
mip_shape_xy = mip_spacing_xy .* mip_shape_ij ;
mip_near_voxel_center_xy = mip_spacing_xy/2 ;
mip_far_voxel_center_xy = mip_shape_xy - mip_spacing_xy/2 ;
%mip_far_corner_xy = mip_origin_xy + mip_spacing_xy .* (mip_shape_ij-1) ;
image(a, 'CData', substack_mip, ...
         'XData', [mip_near_voxel_center_xy(1) mip_far_voxel_center_xy(1)], ...
         'YData', [mip_near_voxel_center_xy(2) mip_far_voxel_center_xy(2)], ...
         'CDataMapping', 'scaled') ;         
xlim([0 mip_shape_xy(1)]) ;
ylim([0 mip_shape_xy(2)]) ;
a.CLim = mip_clim ;
colormap(gray(256)) ;
%set(a, 'DataAspectRatio', [1 1 1]) ;
%axis image    
xlabel('x (um)', 'Color', 'w') ;
ylabel('y (um)', 'Color', 'w') ;
set_figure_to_wysiwyg_printing(f) ;
set_figure_size([13+1/3 7.5]) ;


% plot each target, and each guess, and draw a line between matches
hold on ;
for target_index = 1 : target_count ,
    target_xyz = xyz_from_target_index(target_index,:) ;
    marker_color = fif(is_there_a_matched_guess_from_target_index(target_index), [0 0.5 1], [1 0 0]) ;    
    plot(target_xyz(1)-padded_substack_heckbert_origin_xyz(1), ...
         target_xyz(2)-padded_substack_heckbert_origin_xyz(2), ...
         'Marker', '+', 'Color', marker_color) ;
    %text(target_xyz(1)+5, target_xyz(2)+5, sprintf('t%d', target_index), 'Color', 0.5*[1 1 1]) ;            
end
if do_plot_guesses ,
    for guess_index = 1 : guess_count ,
        guess_xyz = feature_struct_from_guess_index(guess_index).centroidoid_xyz ;
        marker_color = fif(is_there_a_matched_target_from_guess_index(guess_index), [0 0.5 1], [1 0 0]) ;    
        plot(guess_xyz(1)-padded_substack_heckbert_origin_xyz(1), ...
             guess_xyz(2)-padded_substack_heckbert_origin_xyz(2), ...
             'Marker', 'o', 'MarkerSize', 6, 'Color', marker_color) ;
        %text(guess_xyz(1)-5, guess_xyz(2)-5, sprintf('g%d', guess_index), 'Color', 0.5*[1 1 1]) ;            
    end
    for target_index = 1 : target_count ,
        target_xyz = xyz_from_target_index(target_index,:) ;
        if is_there_a_matched_guess_from_target_index(target_index) ,
            guess_index = matching_guess_index_from_target_index(target_index) ;
            guess_xyz = feature_struct_from_guess_index(guess_index).centroidoid_xyz ;
            plot([target_xyz(1) guess_xyz(1)]-padded_substack_heckbert_origin_xyz(1), ...
                 [target_xyz(2) guess_xyz(2)]-padded_substack_heckbert_origin_xyz(2), ...
                 'Color', [0 0.5 1]) ;    
        end
    end
end
hold off ; 


