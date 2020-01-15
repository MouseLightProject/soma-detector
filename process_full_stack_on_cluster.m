do_force_computation = true ;
do_use_bsub = true ;
do_actually_submit = true ;
bsub_options_template = '-P mouselight -n8 -eo find-somata-%d-%d-%d.stdouterr.txt -oo find-somata-%d-%d-%d.stdouterr.txt -W 59 -J find-somata' ;
%bsub_options = '-P mouselight -n8 -eo /dev/null -oo /dev/null -W 59 -J find-somata' ;

sample_date = '2019-10-04' ;
rendered_folder_path = sprintf('/nrs/mouselight/SAMPLES/%s', sample_date) ;
gfp_channel_index = 0 ;
background_channel_index = 1 ;
zoom_level = 4 ;  % The zoom level of the tiles we will analyze
pad_depth_in_um = 50 ; % um
%intensity_threshold = 0.8 * 2^16 ;  % == 52428.8
intensity_threshold = 40000 ;
minimum_volume = 500 ;  % um^3
maximum_volume = 15000 ;  % um^3
maximum_sqrt_condition_number = 10 ;
parameters = struct('intensity_threshold', {intensity_threshold}, ...
                    'minimum_volume', {minimum_volume}, ...
                    'maximum_volume', {maximum_volume}, ...
                    'maximum_sqrt_condition_number', maximum_sqrt_condition_number) ;
                
this_file_path = mfilename('fullpath') ;
this_folder_path = fileparts(this_file_path) ;                
somata_folder_path = fullfile(this_folder_path, 'auto-somata') ;

                
%[shape_xyz, origin_xyz, spacing_xyz] = load_sample_shape_origin_and_spacing(rendered_folder_path) ;
render_parameters_file_path = fullfile(rendered_folder_path, 'calculated_parameters.jl') ;
render_parameters = read_renderer_calculated_parameters_file(render_parameters_file_path) ;
max_zoom_level = render_parameters.level_step_count ;
chunk_shape_ijk = render_parameters.leaf_shape ;  % xyz order, same at all zoom levels, just the chunk count changes
spacing_at_max_zoom_xyz = render_parameters.spacing ;
origin_at_max_zoom_xyz = render_parameters.origin ;
%shape_xyz = 2^max_zoom_level * chunk_shape_ijk ;

spacing_at_zoom_level_0_xyz = 2^max_zoom_level * spacing_at_max_zoom_xyz ;
spacing_at_zoom_level_xyz = spacing_at_zoom_level_0_xyz ./ (2^zoom_level) ;

heckbert_origin_xyz = origin_at_max_zoom_xyz - spacing_at_max_zoom_xyz/2 ;  % this origin does not change with the zoom level
origin_at_zoom_level_xyz = heckbert_origin_xyz + spacing_at_zoom_level_xyz/2 ;

stack_shape_ijk = chunk_shape_ijk * 2^zoom_level ;  % the shape of the full-brain stack, at zoom level zoom_level
stack_shape_xyz = stack_shape_ijk .* spacing_at_zoom_level_xyz  % the shape of the full-brain stack in um, does not change with zoom level
analysis_chunk_shape_ijk = 4*chunk_shape_ijk ;
chunks_per_dimension_ijk = stack_shape_ijk ./ analysis_chunk_shape_ijk ;
if ~exist(somata_folder_path, 'file') ,
    mkdir(somata_folder_path) ;
end
if do_force_computation ,
    system(sprintf('rm -rf %s/*', somata_folder_path)) ;
end
chunks_folder_path = fullfile(somata_folder_path, 'chunks') ;
if ~exist(chunks_folder_path, 'file') ,
    mkdir(chunks_folder_path) ;
end
job_count = prod(chunks_per_dimension_ijk) ;
job_id_array = repmat(-1, [1 job_count]) ;
job_index = 1 ;
for chunk_i = 1 : chunks_per_dimension_ijk(1) ,
    for chunk_j = 1 : chunks_per_dimension_ijk(2) ,
        for chunk_k = 1 : chunks_per_dimension_ijk(3) ,
            somata_mat_file_name = sprintf('somata-for-chunk-%d-%d-%d.mat', chunk_i, chunk_j, chunk_k) ;
            somata_mat_file_path = fullfile(chunks_folder_path, somata_mat_file_name) ;
            chunk_offset_within_chunks_ijk1 = [chunk_i chunk_j chunk_k] 
            chunk_offset_within_stack_ijk1 = (chunk_offset_within_chunks_ijk1-1) .* analysis_chunk_shape_ijk + 1 ;
            if ~exist(somata_mat_file_path, 'file') ,
                if do_use_bsub , 
                    bsub_options = sprintf(bsub_options_template, chunk_i, chunk_j, chunk_k, chunk_i, chunk_j, chunk_k) 
                    job_id_array(job_index) = ...
                        bsub(do_actually_submit, ...
                             bsub_options, ...
                             @pad_and_find_candidate_somata_then_save, ...
                                 somata_mat_file_path, ...
                                 rendered_folder_path, ...
                                 gfp_channel_index, ...
                                 zoom_level, ...
                                 spacing_at_zoom_level_xyz, ...
                                 origin_at_zoom_level_xyz, ...
                                 chunk_offset_within_stack_ijk1, ...
                                 analysis_chunk_shape_ijk, ...
                                 pad_depth_in_um, ...
                                 parameters) ;
                else
                    tic_id = tic() ;
                    pad_and_find_candidate_somata_then_save(somata_mat_file_path, ...
                                                            rendered_folder_path, ...
                                                            gfp_channel_index, ...
                                                            zoom_level, ...
                                                            spacing_at_zoom_level_xyz, ...
                                                            origin_at_zoom_level_xyz, ...
                                                            chunk_offset_within_stack_ijk1, ...
                                                            analysis_chunk_shape_ijk, ...
                                                            pad_depth_in_um, ...
                                                            parameters) ;
                    toc(tic_id) ;                          
                end
            end
            job_index = job_index + 1 ;
        end
    end
end

% Wait for all those jobs to complete
bwait(job_id_array) ;

%%
% Read in all the .mat files
soma_mat_file_name_template = fullfile(chunks_folder_path, 'somata-for-chunk-*.mat') ;
soma_mat_file_names = simple_dir(soma_mat_file_name_template) ;
somata_file_count = length(soma_mat_file_names) ;
%xyz_from_guess_index = zeros(0,3) ;
%feature_struct_from_guess_index = struct_with_shape_and_fields([0 1], somalike_feature_list()) ;
xyz_from_candidate_index = zeros(0,3) ;
feature_struct_from_candidate_index = struct_with_shape_and_fields([0 1], somalike_feature_list()) ;
for i = 1 : somata_file_count ,
    somata_mat_file_name = soma_mat_file_names{i} ;
    somata_mat_file_path = fullfile(chunks_folder_path, somata_mat_file_name) ;
    s = load(somata_mat_file_path) ;        
    %xyz_from_guess_index = vertcat(xyz_from_guess_index, s.xyz_from_guess_index) ;  %#ok<AGROW>
    %feature_struct_from_guess_index = vertcat(feature_struct_from_guess_index, s.feature_struct_from_guess_index) ;  %#ok<AGROW>
    xyz_from_candidate_index = vertcat(xyz_from_candidate_index, s.xyz_from_candidate_index) ;  %#ok<AGROW>
    feature_struct_from_candidate_index = vertcat(feature_struct_from_candidate_index, s.feature_struct_from_candidate_index) ;  %#ok<AGROW>
end
%guess_count = size(xyz_from_guess_index, 1) ;
candidate_count = size(xyz_from_candidate_index, 1) ;

% %%
% % Save as a .swc file, which can hold a forest, it turns out
% forest_name = sprintf('%s-auto-somata', sample_date) ;
% forest_color = [1 0 0] ;
% swc_file_name = horzcat(forest_name, '.swc') ;
% swc_file_path = fullfile(somata_folder_path, swc_file_name) ;
% save_somata_as_single_swc(swc_file_path, xyz_from_guess_index, forest_name, forest_color) ;


% Load the ground-truth (these are all soma locations)
xyz_from_target_index = load_soma_targets() ;
target_count = size(xyz_from_target_index, 1) ;










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

is_match_distance_threshold = 20 ;
[is_target_candidate_match, target_candidate_match_distance] = match_targets_and_guesses(target_candidate_distance_matrix, is_match_distance_threshold) ;

is_there_a_matched_candidate_from_target_index = any(is_target_candidate_match, 2) ;
is_there_a_matched_target_from_candidate_index = any(is_target_candidate_match, 1)' ; % want a col 
[distance_to_matched_candidate_from_target_index, matching_candidate_index_from_target_index] = min(target_candidate_match_distance, [], 2) ;
[distance_to_matched_target_from_candidate_index_as_row, matching_target_index_from_candidate_index_as_row] = min(target_candidate_match_distance, [], 1) ;
distance_to_matched_target_from_candidate_index = distance_to_matched_target_from_candidate_index_as_row' ;
matching_target_index_from_candidate_index = matching_target_index_from_candidate_index_as_row' ;

target_count
candidate_count
candidate_hit_count = sum(is_there_a_matched_candidate_from_target_index)
assert( sum(is_there_a_matched_target_from_candidate_index) == candidate_hit_count ) ;
candidate_miss_count = sum(~is_there_a_matched_candidate_from_target_index)
candidate_chase_count = sum(~is_there_a_matched_target_from_candidate_index)

candidate_precision = candidate_hit_count / candidate_count 
candidate_recall = candidate_hit_count / target_count 


% Plot the MIP image
f = figure('color', 'w', 'name', 'targets-and-candidates') ;
a = axes(f, 'YDir', 'reverse') ;
% image(a, 'CData', substack_mip, ...
%          'XData', [padded_substack_origin_xyz(1) padded_substack_far_corner_xyz(1)], ...
%          'YData', [padded_substack_origin_xyz(2) padded_substack_far_corner_xyz(2)], ...
%          'CDataMapping', 'scaled') ;         
colormap(gray(256)) ;
axis image
xlabel('x (um)') ;
ylabel('y (um)') ;
xlim(heckbert_origin_xyz(1)+[0 stack_shape_xyz(1)]) ;
ylim(heckbert_origin_xyz(2)+[0 stack_shape_xyz(2)]) ;

% plot each target, and each candidate, and draw a line between matches
hold on ;
for target_index = 1 : target_count ,
    target_xyz = xyz_from_target_index(target_index,:) ;
    marker_color = fif(is_there_a_matched_candidate_from_target_index(target_index), [0 0.5 1], [1 0 0]) ;    
    plot(target_xyz(1), target_xyz(2), 'Marker', '+', 'Color', marker_color) ;            
    text(target_xyz(1)+5, target_xyz(2)+5, sprintf('t%d', target_index), 'Color', 0.5*[1 1 1]) ;            
end
% for candidate_index = 1 : candidate_count ,
%     candidate_xyz = xyz_from_candidate_index(candidate_index,:) ;
%     marker_color = fif(is_there_a_matched_target_from_candidate_index(candidate_index), [0 0.5 1], [1 0 0]) ;    
%     plot(candidate_xyz(1), candidate_xyz(2), 'Marker', 'o', 'MarkerSize', 6, 'Color', marker_color) ;
%     %text(candidate_xyz(1)-5, candidate_xyz(2)-5, sprintf('g%d', candidate_index), 'Color', 0.5*[1 1 1]) ;            
% end
% for target_index = 1 : target_count ,
%     target_xyz = xyz_from_target_index(target_index,:) ;
%     if is_there_a_matched_candidate_from_target_index(target_index) ,
%         candidate_index = matching_candidate_index_from_target_index(target_index) ;
%         candidate_xyz = xyz_from_candidate_index(candidate_index,:) ;
%         plot([target_xyz(1) candidate_xyz(1)], [target_xyz(2) candidate_xyz(2)], 'Color', [0 0.5 1]) ;    
%     end
% end
hold off ; 







% 
% 
% 
% 
% volume_per_voxel = prod(spacing_at_zoom_level_xyz) ;
% minimum_volume_in_voxels = minimum_volume / volume_per_voxel
% maximum_volume_in_voxels = maximum_volume / volume_per_voxel
% 
% is_guess_from_candidate_index = ...
%     (minimum_volume_in_voxels < voxel_count_from_candidate_index) & ...
%     (voxel_count_from_candidate_index < maximum_volume_in_voxels) & ...
%     sqrt_condition_number_from_candidate_index < maximum_sqrt_condition_number ;
% voxel_count_from_guess_index = voxel_count_from_candidate_index(is_guess_from_candidate_index) ;
% xyz_from_guess_index = xyz_from_candidate_index(is_guess_from_candidate_index, :) ;    
% sqrt_condition_number_from_guess_index = sqrt_condition_number_from_candidate_index(is_guess_from_candidate_index) ;
% max_intensity_from_guess_index = max_intensity_from_candidate_index(is_guess_from_candidate_index) ;
% guess_count = length(voxel_count_from_guess_index) ;
% 
% % f = figure('color', 'w') ;
% % a = axes(f, 'YDir', 'reverse') ;
% % image(a, 'CData', substack_mip, ...
% %          'XData', [padded_substack_origin_xyz(1) padded_substack_far_corner_xyz(1)], ...
% %          'YData', [padded_substack_origin_xyz(2) padded_substack_far_corner_xyz(2)], ...
% %          'CDataMapping', 'scaled') ;         
% % colormap(gray(256)) ;
% % axis image
% % 
% % hold on ;
% % candidate_count = size(candidate_centroid_xyz_from_label,1) ;
% % for i = 1 : candidate_count ,
% %     centroid_xyz = candidate_centroid_xyz_from_label(i,:) ;
% %     plot(centroid_xyz(1), centroid_xyz(2), '+', 'Color', [0 0.5 1]) ;    
% % end
% % hold off ; 
% % 
% % hold on ;
% % somata_count = size(putative_somata_xyzs,1) ;
% % for i = 1 : somata_count ,
% %     soma_xyz = putative_somata_xyzs(i,:) ;
% %     plot(soma_xyz(1), soma_xyz(2), 'r+') ;    
% % end
% % hold off ; 
% 
% % miss_1_xy = [74096.2049 18331.3971]  % location of a candidate centroid that is *not* classified as a soma, but should be
% % distance_to_miss_1_xy = sqrt(sum((candidate_centroid_xyz_from_label(:,1:2) - miss_1_xy).^2,2)) ;
% % [distance_from_miss_1_to_nearest_candidate_in_xy, label_of_miss_1] = min(distance_to_miss_1_xy)
% % 
% % centroid_xyz_of_miss_1 = candidate_centroid_xyz_from_label(label_of_miss_1, :)
% % voxel_count_of_miss_1 = voxel_count_from_label(label_of_miss_1) 
% % sqrt_condition_number_of_miss_1 = sqrt_condition_number_from_label(label_of_miss_1)
% % is_putative_soma_for_miss_1 = is_putative_soma_from_label(label_of_miss_1)
% % max_intesity_for_miss_1 = max_intensity_from_label(label_of_miss_1)
% % % sqrt condition number is 2.6, voxel_count is 556, which just clears the
% % % current minimum (515)
% 
% 
% 
% 
% % Calculate distance between each target and each guess
% target_guess_distance_matrix = zeros(target_count, guess_count) ;
% for i = 1 : target_count ,
%     target_i_xyz = xyz_from_target_index(i,:) ;
%     for j = 1 : guess_count ,
%         guess_j_xyz = xyz_from_guess_index(j,:) ;
%         distance = sqrt(sum((guess_j_xyz - target_i_xyz).^2)) ;
%         target_guess_distance_matrix(i,j) = distance ;        
%     end
% end
% 
% [is_target_guess_match, target_guess_match_distance] = match_targets_and_guesses(target_guess_distance_matrix, is_match_distance_threshold) ;
% 
% is_there_a_matched_guess_from_target_index = any(is_target_guess_match, 2) 
% is_there_a_matched_target_from_guess_index = any(is_target_guess_match, 1) 
% [distance_to_matched_guess_from_target_index, matching_guess_index_from_target_index] = min(target_guess_match_distance, [], 2) 
% [distance_to_matched_target_from_guess_index, matching_target_index_from_guess_index] = min(target_guess_match_distance, [], 1) 
% 
% is_hit_from_target_index = is_there_a_matched_guess_from_target_index ;
% is_hit_from_guess_index = is_there_a_matched_target_from_guess_index ;
% is_miss_from_target_index = ~is_there_a_matched_guess_from_target_index ;
% is_chase_from_guess_index = ~is_there_a_matched_target_from_guess_index ;
% 
% guess_hit_count = sum(is_hit_from_target_index)
% assert( sum(is_hit_from_guess_index) == guess_hit_count ) ;
% guess_miss_count = sum(is_miss_from_target_index)
% guess_chase_count = sum(is_chase_from_guess_index)
% 
% guess_precision = guess_hit_count / guess_count 
% guess_recall = guess_hit_count / target_count 
% 
% 
% % Plot the MIP image
% f = figure('color', 'w', 'name', 'targets-and-guesses') ;
% a = axes(f, 'YDir', 'reverse') ;
% % image(a, 'CData', substack_mip, ...
% %          'XData', [padded_substack_origin_xyz(1) padded_substack_far_corner_xyz(1)], ...
% %          'YData', [padded_substack_origin_xyz(2) padded_substack_far_corner_xyz(2)], ...
% %          'CDataMapping', 'scaled') ;         
% colormap(gray(256)) ;
% axis image
% xlabel('x (um)') ;
% ylabel('y (um)') ;
% xlim(heckbert_origin_xyz(1)+[0 stack_shape_xyz(1)]) ;
% ylim(heckbert_origin_xyz(2)+[0 stack_shape_xyz(2)]) ;
% 
% % plot each target, and each guess, and draw a line between matches
% hold on ;
% for target_index = 1 : target_count ,
%     target_xyz = xyz_from_target_index(target_index,:) ;
%     marker_color = fif(is_there_a_matched_guess_from_target_index(target_index), [0 0.5 1], [1 0 0]) ;    
%     plot(target_xyz(1), target_xyz(2), 'Marker', '+', 'Color', marker_color) ;            
%     %text(target_xyz(1)+5, target_xyz(2)+5, sprintf('t%d', target_index), 'Color', 0.5*[1 1 1]) ;            
% end
% for guess_index = 1 : guess_count ,
%     guess_xyz = xyz_from_guess_index(guess_index,:) ;
%     marker_color = fif(is_there_a_matched_target_from_guess_index(guess_index), [0 0.5 1], [1 0 0]) ;    
%     plot(guess_xyz(1), guess_xyz(2), 'Marker', 'o', 'MarkerSize', 6, 'Color', marker_color) ;
%     %text(guess_xyz(1)-5, guess_xyz(2)-5, sprintf('g%d', guess_index), 'Color', 0.5*[1 1 1]) ;            
% end
% for target_index = 1 : target_count ,
%     target_xyz = xyz_from_target_index(target_index,:) ;
%     if is_there_a_matched_guess_from_target_index(target_index) ,
%         guess_index = matching_guess_index_from_target_index(target_index) ;
%         guess_xyz = xyz_from_guess_index(guess_index,:) ;
%         plot([target_xyz(1) guess_xyz(1)], [target_xyz(2) guess_xyz(2)], 'Color', [0 0.5 1]) ;    
%     end
% end
% hold off ; 
% 
% 
% 
% 
% 
% % Considering just the cadidates, plot them in feature space, and
% % characterize each as hit/miss/chase/ball.  This makes sense b/c we cast a wide net: every target is
% % also a candidate.
% 
% is_positive_from_candidate_index =  is_there_a_matched_target_from_candidate_index ;
% is_negative_from_candidate_index = ~is_there_a_matched_target_from_candidate_index ;
% positive_count_within_candidates = sum(is_positive_from_candidate_index) 
% negative_count_within_candidates = sum(is_negative_from_candidate_index) 
% 
% is_test_positive_from_candidate_index =  is_guess_from_candidate_index ;
% is_test_negative_from_candidate_index = ~is_guess_from_candidate_index ;
% test_positive_count_within_candidates = sum(is_test_positive_from_candidate_index) 
% test_negative_count_within_candidates = sum(is_test_negative_from_candidate_index) 
% 
% is_hit_from_candidate_index = is_test_positive_from_candidate_index & is_positive_from_candidate_index ;
% is_miss_from_candidate_index = is_test_negative_from_candidate_index & is_positive_from_candidate_index ;
% is_chase_from_candidate_index = is_test_positive_from_candidate_index & is_negative_from_candidate_index ;
% is_ball_from_candidate_index = is_test_negative_from_candidate_index & is_negative_from_candidate_index ;
% 
% hit_count_within_candidates = sum(is_hit_from_candidate_index)
% miss_count_within_candidates = sum(is_miss_from_candidate_index)
% chase_count_within_candidates = sum(is_chase_from_candidate_index)
% ball_count_within_candidates = sum(is_ball_from_candidate_index) 
% 
% precision_within_candidates = hit_count_within_candidates / test_positive_count_within_candidates 
% recall_within_candidates = hit_count_within_candidates / positive_count_within_candidates 
% 
% hit_rate_within_candidates = hit_count_within_candidates / positive_count_within_candidates
% ball_rate_within_candidates = ball_count_within_candidates / negative_count_within_candidates
% 
% 
% 
% % Plot the targets in the (very wimpy) feature space
% f = figure('Color', 'w', 'Name', 'features') ;
% a = axes(f) ;
% h_hit = plot3(max_intensity_from_candidate_index(is_hit_from_candidate_index), ...
%               voxel_count_from_candidate_index(is_hit_from_candidate_index), ...
%               sqrt_condition_number_from_candidate_index(is_hit_from_candidate_index), ...
%               'Marker', '.', 'Color', [0 0.3 1], 'LineStyle', 'none') ;
% hold on ;
% h_miss = plot3(max_intensity_from_candidate_index(is_miss_from_candidate_index), ...
%                voxel_count_from_candidate_index(is_miss_from_candidate_index), ...
%                sqrt_condition_number_from_candidate_index(is_miss_from_candidate_index), ...
%               'Marker', 'd', 'Color', [1 0 0], 'LineStyle', 'none') ;
% h_chase = plot3(max_intensity_from_candidate_index(is_chase_from_candidate_index), ...
%                 voxel_count_from_candidate_index(is_chase_from_candidate_index), ...
%                 sqrt_condition_number_from_candidate_index(is_chase_from_candidate_index), ...
%               'Marker', 'd', 'Color', [138 43 226]/255, 'LineStyle', 'none') ;
% h_ball = plot3(max_intensity_from_candidate_index(is_ball_from_candidate_index), ...
%                voxel_count_from_candidate_index(is_ball_from_candidate_index), ...
%                sqrt_condition_number_from_candidate_index(is_ball_from_candidate_index), ...
%               'Marker', '.', 'Color', [0 0.8 0], 'LineStyle', 'none') ;
% hold off ;
% xlabel('Max GFP signal (counts)') ;
% ylabel('Voxel count') ;
% zlabel('SD ratio') ;
% %xlim([0 2^16]) ;
% %ylim([0 2^16]) ;
% handles = zeros(1,0) ;
% legend_labels = cell(1,0) ;
% if ~isempty(h_hit) ,
%     handles(1, end+1) = h_hit ;
%     legend_labels{1, end+1} = 'hit' ;
% end
% if ~isempty(h_miss) ,
%     handles(1, end+1) = h_miss ;
%     legend_labels{1, end+1} = 'miss' ;
% end
% if ~isempty(h_chase) ,
%     handles(1, end+1) = h_chase ;
%     legend_labels{1, end+1} = 'chase' ;
% end
% if ~isempty(h_ball) ,
%     handles(1, end+1) = h_ball ;
%     legend_labels{1, end+1} = 'ball' ;
% end
% 
% legend(handles, legend_labels, 'Location', 'northwest') ;
% axis vis3d
% grid on
% camproj perspective
% 
% 
% 
