parameters = struct('intensity_threshold', {intensity_threshold}, ...
                    'minimum_volume', {minimum_volume}, ...
                    'maximum_volume', {maximum_volume}, ...
                    'maximum_sqrt_condition_number', maximum_sqrt_condition_number) ;
                
% this_file_path = mfilename('fullpath') ;
% this_folder_path = fileparts(this_file_path) ;                
% output_folder_path = fullfile(this_folder_path, sprintf('auto-somata-%s-%s', sample_date, tag)) ;

                
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
if do_force_computation ,
    system(sprintf('rm -rf %s', output_folder_path)) ;
end
if ~exist(output_folder_path, 'file') ,
    mkdir(output_folder_path) ;
end
chunks_folder_path = fullfile(output_folder_path, 'chunks') ;
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
                                 background_channel_index, ...
                                 zoom_level, ...
                                 origin_at_zoom_level_xyz, ...
                                 spacing_at_zoom_level_xyz, ...
                                 chunk_offset_within_stack_ijk1, ...
                                 analysis_chunk_shape_ijk, ...
                                 pad_depth_in_um, ...
                                 parameters) ;
                else
                    tic_id = tic() ;
                    pad_and_find_candidate_somata_then_save(somata_mat_file_path, ...
                                                            rendered_folder_path, ...
                                                            gfp_channel_index, ...
                                                            background_channel_index, ...
                                                            zoom_level, ...
                                                            origin_at_zoom_level_xyz, ...
                                                            spacing_at_zoom_level_xyz, ...
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
feature_struct_from_candidate_index = compute_derived_component_features(zeros(0,1)) ;
for i = 1 : somata_file_count ,
    somata_mat_file_name = soma_mat_file_names{i} ;
    somata_mat_file_path = fullfile(chunks_folder_path, somata_mat_file_name) ;
    s = load(somata_mat_file_path) ;        
    %xyz_from_guess_index = vertcat(xyz_from_guess_index, s.xyz_from_guess_index) ;  %#ok<AGROW>
    %feature_struct_from_guess_index = vertcat(feature_struct_from_guess_index, s.feature_struct_from_guess_index) ;  %#ok<AGROW>
    feature_struct_from_candidate_index = vertcat(feature_struct_from_candidate_index, s.feature_struct_from_candidate_index) ;  %#ok<AGROW>
end
%guess_count = size(xyz_from_guess_index, 1) ;
candidate_count = length(feature_struct_from_candidate_index) ;

% break out the individual fields
voxel_count_from_candidate_index = [feature_struct_from_candidate_index.voxel_count]' ;
sqrt_condition_number_from_candidate_index = [feature_struct_from_candidate_index.sqrt_condition_number]' ;
max_intensity_from_candidate_index = [feature_struct_from_candidate_index.max_intensity]' ;
max_background_intensity_from_candidate_index = [feature_struct_from_candidate_index.max_background_intensity]' ;

% Get a mip of the whole brain, for visualization
overview_zoom_level = 2 ;
stack_shape_at_overview_zoom_level_ijk = chunk_shape_ijk * 2^overview_zoom_level ;
spacing_at_overview_zoom_level_xyz = spacing_at_zoom_level_0_xyz ./ (2^overview_zoom_level) ;
origin_at_overview_zoom_level_xyz = heckbert_origin_xyz + spacing_at_overview_zoom_level_xyz/2 ;
whole_brain_at_overview_zoom_level_yxz = ...
   get_mouselight_rendered_substack(rendered_folder_path, gfp_channel_index, [1 1 1], stack_shape_at_overview_zoom_level_ijk, overview_zoom_level) ;  
overview_mip_yx = max(whole_brain_at_overview_zoom_level_yxz, [], 3) ;
mip_origin_xy = origin_at_overview_zoom_level_xyz(1:2) ;
mip_spacing_xy = spacing_at_overview_zoom_level_xyz(1:2) ;

% % Load the ground-truth (these are all soma locations)
% %[xyz_from_target_index, name_from_target_index] = load_traceable_soma_targets_from_tracers() ;
% xyz_from_target_index = zeros(0,3) ;
% name_from_target_index = cell(0,1) ;
% target_count = size(xyz_from_target_index, 1) ;
% 
% % Report the perf of the candidates
% do_plot_candidates = false ;
% [matching_candidate_index_from_target_index, matching_target_index_from_candidate_index] = ...
%     print_performace_statistics_and_plot('candidates', ...
%                                          xyz_from_target_index, ...
%                                          feature_struct_from_candidate_index, ...
%                                          heckbert_origin_xyz, ...
%                                          stack_shape_xyz, ...
%                                          spacing_at_zoom_level_xyz, ...
%                                          overview_mip_yx, ...
%                                          mip_origin_xy, ...
%                                          mip_spacing_xy, ...
%                                          [], ...
%                                          do_plot_candidates) ;




volume_per_voxel = prod(spacing_at_zoom_level_xyz) ;
minimum_volume_in_voxels = minimum_volume / volume_per_voxel
maximum_volume_in_voxels = maximum_volume / volume_per_voxel

is_guess_from_candidate_index = ...
    (minimum_volume_in_voxels < voxel_count_from_candidate_index) & ...
    (voxel_count_from_candidate_index < maximum_volume_in_voxels) & ...
    sqrt_condition_number_from_candidate_index < maximum_sqrt_condition_number & ...
    max_intensity_from_candidate_index > max_background_intensity_from_candidate_index ;
feature_struct_from_guess_index = feature_struct_from_candidate_index(is_guess_from_candidate_index) ;
guess_count = length(feature_struct_from_guess_index) ;

% % Save auto-somata as a .swc file, which can hold a forest, it turns out
% forest_name = sprintf('%s-%s-auto-somata', sample_date, tag) ;
% forest_color = [1 0 0] ;
% swc_file_name = horzcat(forest_name, '.swc') ;
% swc_file_path = fullfile(output_folder_path, swc_file_name) ;
% save_somata_as_single_swc(swc_file_path, xyz_from_guess_index, forest_name, forest_color) ;
% save_somata_as_multiple_swcs(


% output the predictions
centroidoid_xyz_from_guess_index = reshape([feature_struct_from_guess_index(:).centroidoid_xyz], [3 guess_count])' ;
name_template = 'soma-prediction-%d' ;
swc_folder_path = fullfile(output_folder_path, 'swcs') ;
if ~exist(swc_folder_path, 'file') ,
    mkdir(swc_folder_path) ,
end
output_swc_file_name_template = fullfile(swc_folder_path, 'soma-prediction-%d.swc') ;
color = [1 0 1] ;  % magenta
%save_somata_as_single_swc(output_swc_file_name, centroidoid_xyz_from_guess_index, name, color)
save_somata_as_multiple_swcs(output_swc_file_name_template, centroidoid_xyz_from_guess_index, name_template, color) ;


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


% % Report the perf of the guesses
% do_plot_candidates = true ;
% mip_clim = [11000 2^16-1] ;
% matching_guess_index_from_target_index = ...
% print_performace_statistics_and_plot('guesses', ...
%                                      xyz_from_target_index, ...
%                                      feature_struct_from_guess_index, ...
%                                      heckbert_origin_xyz, ...
%                                      stack_shape_xyz, ...
%                                      spacing_at_zoom_level_xyz, ...
%                                      overview_mip_yx, ...
%                                      mip_origin_xy, ...
%                                      mip_spacing_xy, ...
%                                      mip_clim, ...
%                                      do_plot_candidates) ;
%                                  
%                                  
% %%                                 
% % Considering just the cadidates, plot them in feature space, and
% % characterize each as hit/miss/chase/ball.  This makes sense b/c we cast a wide net: every target is
% % also a candidate.
% 
% is_there_a_matched_candidate_from_target_index = isfinite(matching_candidate_index_from_target_index) ;
% is_there_a_matched_target_from_candidate_index = isfinite(matching_target_index_from_candidate_index) ;
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
% % break out the individual fields
% voxel_count_from_guess_index = [feature_struct_from_guess_index.voxel_count]' ;
% sqrt_condition_number_from_guess_index = [feature_struct_from_guess_index.sqrt_condition_number]' ;
% max_intensity_from_guess_index = [feature_struct_from_guess_index.max_intensity]' ;
% 
% 
% % Plot the targets in a projection of the feature space
% f = figure('Color', 'w', 'Name', 'features') ;
% a = axes(f) ;
% h_ball = [] ;
% % h_ball = plot3(max_intensity_from_candidate_index(is_ball_from_candidate_index), ...
% %                voxel_count_from_candidate_index(is_ball_from_candidate_index), ...
% %                sqrt_condition_number_from_candidate_index(is_ball_from_candidate_index), ...
% %               'Marker', '.', 'Color', [0 0.8 0], 'LineStyle', 'none') ;
% h_hit = plot3(max_intensity_from_candidate_index(is_hit_from_candidate_index), ...
%               voxel_count_from_candidate_index(is_hit_from_candidate_index), ...
%               sqrt_condition_number_from_candidate_index(is_hit_from_candidate_index), ...
%               'Marker', '.', 'Color', [0 0.3 1], 'LineStyle', 'none') ;
% hold on ;
% h_chase = plot3(max_intensity_from_candidate_index(is_chase_from_candidate_index), ...
%                 voxel_count_from_candidate_index(is_chase_from_candidate_index), ...
%                 sqrt_condition_number_from_candidate_index(is_chase_from_candidate_index), ...
%               'Marker', 'd', 'Color', [138 43 226]/255, 'LineStyle', 'none') ;
% h_miss = plot3(max_intensity_from_candidate_index(is_miss_from_candidate_index), ...
%                voxel_count_from_candidate_index(is_miss_from_candidate_index), ...
%                sqrt_condition_number_from_candidate_index(is_miss_from_candidate_index), ...
%               'Marker', 'd', 'Color', [1 0 0], 'LineStyle', 'none') ;
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
% %axis vis3d
% grid on
% %camproj perspective
% 
% 
% 
% % Plot the targets another projection of the feature space
% f = figure('Color', 'w', 'Name', 'features') ;
% a = axes(f) ;
% h_hit = plot(max_intensity_from_candidate_index(is_hit_from_candidate_index), ...
%              max_background_intensity_from_candidate_index(is_hit_from_candidate_index), ...
%               'Marker', '.', 'Color', [0 0.3 1], 'LineStyle', 'none') ;
% hold on ;
% h_miss = plot(max_intensity_from_candidate_index(is_miss_from_candidate_index), ...
%               max_background_intensity_from_candidate_index(is_miss_from_candidate_index), ...
%               'Marker', 'd', 'Color', [1 0 0], 'LineStyle', 'none') ;
% h_chase = plot(max_intensity_from_candidate_index(is_chase_from_candidate_index), ...
%                max_background_intensity_from_candidate_index(is_chase_from_candidate_index), ...
%               'Marker', 'd', 'Color', [138 43 226]/255, 'LineStyle', 'none') ;
% h_ball = []  ;        
% % h_ball = plot(max_intensity_from_candidate_index(is_ball_from_candidate_index), ...
% %               max_background_intensity_from_candidate_index(is_ball_from_candidate_index), ...
% %               'Marker', '.', 'Color', [0 0.8 0], 'LineStyle', 'none') ;
% hold off ;
% xlabel('Max GFP signal (counts)') ;
% ylabel('Max background signal (counts)') ;
% axis equal
% xlim([0 2^16]) ;
% ylim([0 2^16]) ;
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
% 
% grid on
% 
% 
