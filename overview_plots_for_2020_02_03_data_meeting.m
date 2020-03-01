do_force_computation = false ;
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

% intensity_threshold = 40000 ;
% minimum_volume = 500 ;  % um^3
% maximum_volume = 15000 ;  % um^3
% maximum_sqrt_condition_number = 10 ;

intensity_threshold = 0.75 * 2^16 ;
minimum_volume = 500 ;  % um^3
%maximum_volume = 15000 ;  % um^3
maximum_volume = 25000 ;  % um^3
%maximum_sqrt_condition_number = 10 ;
maximum_sqrt_condition_number = 20 ;

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
stack_shape_xyz = stack_shape_ijk .* spacing_at_zoom_level_xyz  %#ok<NOPTS> 
  % the shape of the full-brain stack in um, does not change with zoom level
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
tic_id = tic() ;
job_count = prod(chunks_per_dimension_ijk) ;
job_id_array = repmat(-1, [1 job_count]) ;
job_index = 1 ;
for chunk_i = 1 : chunks_per_dimension_ijk(1) ,
    for chunk_j = 1 : chunks_per_dimension_ijk(2) ,
        for chunk_k = 1 : chunks_per_dimension_ijk(3) ,
            somata_mat_file_name = sprintf('somata-for-chunk-%d-%d-%d.mat', chunk_i, chunk_j, chunk_k) ;
            somata_mat_file_path = fullfile(chunks_folder_path, somata_mat_file_name) ;
            chunk_offset_within_chunks_ijk1 = [chunk_i chunk_j chunk_k]   %#ok<NOPTS>
            chunk_offset_within_stack_ijk1 = (chunk_offset_within_chunks_ijk1-1) .* analysis_chunk_shape_ijk + 1 ;
            if ~exist(somata_mat_file_path, 'file') ,
                if do_use_bsub , 
                    bsub_options = sprintf(bsub_options_template, chunk_i, chunk_j, chunk_k, chunk_i, chunk_j, chunk_k)   %#ok<NOPTS>
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

toc(tic_id) ;

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

% %%
% % Save as a .swc file, which can hold a forest, it turns out
% forest_name = sprintf('%s-auto-somata', sample_date) ;
% forest_color = [1 0 0] ;
% swc_file_name = horzcat(forest_name, '.swc') ;
% swc_file_path = fullfile(somata_folder_path, swc_file_name) ;
% save_somata_as_single_swc(swc_file_path, xyz_from_guess_index, forest_name, forest_color) ;

% Get a mip of the whole brain, for visualization
overview_zoom_level = 2 ;
stack_shape_at_overview_zoom_level_ijk = chunk_shape_ijk * 2^overview_zoom_level ;
spacing_at_overview_zoom_level_xyz = spacing_at_zoom_level_0_xyz ./ (2^overview_zoom_level) ;
origin_at_overview_zoom_level_xyz = heckbert_origin_xyz + spacing_at_overview_zoom_level_xyz/2 ;
whole_brain_at_overview_zoom_level_yxz = ...
   get_mouselight_rendered_substack(rendered_folder_path, gfp_channel_index, [1 1 1], stack_shape_at_overview_zoom_level_ijk, overview_zoom_level) ;  
overview_mip_yx = max(whole_brain_at_overview_zoom_level_yxz, [], 3) ;
%overview_mip_yx = whole_brain_at_overview_zoom_level_yxz(:,:,round((end+1)/2)) ;
mip_origin_xy = origin_at_overview_zoom_level_xyz(1:2) ;
mip_spacing_xy = spacing_at_overview_zoom_level_xyz(1:2) ;
mip_stack_heckbert_origin_xyz = heckbert_origin_xyz ;


% Load the ground-truth (these are all soma locations)
[xyz_from_target_index, name_from_target_index] = load_traceable_soma_targets_from_tracers() ;
% % Patch for bad G-040 soma
% if norm(xyz_from_target_index(40,:)-[ 68289.340023  18112.354727  36828.363892 ]) < 1e-6 ,
%     xyz_from_target_index(40,:) = [69232.8, 18467.3, 36351.0] ;
% end
target_count = size(xyz_from_target_index, 1) ;

% Compute the guesses
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
do_plot_candidates = true ;
mip_clim = [11000 2^16-1] ;
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
                                 

% Figure out what targets are in what guesss
%target_count = size(xyz_from_target_index, 1) ;
%guess_count = length(feature_struct_from_guess_index) ;
origin_at_zoom_level_xyz = heckbert_origin_xyz + spacing_at_zoom_level_xyz/2 ;
matching_guess_index_from_target_index = ...
    match_targets_and_components(xyz_from_target_index, feature_struct_from_guess_index, origin_at_zoom_level_xyz, spacing_at_zoom_level_xyz) ;
matching_target_index_from_guess_index = invert_partial_map_array(matching_guess_index_from_target_index, guess_count) ;
is_there_a_matched_guess_from_target_index = isfinite(matching_guess_index_from_target_index) ;
is_there_a_matched_target_from_guess_index = isfinite(matching_target_index_from_guess_index) ;

% Compute hit counts, precision, recall for the guesss
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
f = figure('color', 'k', 'name', 'overview-mip') ;
a = axes(f, 'YDir', 'reverse', 'XColor', 'w', 'YColor', 'w', 'Layer', 'top', 'Box', 'on', 'DataAspectRatio', [1 1 1]) ;
mip_shape_ji = size(overview_mip_yx) ;
mip_shape_ij = mip_shape_ji([2 1]) ;
mip_shape_xy = mip_spacing_xy .* mip_shape_ij ;
mip_near_voxel_center_xy = mip_spacing_xy/2 ;
mip_far_voxel_center_xy = mip_shape_xy - mip_spacing_xy/2 ;
%mip_far_corner_xy = mip_origin_xy + mip_spacing_xy .* (mip_shape_ij-1) ;
image(a, 'CData', overview_mip_yx, ...
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
f = figure('color', 'k', 'name', 'overview-targets-and-guesses-mip') ;
a = axes(f, 'YDir', 'reverse', 'XColor', 'w', 'YColor', 'w', 'Layer', 'top', 'Box', 'on', 'DataAspectRatio', [1 1 1]) ;
mip_shape_ji = size(overview_mip_yx) ;
mip_shape_ij = mip_shape_ji([2 1]) ;
mip_shape_xy = mip_spacing_xy .* mip_shape_ij ;
mip_near_voxel_center_xy = mip_spacing_xy/2 ;
mip_far_voxel_center_xy = mip_shape_xy - mip_spacing_xy/2 ;
%mip_far_corner_xy = mip_origin_xy + mip_spacing_xy .* (mip_shape_ij-1) ;
image(a, 'CData', overview_mip_yx, ...
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
    plot(target_xyz(1)-mip_stack_heckbert_origin_xyz(1), ...
         target_xyz(2)-mip_stack_heckbert_origin_xyz(2), ...
         'Marker', '+', 'Color', marker_color) ;
    %text(target_xyz(1)+5, target_xyz(2)+5, sprintf('t%d', target_index), 'Color', 0.5*[1 1 1]) ;            
end
if true ,
    for guess_index = 1 : guess_count ,
        guess_xyz = feature_struct_from_guess_index(guess_index).centroidoid_xyz ;
        marker_color = fif(is_there_a_matched_target_from_guess_index(guess_index), [0 0.5 1], [1 0 0]) ;    
        plot(guess_xyz(1)-mip_stack_heckbert_origin_xyz(1), ...
             guess_xyz(2)-mip_stack_heckbert_origin_xyz(2), ...
             'Marker', 'o', 'MarkerSize', 6, 'Color', marker_color) ;
        %text(guess_xyz(1)-5, guess_xyz(2)-5, sprintf('g%d', guess_index), 'Color', 0.5*[1 1 1]) ;            
    end
    for target_index = 1 : target_count ,
        target_xyz = xyz_from_target_index(target_index,:) ;
        if is_there_a_matched_guess_from_target_index(target_index) ,
            guess_index = matching_guess_index_from_target_index(target_index) ;
            guess_xyz = feature_struct_from_guess_index(guess_index).centroidoid_xyz ;
            plot([target_xyz(1) guess_xyz(1)]-mip_stack_heckbert_origin_xyz(1), ...
                 [target_xyz(2) guess_xyz(2)]-mip_stack_heckbert_origin_xyz(2), ...
                 'Color', [0 0.5 1]) ;    
        end
    end
end
hold off ; 



%%                                 
%
% Show a handful of examples of hits, misses, chases
%
maximum_example_count_per_class = 4 ;
desired_example_stack_shape_xyz = [150 150 150] ;
example_stack_half_shape_ijk = round((desired_example_stack_shape_xyz/2) ./ spacing_at_max_zoom_xyz) ;
example_stack_shape_ijk = 2 * example_stack_half_shape_ijk + 1 ;
example_stack_shape_xyz = example_stack_shape_ijk .* spacing_at_max_zoom_xyz ;

% Hits
target_index_from_hit_index = find(is_there_a_matched_guess_from_target_index) ;
target_index_from_example_hit_index = target_index_from_hit_index(randperm(hit_count, maximum_example_count_per_class)) ;
xyz_from_example_hit_index = xyz_from_target_index(target_index_from_example_hit_index, :) ;
mip_clim = [11000 2^16-1] ;

for i = 1 : size(xyz_from_example_hit_index, 1) ,
    center_xyz = xyz_from_example_hit_index(i, :) ;
    center_ijk1 = round((center_xyz-origin_at_max_zoom_xyz) ./ spacing_at_max_zoom_xyz)+1 ;
    corner_ijk1 = center_ijk1 - example_stack_half_shape_ijk ;
    example_substack_yxz = ...
        get_mouselight_rendered_substack(rendered_folder_path, gfp_channel_index, corner_ijk1, example_stack_shape_ijk, max_zoom_level) ;  
    example_mip_yx = max(example_substack_yxz, [], 3) ;
    mip_spacing_xy = spacing_at_max_zoom_xyz(1:2) ;
    
    
    % Plot the example MIP image
    f = figure('color', 'k', 'name', sprintf('hit-example-%d', i)) ;
    a = axes(f, 'YDir', 'reverse', 'XColor', 'w', 'YColor', 'w', 'Layer', 'top', 'Box', 'on', 'DataAspectRatio', [1 1 1]) ;
    mip_shape_ji = size(example_mip_yx) ;
    mip_shape_ij = mip_shape_ji([2 1]) ;
    mip_shape_xy = mip_spacing_xy .* mip_shape_ij ;
    mip_near_voxel_center_xy = mip_spacing_xy/2 ;
    mip_far_voxel_center_xy = mip_shape_xy - mip_spacing_xy/2 ;
    %mip_far_corner_xy = mip_origin_xy + mip_spacing_xy .* (mip_shape_ij-1) ;
    image(a, 'CData', example_mip_yx, ...
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
end


% Misses
target_index_from_miss_index = find(~is_there_a_matched_guess_from_target_index) ;
target_index_from_example_miss_index = target_index_from_miss_index(randperm(miss_count, min(miss_count, maximum_example_count_per_class))) ;
xyz_from_example_miss_index = xyz_from_target_index(target_index_from_example_miss_index, :) ;
mip_clim = [11000 2^16-1] ;

for i = 1 : size(xyz_from_example_miss_index, 1) ,
    center_xyz = xyz_from_example_miss_index(i, :) ;
    center_ijk1 = round((center_xyz-origin_at_max_zoom_xyz) ./ spacing_at_max_zoom_xyz)+1 ;
    corner_ijk1 = center_ijk1 - example_stack_half_shape_ijk ;
    example_substack_yxz = ...
        get_mouselight_rendered_substack(rendered_folder_path, gfp_channel_index, corner_ijk1, example_stack_shape_ijk, max_zoom_level) ;  
    example_mip_yx = max(example_substack_yxz, [], 3) ;
    mip_spacing_xy = spacing_at_max_zoom_xyz(1:2) ;
    
    
    % Plot the example MIP image
    f = figure('color', 'k', 'name', sprintf('miss-example-%d', i)) ;
    a = axes(f, 'YDir', 'reverse', 'XColor', 'w', 'YColor', 'w', 'Layer', 'top', 'Box', 'on', 'DataAspectRatio', [1 1 1]) ;
    mip_shape_ji = size(example_mip_yx) ;
    mip_shape_ij = mip_shape_ji([2 1]) ;
    mip_shape_xy = mip_spacing_xy .* mip_shape_ij ;
    mip_near_voxel_center_xy = mip_spacing_xy/2 ;
    mip_far_voxel_center_xy = mip_shape_xy - mip_spacing_xy/2 ;
    %mip_far_corner_xy = mip_origin_xy + mip_spacing_xy .* (mip_shape_ij-1) ;
    image(a, 'CData', example_mip_yx, ...
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
end



                                 
% Chases
xyz_from_guess_index = reshape([feature_struct_from_guess_index.centroidoid_xyz], [3 guess_count])' ;
guess_index_from_miss_index = find(~is_there_a_matched_target_from_guess_index) ;
guess_index_from_example_miss_index = guess_index_from_miss_index(randperm(chase_count, min(chase_count, maximum_example_count_per_class))) ;
xyz_from_example_chase_index = xyz_from_guess_index(guess_index_from_example_miss_index, :) ;
mip_clim = [11000 2^16-1] ;

for i = 1 : size(xyz_from_example_chase_index, 1) ,
    center_xyz = xyz_from_example_chase_index(i, :) ;
    center_ijk1 = round((center_xyz-origin_at_max_zoom_xyz) ./ spacing_at_max_zoom_xyz)+1 ;
    corner_ijk1 = center_ijk1 - example_stack_half_shape_ijk ;
    example_substack_yxz = ...
        get_mouselight_rendered_substack(rendered_folder_path, gfp_channel_index, corner_ijk1, example_stack_shape_ijk, max_zoom_level) ;  
    example_mip_yx = max(example_substack_yxz, [], 3) ;
    mip_spacing_xy = spacing_at_max_zoom_xyz(1:2) ;
    
    
    % Plot the example MIP image
    f = figure('color', 'k', 'name', sprintf('chase-example-%d', i)) ;
    a = axes(f, 'YDir', 'reverse', 'XColor', 'w', 'YColor', 'w', 'Layer', 'top', 'Box', 'on', 'DataAspectRatio', [1 1 1]) ;
    mip_shape_ji = size(example_mip_yx) ;
    mip_shape_ij = mip_shape_ji([2 1]) ;
    mip_shape_xy = mip_spacing_xy .* mip_shape_ij ;
    mip_near_voxel_center_xy = mip_spacing_xy/2 ;
    mip_far_voxel_center_xy = mip_shape_xy - mip_spacing_xy/2 ;
    %mip_far_corner_xy = mip_origin_xy + mip_spacing_xy .* (mip_shape_ij-1) ;
    image(a, 'CData', example_mip_yx, ...
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
end



                                 
