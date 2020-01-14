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
intensity_threshold = 0.8 * 2^16 ;
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
                             @pad_and_find_somata_then_save, ...
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
                    pad_and_find_somata_then_save(somata_mat_file_path, ...
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
xyz_from_guess_index = zeros(0,3) ;
for i = 1 : somata_file_count ,
    somata_mat_file_name = soma_mat_file_names{i} ;
    somata_mat_file_path = fullfile(chunks_folder_path, somata_mat_file_name) ;
    s = load(somata_mat_file_path) ;    
    xyz_from_guess_index = vertcat(xyz_from_guess_index, s.xyz_from_guess_index) ;  %#ok<AGROW>
end
guess_count = size(xyz_from_guess_index, 1) ;

%%
% Save as a .swc file, which can hold a forest, it turns out
forest_name = sprintf('%s-auto-somata', sample_date) ;
forest_color = [1 0 0] ;
swc_file_name = horzcat(forest_name, '.swc') ;
swc_file_path = fullfile(somata_folder_path, swc_file_name) ;
save_somata_as_single_swc(swc_file_path, xyz_from_guess_index, forest_name, forest_color) ;


% Load the ground-truth (these are all soma locations)
xyz_from_target_index = load_soma_targets() ;
target_count = size(xyz_from_target_index, 1) ;

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

is_match_distance_threshold = 10 ;  % um
[is_target_guess_match, target_guess_match_distance] = match_targets_and_guesses(target_guess_distance_matrix, is_match_distance_threshold) ;

is_there_a_matched_guess_from_target_index = any(is_target_guess_match, 2) ;
is_there_a_matched_target_from_guess_index = any(is_target_guess_match, 1) ;
[distance_to_matched_guess_from_target_index, matching_guess_index_from_target_index] = min(target_guess_match_distance, [], 2) ;
[distance_to_matched_target_from_guess_index, matching_target_index_from_guess_index] = min(target_guess_match_distance, [], 1) ;

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


% Make an axes that would accomodate the a MIP of the full stack
f = figure('color', 'w', 'name', 'targets-and-guesses') ;
a = axes(f, 'YDir', 'reverse') ;
% image(a, 'CData', substack_mip, ...
%          'XData', [padded_substack_origin_xyz(1) padded_substack_far_corner_xyz(1)], ...
%          'YData', [padded_substack_origin_xyz(2) padded_substack_far_corner_xyz(2)], ...
%          'CDataMapping', 'scaled') ;         
colormap(gray(256)) ;
axis image
xlim(heckbert_origin_xyz(1)+[0 stack_shape_xyz(1)]) ;
ylim(heckbert_origin_xyz(2)+[0 stack_shape_xyz(2)]) ;
xlabel('x (um)') ;
ylabel('y (um)') ;

% plot each target, and each guess, and draw a line between matches
hold on ;
for target_index = 1 : target_count ,
    target_xyz = xyz_from_target_index(target_index,:) ;
    marker_color = fif(is_there_a_matched_guess_from_target_index(target_index), [0 0.5 1], [1 0 0]) ;    
    plot(target_xyz(1), target_xyz(2), 'Marker', '+', 'Color', marker_color) ;            
    text(target_xyz(1)+5, target_xyz(2)+5, sprintf('t%d', target_index), 'Color', 0.5*[1 1 1]) ;            
end
for guess_index = 1 : guess_count ,
    guess_xyz = xyz_from_guess_index(guess_index,:) ;
    marker_color = fif(is_there_a_matched_target_from_guess_index(guess_index), [0 0.5 1], [1 0 0]) ;    
    plot(guess_xyz(1), guess_xyz(2), 'Marker', 'o', 'MarkerSize', 6, 'Color', marker_color) ;
    text(guess_xyz(1)-5, guess_xyz(2)-5, sprintf('g%d', guess_index), 'Color', 0.5*[1 1 1]) ;            
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
