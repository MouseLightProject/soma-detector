close all
clear

this_file_path = mfilename('fullpath') ;
this_folder_path = fileparts(this_file_path) ;
swc_folder_path = fullfile(this_folder_path, '2019-10-04_adam_AT') ;

soma_swc_file_name_template = 'Neuron*.swc' ;
soma_xyzs = load_somata_xyzs_from_swc_files(swc_folder_path, soma_swc_file_name_template) ;

distractor_swc_file_name_template = '2019-10-04-auto-somata*.swc' ;
distractor_xyzs = load_somata_xyzs_from_swc_files(swc_folder_path, distractor_swc_file_name_template) ;

% Put them all into one array
target_xyzs = vertcat(soma_xyzs, ...
                      distractor_xyzs) ;
soma_count = size(soma_xyzs, 1) ;
distractor_count = size(distractor_xyzs, 1) ;
target_labels = vertcat(ones(soma_count, 1) , ...
                        zeros(distractor_count, 1)) ;
target_count = length(target_labels) ;
                    
% load the classifier's detected somata
auto_somata_file_path = fullfile(this_folder_path, 'auto-somata', '2019-10-04-auto-somata.swc') ;
auto_somata_swc_array = load_swc(auto_somata_file_path) ;
auto_soma_xyzs = auto_somata_swc_array(:, 3:5) ;
auto_soma_count = size(auto_soma_xyzs, 1) ;



%%
% Calculate distances between all pairs
distance_matrix = zeros(target_count, auto_soma_count) ;
for i = 1 : target_count ,
    target_i_xyz = target_xyzs(i,:) ;
    for j = 1 : auto_soma_count ,
        auto_soma_j_xyz = auto_soma_xyzs(j,:) ;
        distance = sqrt(sum((auto_soma_j_xyz - target_i_xyz).^2)) ;
        distance_matrix(i,j) = distance ;        
    end
end

is_close_matrix = (distance_matrix<10) ;
is_test_positive_from_target_index = any(is_close_matrix,2) ;

is_hit = target_labels & is_test_positive_from_target_index ;
is_miss = target_labels & ~is_test_positive_from_target_index ;
is_ball = ~target_labels & ~is_test_positive_from_target_index ;  % there won't be any of these, b/c of the way the targets were consructed
is_chase = ~target_labels & is_test_positive_from_target_index ;
is_positive = logical(target_labels) ;
is_negative = ~target_labels ;


hit_rate = sum(is_hit) / sum(is_positive) 
miss_rate = sum(is_miss) / sum(is_positive) 
ball_rate  = sum(is_ball) / sum(is_negative) 
chase_rate = sum(is_chase) / sum(is_negative)


synthetic_ball_count = 1000 ;
hit_rate = sum(is_hit) / sum(is_positive) 
miss_rate = sum(is_miss) / sum(is_positive) 
ball_rate  = (sum(is_ball) + synthetic_ball_count) / (sum(is_negative) + synthetic_ball_count)
chase_rate = sum(is_chase) / (sum(is_negative) + synthetic_ball_count)


%
% Look at the image data at the targets
sample_date = '2019-10-04' ;
rendered_folder_path = sprintf('/nrs/mouselight/SAMPLES/%s', sample_date) ;
gfp_channel_index = 0 ;
zoom_level = 4 ;  % The zoom level of the tiles we will analyze
pad_depth_in_um = 50 ; % um
intensity_threshold = 0.8 * 2^16 ;
minimum_volume = 3000 ;  % um^3
maximum_volume = 5*minimum_volume ;  % um^3
maximum_condition_number = 5 ;
parameters = struct('intensity_threshold', {intensity_threshold}, ...
                    'minimum_volume', {minimum_volume}, ...
                    'maximum_volume', {maximum_volume}, ...
                    'maximum_condition_number', maximum_condition_number) ;
                
this_file_path = mfilename('fullpath') ;
this_folder_path = fileparts(this_file_path) ;                
somata_folder_path = sprintf(this_folder_path, 'auto-somata') ;

                
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
analysis_chunk_shape_ijk = 4*chunk_shape_ijk ;
chunks_per_dimension_ijk = stack_shape_ijk ./ analysis_chunk_shape_ijk ;
stack_shape_xyz = stack_shape_ijk .* spacing_at_zoom_level_xyz ;  % the shape of the stack, in real units (doesn't change with zoom)




target_fluorescence_values = zeros(target_count, 2) ;
for i = 1 : target_count ,
    target_xyz = target_xyzs(i,:) ;
    target_ijk1 = round((target_xyz - origin_at_zoom_level_xyz) ./ spacing_at_zoom_level_xyz) + 1 ;
    target_fluorescence_values(i,1) = ...
        get_mouselight_rendered_substack(rendered_folder_path, 0, target_ijk1, [1 1 1], zoom_level) ;  
    target_fluorescence_values(i,2) = ...
        get_mouselight_rendered_substack(rendered_folder_path, 1, target_ijk1, [1 1 1], zoom_level) ;  
end

% plot those
f = figure('Color', 'w', 'Name', 'counts') ;
a = axes(f) ;
h_hit = plot(target_fluorescence_values(is_hit,1), target_fluorescence_values(is_hit,2), 'b+') ;
hold on ;
h_miss = plot(target_fluorescence_values(is_miss,1), target_fluorescence_values(is_miss,2), 'bo') ;
h_ball = plot(target_fluorescence_values(is_ball,1), target_fluorescence_values(is_ball,2), 'r+') ;
h_chase = plot(target_fluorescence_values(is_chase,1), target_fluorescence_values(is_chase,2), 'ro') ;
plot(intensity_threshold*[1 1], [0 2^16], ':k') ;
hold off ;
xlabel('GFP signal (counts)') ;
ylabel('Background signal (counts)') ;
axis equal ;
xlim([0 2^16]) ;
ylim([0 2^16]) ;
legend([h_hit h_miss h_chase], {'hit' 'miss', 'chase'}, 'Location', 'southeast') ;

% What fraction of the positive targets are above threshold?
fraction_positive_targets_above_threshold = mean(target_fluorescence_values(is_positive,1)>intensity_threshold)



    