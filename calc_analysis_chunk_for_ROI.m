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

xyz_of_interest = [73526.6, 18412.8, 36071.0] 
ijk1_of_interest = round((xyz_of_interest-origin_at_zoom_level_xyz) ./ spacing_at_zoom_level_xyz) +1 ;
analysis_chunk_of_interest_ijk0_fractional = (ijk1_of_interest-1) ./ analysis_chunk_shape_ijk
analysis_chunk_of_interest_ijk0 = floor(analysis_chunk_of_interest_ijk0_fractional) 
analysis_chunk_of_interest_ijk1 = analysis_chunk_of_interest_ijk0 + 1
% chunk 3-3-2 is the one
