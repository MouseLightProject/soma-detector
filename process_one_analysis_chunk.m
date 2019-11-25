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
if ~exist(somata_folder_path, 'file') ,
    mkdir(somata_folder_path) ;
end
% Just analyze a chunk with a good number of somata
for chunk_i = 3 ,
    for chunk_j = 3 ,
        for chunk_k = 2 ,
            somata_mat_file_name = sprintf('somata-for-chunk-%d-%d-%d.mat', chunk_i, chunk_j, chunk_k) ;
            somata_mat_file_path = fullfile(somata_folder_path, 'chunks', somata_mat_file_name) ;
            chunk_offset_within_chunks_ijk1 = [chunk_i chunk_j chunk_k] 
            chunk_offset_within_stack_ijk1 = (chunk_offset_within_chunks_ijk1-1) .* analysis_chunk_shape_ijk + 1 ;
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
end

