sample_date = '2019-10-04' ;
rendered_folder_path = sprintf('/nrs/mouselight/SAMPLES/%s', sample_date) ;
example_soma_location_xyz = [74532.9, 17809.5, 36575.8] ;  % um
pad_depth_in_um = 50 ; % um
zoom_level = 4 ;

%[shape_xyz, origin_xyz, spacing_xyz] = load_sample_shape_origin_and_spacing(rendered_folder_path) ;
parameter_file_path = fullfile(rendered_folder_path, 'calculated_parameters.jl') ;
parameters = read_renderer_calculated_parameters_file(parameter_file_path) ;
max_zoom_level = parameters.level_step_count ;
chunk_shape_ijk = parameters.leaf_shape ;  % xyz order, same at all zoom levels, just the chunk count changes
spacing_at_max_zoom_xyz = parameters.spacing ;
origin_at_max_zoom_xyz = parameters.origin ;
%shape_xyz = 2^max_zoom_level * chunk_shape_ijk ;

spacing_at_zoom_level_0_xyz = 2^max_zoom_level * spacing_at_max_zoom_xyz ;
spacing_at_zoom_level_xyz = spacing_at_zoom_level_0_xyz ./ (2^zoom_level) ;

heckbert_origin_xyz = origin_at_max_zoom_xyz - spacing_at_max_zoom_xyz/2 ;  % this origin does not change with the zoom level
origin_at_zoom_level_xyz = heckbert_origin_xyz + spacing_at_zoom_level_xyz/2 ;



substack_half_shape_xyz = 100*[1 1 1] ; % um
substack_desired_shape_xyz = 2 * substack_half_shape_xyz 
substack_desired_origin_xyz = example_soma_location_xyz - substack_half_shape_xyz 
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



gfp_channel_index = 0 ;
padded_substack_yxz = ...
    get_mouselight_rendered_substack(rendered_folder_path, gfp_channel_index, padded_substack_origin_ijk1, padded_substack_shape_ijk, zoom_level) ;  

substack_mip = max(padded_substack_yxz,[], 3) ;

figure('color', 'w') ;
imagesc(substack_mip) ;
colormap(gray(256)) ;
axis image

padded_substack_origin_xyz = spacing_at_zoom_level_xyz .* (padded_substack_origin_ijk1-1) + origin_at_zoom_level_xyz ;  
    % coord of voxel center of lowest-indices voxel

x_line = padded_substack_origin_xyz(1) + spacing_at_zoom_level_xyz(1) * (0:(padded_substack_shape_ijk(1)-1)) ;
x_grid_yxz = repmat(x_line, [padded_substack_shape_ijk(2) 1 padded_substack_shape_ijk(3)]) ;

y_line = padded_substack_origin_xyz(2) + spacing_at_zoom_level_xyz(2) * (0:(padded_substack_shape_ijk(2)-1))' ;
y_grid_yxz = repmat(y_line, [1 padded_substack_shape_ijk(1) padded_substack_shape_ijk(3)]) ;

z_line = padded_substack_origin_xyz(3) + spacing_at_zoom_level_xyz(3) * reshape(0:(padded_substack_shape_ijk(3)-1), [1 1 padded_substack_shape_ijk(3)]) ;
z_grid_yxz = repmat(z_line, [padded_substack_shape_ijk(2) padded_substack_shape_ijk(1) 1]) ;

somata_xyzs = find_somata_in_uint16_stack(padded_substack_yxz, spacing_at_zoom_level_xyz, x_grid_yxz, y_grid_yxz, z_grid_yxz)
    % these are in the coordinate system of the full stack

hold on ;
somata_count = size(somata_xyzs,1) ;
for i = 1 : somata_count ,
    soma_xyz = somata_xyzs(i,:) ;
    % have to translate to image coords in the mip
    soma_ijk1_fractional = (soma_xyz - padded_substack_origin_xyz) ./ spacing_at_zoom_level_xyz + 1 ;
    plot(soma_ijk1_fractional(1), soma_ijk1_fractional(2), 'r+') ;    
end
hold off ; 
