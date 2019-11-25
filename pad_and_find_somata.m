function somata_xyzs = pad_and_find_somata(rendered_folder_path, ...
                                           gfp_channel_index, ...
                                           zoom_level, ...
                                           spacing_at_zoom_level_xyz, ...
                                           origin_at_zoom_level_xyz, ...
                                           substack_origin_ijk1, ...
                                           substack_shape_ijk, ...
                                           pad_depth_in_um, ...
                                           parameters)
    % substack_origin_ijk1 is the origin of the substack within the full stack at this zoom

    % substack_shape_xyz = spacing_at_zoom_level_xyz .* substack_shape_ijk 
    % substack_origin_xyz = origin_at_zoom_level_xyz + spacing_at_zoom_level_xyz .* (substack_origin_ijk1-1) 

    pad_depth_ijk = ceil(pad_depth_in_um ./ spacing_at_zoom_level_xyz) ;  % 3 x 1
    padded_substack_origin_ijk1 = substack_origin_ijk1 - pad_depth_ijk ;
    padded_substack_shape_ijk = substack_shape_ijk + 2*pad_depth_ijk ;
%     padded_substack_shape_xyz = spacing_at_zoom_level_xyz .* padded_substack_shape_ijk 
%     padded_substack_origin_xyz = origin_at_zoom_level_xyz + spacing_at_zoom_level_xyz .* (padded_substack_origin_ijk1-1) 

    padded_substack_yxz = ...
        get_mouselight_rendered_substack(rendered_folder_path, gfp_channel_index, padded_substack_origin_ijk1, padded_substack_shape_ijk, zoom_level) ;  

%     substack_mip = max(padded_substack_yxz,[], 3) ;
% 
%     figure('color', 'w') ;
%     imagesc(substack_mip) ;
%     colormap(gray(256)) ;
%     axis image

    padded_substack_origin_xyz = spacing_at_zoom_level_xyz .* (padded_substack_origin_ijk1-1) + origin_at_zoom_level_xyz ;  
        % coord of voxel center of lowest-indices voxel

    somata_xyzs = find_somata_in_uint16_stack(padded_substack_yxz, padded_substack_origin_xyz, spacing_at_zoom_level_xyz) ;
        % these are in the coordinate system of the full stack
end
