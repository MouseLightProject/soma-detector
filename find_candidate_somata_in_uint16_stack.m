function component_features_from_label = ...
        find_candidate_somata_in_uint16_stack(padded_substack_yxz, ...
                                              padded_substack_origin_xyz, ...
                                              origin_at_zoom_level_xyz, ...
                                              spacing_at_zoom_level_xyz, ...
                                              substack_origin_xyz, ...
                                              substack_shape_xyz, ...
                                              parameters)
                                          
%     find_candidate_somata_in_uint16_stack(...
%         padded_substack_yxz, ...
%         padded_substack_origin_xyz, ...
%         spacing_at_zoom_level_xyz, ...
%         substack_origin_xyz, ...
%         substack_shape_xyz, ...
%         parameters) ;
                                          
                                          
    % Note that origin_xyz and shape_xyz are the origin and shape of the
    % pre-padding stack.  These determine the bounds of the centroids that
    % will be returned.  Note that origin_xyz is the coordinate of the
    % lowest-index voxel *center* of the unpadded stack.  (Just as
    % padded_origin_xyz is the coordinate of the lowest-index voxel
    % *center* of the padded stack.)
    
%     xyz_of_interest = [72135.185021 17489.210516 32651.691564] ;
%     ijk1_of_interest = round( (xyz_of_interest - padded_origin_xyz) ./ spacing_xyz ) + 1 ;
%     jik1_of_interest = [ijk1_of_interest(2) ijk1_of_interest(1) ijk1_of_interest(3)] ;
%     serial_voxel_index_of_interest = sub2ind(size(padded_stack_yxz), jik1_of_interest(1), jik1_of_interest(2), jik1_of_interest(3)) ;
    
    % Unpack parameters            
    intensity_threshold = parameters.intensity_threshold ;
    %minimum_volume = parameters.minimum_volume ;
    %maximum_volume = parameters.maximum_volume ;
    %maximum_sqrt_condition_number = parameters.maximum_sqrt_condition_number ;
    
    % Threshold stack and extract components    
    is_bright_enough = logical(padded_substack_yxz > intensity_threshold) ;
    component_struct = bwconncomp(is_bright_enough) ;
    
    % Compute various useful properties of each componenent
    component_properties_as_table = regionprops3(component_struct, 'Volume', 'Centroid', 'BoundingBox', 'Image', 'VoxelIdxList') ;
    component_properties_as_struct = table2struct(component_properties_as_table) ;
        
    % Create grids, used for calcualating covariance matrix of each
    % component
    stack_shape_jik = size(padded_substack_yxz) ;
    stack_shape_ijk = stack_shape_jik([2 1 3]) ;   
    x_line = padded_substack_origin_xyz(1) + spacing_at_zoom_level_xyz(1) * (0:(stack_shape_ijk(1)-1)) ;
    x_grid_yxz = repmat(x_line, [stack_shape_ijk(2) 1 stack_shape_ijk(3)]) ;
    y_line = padded_substack_origin_xyz(2) + spacing_at_zoom_level_xyz(2) * (0:(stack_shape_ijk(2)-1))' ;
    y_grid_yxz = repmat(y_line, [1 stack_shape_ijk(1) stack_shape_ijk(3)]) ;
    z_line = padded_substack_origin_xyz(3) + spacing_at_zoom_level_xyz(3) * reshape(0:(stack_shape_ijk(3)-1), [1 1 stack_shape_ijk(3)]) ;
    z_grid_yxz = repmat(z_line, [stack_shape_ijk(2) stack_shape_ijk(1) 1]) ;

    % Compute the component features we really want
    if isempty(component_properties_as_struct) ,
        raw_component_features_from_label = compute_derived_component_features(component_properties_as_struct) ;
    else
        raw_component_features_from_label = ...
            arrayfun(...
                @(s)(compute_derived_component_features(s, padded_substack_yxz, padded_substack_origin_xyz, origin_at_zoom_level_xyz, spacing_at_zoom_level_xyz, x_grid_yxz, y_grid_yxz, z_grid_yxz)), ...
                component_properties_as_struct) ;
    end

    % Extract the centroidoid for each component
    component_count = length(raw_component_features_from_label) ;
    component_centroidoid_xyz_from_label = ...
        reshape([raw_component_features_from_label.centroidoid_xyz], [3 component_count])' ;  % force empty to be right shape

    % Compute the bounding box for the unpadded stack
    bounding_box_lower_corner_xyz = substack_origin_xyz - spacing_at_zoom_level_xyz/2 ;
    bounding_box_upper_corner_xyz = bounding_box_lower_corner_xyz + substack_shape_xyz ;
    
    % Only keep components with centroids within the bounding box of the
    % unpadded stack
    is_within_bounding_box = ...
        all(bounding_box_lower_corner_xyz <= component_centroidoid_xyz_from_label & component_centroidoid_xyz_from_label < bounding_box_upper_corner_xyz,2) ;
    component_features_from_label = raw_component_features_from_label(is_within_bounding_box) ;
end
