function [voxel_count_from_label, ...
          component_centroid_xyz_from_label, ...
          sqrt_condition_number_from_label, ...
          max_intensity_from_label, ...
          is_soma_from_label, ...
          is_within_bounding_box] = ...
        find_somata_in_uint16_stack(padded_stack_yxz, padded_origin_xyz, spacing_xyz, origin_xyz, shape_xyz, parameters)
    % Note that origin_xyz and shape_xyz are the origin and shape of the
    % pre-padding stack.  These determine the bounds of the centroids that
    % will be returned.  Note that origin_xyz is the coordinate of the
    % lowest-index voxel *center* of the unpadded stack.  (Just as
    % padded_origin_xyz is the coordinate of the lowest-index voxel
    % *center* of the padded stack.)
    
    % Unpack parameters            
    intensity_threshold = parameters.intensity_threshold ;
    minimum_volume = parameters.minimum_volume ;
    maximum_volume = parameters.maximum_volume ;
    maximum_sqrt_condition_number = parameters.maximum_sqrt_condition_number ;
    
    volume_per_voxel = prod(spacing_xyz) ;
    minimum_volume_in_voxels = minimum_volume / volume_per_voxel 
    maximum_volume_in_voxels = maximum_volume / volume_per_voxel 
    
    is_bright_enough = logical(padded_stack_yxz>intensity_threshold) ;
    [label_stack_yxz, component_count] = bwlabeln(is_bright_enough) ;

    %stats = regionprops3(label_stack_yxz, {'Volume', 'Centroid',
    %'EigenValues'}) ;  % not sure how to use in presence of anisotropy
        
    % Create grids
    stack_shape_jik = size(padded_stack_yxz) ;
    stack_shape_ijk = stack_shape_jik([2 1 3]) ;   
    x_line = padded_origin_xyz(1) + spacing_xyz(1) * (0:(stack_shape_ijk(1)-1)) ;
    x_grid_yxz = repmat(x_line, [stack_shape_ijk(2) 1 stack_shape_ijk(3)]) ;
    y_line = padded_origin_xyz(2) + spacing_xyz(2) * (0:(stack_shape_ijk(2)-1))' ;
    y_grid_yxz = repmat(y_line, [1 stack_shape_ijk(1) stack_shape_ijk(3)]) ;
    z_line = padded_origin_xyz(3) + spacing_xyz(3) * reshape(0:(stack_shape_ijk(3)-1), [1 1 stack_shape_ijk(3)]) ;
    z_grid_yxz = repmat(z_line, [stack_shape_ijk(2) stack_shape_ijk(1) 1]) ;
    
    voxel_count_from_label = nan(component_count, 1) ;
    component_centroid_xyz_from_label = nan(component_count, 3) ;
    sqrt_condition_number_from_label = nan(component_count, 1) ;
    max_intensity_from_label = nan(component_count, 1) ;
    is_soma_from_label = false(component_count, 1) ;
    Sigma_regularizer = diag(repmat((0.5*spacing_xyz(1))^2, [3 1])) ;
      % Same approximate effect as addsing some isotropic noise to the voxel positions, so we
      % don't get infinite condition numbers
    for label = 1 : component_count ,
        is_in_this_component_yxz = (label_stack_yxz == label) ;
        voxel_count = sum(sum(sum(is_in_this_component_yxz))) ;
        voxel_count_from_label(label) = voxel_count ;
        intensities = padded_stack_yxz(is_in_this_component_yxz) ;
        max_intensity = max(intensities) ;
        max_intensity_from_label(label) = max_intensity ;        
        if voxel_count < minimum_volume_in_voxels/10 ,
            component_centroid_xyz_from_label(label,:) = [nan nan nan] ;
            sqrt_condition_number_from_label(label) = 1 ;
        else            
            % Compute first & second moments, and then condition number
            xs = x_grid_yxz(is_in_this_component_yxz) ;
            ys = y_grid_yxz(is_in_this_component_yxz) ;
            zs = z_grid_yxz(is_in_this_component_yxz) ;
            rs = [xs ys zs] ;  % want position vectors in rows
            centroid = mean(rs,1) ;
            component_centroid_xyz_from_label(label,:) = centroid ;
            if voxel_count == 1 ,
                condition_number = 1 ;
            else
                Sigma_raw = cov(rs) ;
                Sigma = Sigma_raw + Sigma_regularizer ;
                condition_number = cond(Sigma) ;
            end
            sqrt_condition_number = sqrt(condition_number) ;
            sqrt_condition_number_from_label(label) = sqrt_condition_number ;  % want ratio of SDs, not ratio of variances
            if (minimum_volume_in_voxels<voxel_count) && (voxel_count<maximum_volume_in_voxels) && sqrt_condition_number < maximum_sqrt_condition_number ,
                is_soma_from_label(label) = true ;
            end
       end
    end
    
    bounding_box_lower_corner_xyz = origin_xyz - spacing_xyz/2 ;
    bounding_box_upper_corner_xyz = bounding_box_lower_corner_xyz + shape_xyz ;
    is_within_bounding_box = ...
        all(bounding_box_lower_corner_xyz <= component_centroid_xyz_from_label & component_centroid_xyz_from_label < bounding_box_upper_corner_xyz,2) ;
end
