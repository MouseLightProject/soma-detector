function [component_centroid_xyz_from_label, ...
          feature_struct_from_label] = ...
        find_candidate_somata_in_uint16_stack(padded_stack_yxz, ...
                                              padded_origin_xyz, ...
                                              spacing_xyz, ...
                                              origin_xyz, ...
                                              shape_xyz, ...
                                              parameters)
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
    
    %volume_per_voxel = prod(spacing_xyz) ;
    %minimum_volume_in_voxels = minimum_volume / volume_per_voxel ;
    %maximum_volume_in_voxels = maximum_volume / volume_per_voxel 
    
    is_bright_enough = logical(padded_stack_yxz > intensity_threshold) ;
    %[label_stack_yxz, component_count] = bwlabeln(is_bright_enough) ;
    component_struct = bwconncomp(is_bright_enough) ;
    serial_voxel_indices_from_label = component_struct.PixelIdxList ;
    component_count = length(serial_voxel_indices_from_label) ;
    component_properties = regionprops3(component_struct, 'Volume', 'Centroid', 'BoundingBox') ;
    %component_properties = regionprops3(component_struct, 'Volume', 'Centroid', 'BoundingBox', 'VoxelList') ;
    voxel_count_from_label = reshape(component_properties.Volume, [component_count 1]) ;  % force empty to be right shape
    component_centroid_ijk1_from_label = reshape(component_properties.Centroid, [component_count 3]) ;  % force empty to be right shape
    %funky_bounding_box_from_label = component_properties.BoundingBox ;
    component_centroid_xyz_from_label = padded_origin_xyz + spacing_xyz .* (component_centroid_ijk1_from_label-1) ;
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
    
    sqrt_condition_number_from_label = nan(component_count, 1) ;
    max_intensity_from_label = nan(component_count, 1) ;
    Sigma_regularizer = diag(repmat((0.5*spacing_xyz(1))^2, [3 1])) ;
      % Same approximate effect as addsing some isotropic noise to the voxel positions, so we
      % don't get infinite condition numbers
    for label = 1 : component_count ,
        %is_in_this_component_yxz = (label_stack_yxz == label) ;
        serial_voxel_indices = serial_voxel_indices_from_label{label} ;
%         if ismember(serial_voxel_index_of_interest, serial_voxel_indices) ,
%             keyboard
%         end
        intensities = padded_stack_yxz(serial_voxel_indices) ;
        max_intensity = max(intensities) ;
        max_intensity_from_label(label) = max_intensity ;        
        if voxel_count_from_label(label) < 10 ,
            sqrt_condition_number_from_label(label) = 1 ;
        else            
            % Compute first & second moments, and then condition number
            xs = x_grid_yxz(serial_voxel_indices) ;
            ys = y_grid_yxz(serial_voxel_indices) ;
            zs = z_grid_yxz(serial_voxel_indices) ;
            rs = [xs ys zs] ;  % want position vectors in rows
            Sigma_raw = cov(rs) ;
            Sigma = Sigma_raw + Sigma_regularizer ;
            condition_number = cond(Sigma) ;
            sqrt_condition_number = sqrt(condition_number) ;
            sqrt_condition_number_from_label(label) = sqrt_condition_number ;  % want ratio of SDs, not ratio of variances
       end
    end
    
    bounding_box_lower_corner_xyz = origin_xyz - spacing_xyz/2 ;
    bounding_box_upper_corner_xyz = bounding_box_lower_corner_xyz + shape_xyz ;
    is_within_bounding_box = ...
        all(bounding_box_lower_corner_xyz <= component_centroid_xyz_from_label & component_centroid_xyz_from_label < bounding_box_upper_corner_xyz,2) ;
    
    % Only keep components with centroids within the bounding box of the
    % unpadded stack
    component_centroid_xyz_from_label = component_centroid_xyz_from_label(is_within_bounding_box, :) ;
    voxel_count_from_label = voxel_count_from_label(is_within_bounding_box, :) ;
    sqrt_condition_number_from_label = sqrt_condition_number_from_label(is_within_bounding_box, :) ;
    max_intensity_from_label = max_intensity_from_label(is_within_bounding_box, :) ;    
    feature_struct_from_label = struct('voxel_count', num2cell(voxel_count_from_label), ...
                                       'sqrt_condition_number', num2cell(sqrt_condition_number_from_label), ...
                                       'max_intensity', num2cell(max_intensity_from_label)) ;
end
