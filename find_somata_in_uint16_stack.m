function result = find_somata_in_uint16_stack(stack_yxz, origin_xyz, spacing_xyz, parameters)    
    % Unpack parameters            
    intensity_threshold = parameters.intensity_threshold ;
    minimum_volume = parameters.minimum_volume ;
    maximum_volume = parameters.maximum_volume ;
    maximum_condition_number = parameters.maximum_condition_number ;
    
    volume_per_voxel = prod(spacing_xyz) ;
    minimum_volume_in_voxels = minimum_volume / volume_per_voxel ;
    maximum_volume_in_voxels = maximum_volume / volume_per_voxel ;
    
    is_bright_enough = logical(stack_yxz>intensity_threshold) ;
    [label_stack_yxz, component_count] = bwlabeln(is_bright_enough) ;

    %stats = regionprops3(label_stack_yxz, {'Volume', 'Centroid',
    %'EigenValues'}) ;  % not sure how to use in presence of anisotropy
        
    % Create grids
    stack_shape_jik = size(stack_yxz) ;
    stack_shape_ijk = stack_shape_jik([2 1 3]) ;   
    x_line = origin_xyz(1) + spacing_xyz(1) * (0:(stack_shape_ijk(1)-1)) ;
    x_grid_yxz = repmat(x_line, [stack_shape_ijk(2) 1 stack_shape_ijk(3)]) ;
    y_line = origin_xyz(2) + spacing_xyz(2) * (0:(stack_shape_ijk(2)-1))' ;
    y_grid_yxz = repmat(y_line, [1 stack_shape_ijk(1) stack_shape_ijk(3)]) ;
    z_line = origin_xyz(3) + spacing_xyz(3) * reshape(0:(stack_shape_ijk(3)-1), [1 1 stack_shape_ijk(3)]) ;
    z_grid_yxz = repmat(z_line, [stack_shape_ijk(2) stack_shape_ijk(1) 1]) ;
    
    voxel_count_from_label = nan(1, component_count) ;
    condition_number_from_label = nan(1, component_count) ;
    soma_centroid_xyz = nan(component_count, 3) ;
    is_soma = false(1, component_count) ;
    for label = 1 : component_count ,
        is_in_this_component_yxz = (label_stack_yxz == label) ;
        voxel_count = sum(sum(sum(is_in_this_component_yxz))) ;
        voxel_count_from_label(label) = voxel_count ;
        if (minimum_volume_in_voxels<voxel_count) && (voxel_count<maximum_volume_in_voxels) ,
            % Compute first & second moments, and then condition number
            xs = x_grid_yxz(is_in_this_component_yxz) ;
            ys = y_grid_yxz(is_in_this_component_yxz) ;
            zs = z_grid_yxz(is_in_this_component_yxz) ;
            rs = [xs ys zs] ;  % want position vectors in rows
            centroid = mean(rs) ;
            soma_centroid_xyz(label,:) = centroid ;
            Sigma = cov(rs) ;
%             [~, Lambda] = eig(Sigma) ;
%             lambda = diag(Lambda) ;
%             k_maybe = lambda(end)/lambda(1) 
            condition_number = cond(Sigma) ;
            condition_number_from_label(label) = condition_number ;
            if condition_number < maximum_condition_number ,
                is_soma(label) = true ;
            end
        end
    end
    
    result = soma_centroid_xyz(is_soma, :) ;
end
