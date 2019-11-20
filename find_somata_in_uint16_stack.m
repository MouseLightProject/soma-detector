function result = find_somata_in_uint16_stack(stack_yxz, spacing_xyz, x_grid_yxz, y_grid_yxz, z_grid_yxz)
    intensity_threshold = 0.8 * 2^16 ;
%     minSize = 800; //1500
%     maxSize = minSize * 5;
    
    minimum_volume = 3000 ;  % um^3
    maximum_volume = 5*minimum_volume ;  % um^3
    volume_per_voxel = prod(spacing_xyz) ;
    minimum_volume_in_voxels = minimum_volume / volume_per_voxel 
    maximum_volume_in_voxels = maximum_volume / volume_per_voxel 
    
    is_bright_enough = logical(stack_yxz>intensity_threshold) ;
    [label_stack_yxz, component_count] = bwlabeln(is_bright_enough) ;

    %stats = regionprops3(label_stack_yxz, {'Volume', 'Centroid',
    %'EigenValues'}) ;  % not sure how to use in presence of anisotropy
        
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
            condition_number = cond(Sigma) 
            condition_number_from_label(label) = condition_number ;
            if condition_number < 5 ,
                is_soma(label) = true ;
            end
        end
    end
    
    result = soma_centroid_xyz(is_soma, :) ;
end
