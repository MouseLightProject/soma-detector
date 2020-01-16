function result = ...
        compute_derived_component_features(this_component_properties_as_struct, ...
                                           padded_stack_yxz, ...
                                           padded_stack_origin_xyz, ...
                                           spacing_xyz, ...
                                           origin_xyz, ...
                                           x_grid_yxz, ...
                                           y_grid_yxz, ...
                                           z_grid_yxz)
                                       
    voxel_count = this_component_properties_as_struct.Volume ;
    centroid_ijk1 = this_component_properties_as_struct.Centroid ;  % this is within the padded substack
    centroid_xyz = padded_stack_origin_xyz + spacing_xyz .* (centroid_ijk1-1) ;  % this is within the full stack

    Sigma_regularizer = diag(repmat((0.5*spacing_xyz(1))^2, [3 1])) ;
      % Same approximate effect as addsing some isotropic noise to the voxel positions, so we
      % don't get infinite condition numbers
    serial_voxel_indices = this_component_properties_as_struct.VoxelIdxList ;
    intensities = padded_stack_yxz(serial_voxel_indices) ;
    max_intensity = max(intensities) ;
    if voxel_count < 10 ,
        sqrt_condition_number = 1 ;
    else            
        % Compute first & second moments, and then condition number
        xs = x_grid_yxz(serial_voxel_indices) ;
        ys = y_grid_yxz(serial_voxel_indices) ;
        zs = z_grid_yxz(serial_voxel_indices) ;
        rs = [xs ys zs] ;  % want position vectors in rows
        Sigma_raw = cov(rs) ;
        Sigma = Sigma_raw + Sigma_regularizer ;
        condition_number = cond(Sigma) ;
        sqrt_condition_number = sqrt(condition_number) ;  % want ratio of SDs, not ratio of variances
    end

    component_stack = this_component_properties_as_struct.Image ;
%    component_mip = max(component_stack, [], 3) ;

%     % Plot the MIP image
%     padded_stack_shape_jik = size(padded_stack_yxz) ;
%     padded_stack_shape_ijk = padded_stack_shape_jik([2 1 3]) ;
%     padded_stack_far_corner_xyz = padded_stack_origin_xyz + spacing_xyz .* (padded_stack_shape_ijk-1) ;
%     f = figure('color', 'w', 'name', 'component-and-centroid') ;
%     a = axes(f, 'YDir', 'reverse') ;
%     image(a, 'CData', component_mip, ...
%              'XData', [padded_stack_origin_xyz(1) padded_stack_far_corner_xyz(1)], ...
%              'YData', [padded_stack_origin_xyz(2) padded_stack_far_corner_xyz(2)], ...
%              'CDataMapping', 'scaled') ;         
%     colormap(gray(256)) ;
%     axis image
%     xlabel('x (um)') ;
%     ylabel('y (um)') ;
%     line(a, ...
%          'XData', centroid_xyz(1), ...
%          'YData', centroid_xyz(2), ...
%          'Marker', 'o', ...
%          'Color', [1 0 0]) ;

    component_bounding_box = this_component_properties_as_struct.BoundingBox ;  % this is within the padded substack
    component_bounding_box_heckbertian_origin_ijk1 = component_bounding_box(1:3) ;  % within the padded substack
    %bounding_box_shape_ijk = component_bounding_box(4:6) ;
    component_bounding_box_origin_ijk1 = component_bounding_box_heckbertian_origin_ijk1 + 0.5 ;  
      % index of lowest-index voxel in the bounding box, within the padded substack

    padded_origin_ijk1 = round( (padded_stack_origin_xyz - origin_xyz) ./ spacing_xyz ) + 1 ;  % lowest-index voxel of the substack, within the full stack

    component_stack_origin_ijk1 = padded_origin_ijk1 + component_bounding_box_origin_ijk1 - 1 ;
      % lowest-index voxel of the component substack, within the full stack

    % Compute the distance transform of the component
    component_distance_stack = double(bwdist(~component_stack)) ; 

%     % Plot the MIP image
%     component_distance_mip = max(component_distance_stack, [], 3) ;
%     f = figure('color', 'w', 'name', 'component-distance-mip') ;
%     a = axes(f, 'YDir', 'reverse') ;
%     image(a, 'CData', component_distance_mip, ...
%              'XData', [padded_stack_origin_xyz(1) padded_stack_far_corner_xyz(1)], ...
%              'YData', [padded_stack_origin_xyz(2) padded_stack_far_corner_xyz(2)], ...
%              'CDataMapping', 'scaled') ;         
%     colormap(parula(256)) ;
%     axis image
%     xlabel('x (um)') ;
%     ylabel('y (um)') ;
%     colorbar;
% 
% %     % Plot slice
% %     f = figure('color', 'w', 'name', 'foo') ;
% %     a = axes(f, 'YDir', 'reverse') ;
% %     image(a, 'CData', component_stack(:,:,35), ...
% %              'XData', [padded_stack_origin_xyz(1) padded_stack_far_corner_xyz(1)], ...
% %              'YData', [padded_stack_origin_xyz(2) padded_stack_far_corner_xyz(2)], ...
% %              'CDataMapping', 'scaled') ;         
% %     colormap(gray(256)) ;
% %     axis image
% %     xlabel('x (um)') ;
% %     ylabel('y (um)') ;
% %     colorbar;
% % 
% %     % Plot slice
% %     f = figure('color', 'w', 'name', 'foo') ;
% %     a = axes(f, 'YDir', 'reverse') ;
% %     image(a, 'CData', component_distance_stack(:,:,35), ...
% %              'XData', [padded_stack_origin_xyz(1) padded_stack_far_corner_xyz(1)], ...
% %              'YData', [padded_stack_origin_xyz(2) padded_stack_far_corner_xyz(2)], ...
% %              'CDataMapping', 'scaled') ;         
% %     colormap(parula(256)) ;
% %     axis image
% %     xlabel('x (um)') ;
% %     ylabel('y (um)') ;
% %     colorbar;

    maximum_distance = max(max(max(component_distance_stack))) ;
    i_max_serial = find(component_distance_stack==maximum_distance,1) ;
    [j,i,k] = ind2sub(size(component_distance_stack), i_max_serial) ;
    centroidoid_ijk1_within_component = [i j k] ;
    centroidoid_ijk1 = component_stack_origin_ijk1 + centroidoid_ijk1_within_component - 1 ; % within 
    centroidoid_xyz = origin_xyz + spacing_xyz .* (centroidoid_ijk1-1) ;
    
    % Package things into a result struct
    result = struct('centroid_xyz', centroid_xyz, ...
                    'voxel_count', voxel_count, ...
                    'max_intensity', max_intensity, ...
                    'sqrt_condition_number', sqrt_condition_number, ...
                    'bounding_box_origin_ijk1', component_bounding_box_origin_ijk1, ...
                    'component_stack', component_stack, ...
                    'centroidoid_xyz', centroidoid_xyz) ;
end
