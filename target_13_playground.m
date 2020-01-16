load target-13-input.mat
%     padded_stack_yxz, ...
%     padded_origin_xyz, ...
%     spacing_xyz, ...
%     origin_xyz, ...
%     shape_xyz, ...
%     parameters

padded_shape_jik = size(padded_stack_yxz) ;
padded_shape_ijk = padded_shape_jik([2 1 3]) ;
padded_far_corner_xyz = padded_origin_xyz + (padded_shape_ijk-1) .* spacing_xyz ;
    % coord of voxel center of highest-indices voxel

xyz_of_interest = [72135.185021 17489.210516 32651.691564] ;
ijk1_of_interest = round( (xyz_of_interest - padded_origin_xyz) ./ spacing_xyz ) + 1 ;
jik1_of_interest = [ijk1_of_interest(2) ijk1_of_interest(1) ijk1_of_interest(3)] ;
serial_voxel_index_of_interest = sub2ind(size(padded_stack_yxz), jik1_of_interest(1), jik1_of_interest(2), jik1_of_interest(3)) ;

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
component_properties_as_table = regionprops3(component_struct, 'Volume', 'Centroid', 'BoundingBox', 'Image', 'VoxelIdxList') ;
component_properties_as_struct = table2struct(component_properties_as_table) ;
voxel_count_from_label = reshape([component_properties_as_struct.Volume], [component_count 1]) ;  % force empty to be right shape
component_centroid_ijk1_from_label = reshape([component_properties_as_struct.Centroid], [3 component_count])' ;  % force empty to be right shape
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

for label = 1 : component_count ,
    %is_in_this_component_yxz = (label_stack_yxz == label) ;
    serial_voxel_indices = serial_voxel_indices_from_label{label} ;
    if ismember(serial_voxel_index_of_interest, serial_voxel_indices) ,
        label_of_interest = label ;
        break ;
    end
end

this_component_properties_as_struct = component_properties_as_struct(label_of_interest) ;

voxel_count = this_component_properties_as_struct.Volume 
centroid_ijk1 = this_component_properties_as_struct.Centroid ;  % this is within the padded substack
centroid_xyz = padded_origin_xyz + spacing_xyz .* (centroid_ijk1-1)  % this is within the full stack

Sigma_regularizer = diag(repmat((0.5*spacing_xyz(1))^2, [3 1])) ;
  % Same approximate effect as addsing some isotropic noise to the voxel positions, so we
  % don't get infinite condition numbers
serial_voxel_indices = this_component_properties_as_struct.VoxelIdxList ;
intensities = padded_stack_yxz(serial_voxel_indices) ;
max_intensity = max(intensities) ;
max_intensity_from_label(label) = max_intensity ;        
if voxel_count < 10 ,
    sqrt_condition_number = 1 
else            
    % Compute first & second moments, and then condition number
    xs = x_grid_yxz(serial_voxel_indices) ;
    ys = y_grid_yxz(serial_voxel_indices) ;
    zs = z_grid_yxz(serial_voxel_indices) ;
    rs = [xs ys zs] ;  % want position vectors in rows
    Sigma_raw = cov(rs) ;
    Sigma = Sigma_raw + Sigma_regularizer ;
    condition_number = cond(Sigma) ;
    sqrt_condition_number = sqrt(condition_number)   % want ratio of SDs, not ratio of variances
end

component_stack = this_component_properties_as_struct.Image ;
component_mip = max(component_stack, [], 3) ;

% Plot the MIP image
f = figure('color', 'w', 'name', 'component-and-centroid') ;
a = axes(f, 'YDir', 'reverse') ;
image(a, 'CData', component_mip, ...
         'XData', [padded_origin_xyz(1) padded_far_corner_xyz(1)], ...
         'YData', [padded_origin_xyz(2) padded_far_corner_xyz(2)], ...
         'CDataMapping', 'scaled') ;         
colormap(gray(256)) ;
axis image
xlabel('x (um)') ;
ylabel('y (um)') ;
line(a, ...
     'XData', centroid_xyz(1), ...
     'YData', centroid_xyz(2), ...
     'Marker', 'o', ...
     'Color', [1 0 0]) ;
 
component_bounding_box = this_component_properties_as_struct.BoundingBox  % this is within the padded substack
component_bounding_box_heckbertian_origin_ijk1 = component_bounding_box(1:3) ;  % within the padded substack
bounding_box_shape_ijk = component_bounding_box(4:6) ;
component_bounding_box_origin_ijk1 = component_bounding_box_heckbertian_origin_ijk1 + 0.5 ;  
  % index of lowest-index voxel in the bounding box, within the padded substack

padded_origin_ijk1 = round( (padded_origin_xyz - origin_xyz) ./ spacing_xyz ) + 1 ;  % lowest-index voxel of the substack, within the full stack

component_stack_origin_ijk1 = padded_origin_ijk1 + component_bounding_box_origin_ijk1 - 1 
  % lowest-index voxel of the component substack, within the full stack

% Compute the distance transform of the component
component_distance_stack = double(bwdist(~component_stack)) ; 

% Plot the MIP image
component_distance_mip = max(component_distance_stack, [], 3) ;
f = figure('color', 'w', 'name', 'component-distance-mip') ;
a = axes(f, 'YDir', 'reverse') ;
image(a, 'CData', component_distance_mip, ...
         'XData', [padded_origin_xyz(1) padded_far_corner_xyz(1)], ...
         'YData', [padded_origin_xyz(2) padded_far_corner_xyz(2)], ...
         'CDataMapping', 'scaled') ;         
colormap(parula(256)) ;
axis image
xlabel('x (um)') ;
ylabel('y (um)') ;
colorbar;


% Plot slice
f = figure('color', 'w', 'name', 'foo') ;
a = axes(f, 'YDir', 'reverse') ;
image(a, 'CData', component_stack(:,:,35), ...
         'XData', [padded_origin_xyz(1) padded_far_corner_xyz(1)], ...
         'YData', [padded_origin_xyz(2) padded_far_corner_xyz(2)], ...
         'CDataMapping', 'scaled') ;         
colormap(gray(256)) ;
axis image
xlabel('x (um)') ;
ylabel('y (um)') ;
colorbar;

% Plot slice
f = figure('color', 'w', 'name', 'foo') ;
a = axes(f, 'YDir', 'reverse') ;
image(a, 'CData', component_distance_stack(:,:,35), ...
         'XData', [padded_origin_xyz(1) padded_far_corner_xyz(1)], ...
         'YData', [padded_origin_xyz(2) padded_far_corner_xyz(2)], ...
         'CDataMapping', 'scaled') ;         
colormap(parula(256)) ;
axis image
xlabel('x (um)') ;
ylabel('y (um)') ;
colorbar;

maximum_distance = max(max(max(component_distance_stack)))
i_max_serial = find(component_distance_stack==maximum_distance,1) ;
[j,i,k] = ind2sub(size(component_distance_stack), i_max_serial) ;
centroidoid_ijk1_within_component = [i j k] 
centroidoid_ijk1 = component_stack_origin_ijk1 + centroidoid_ijk1_within_component - 1  % within 
centroidoid_xyz = origin_xyz + spacing_xyz .* (centroidoid_ijk1-1) 



