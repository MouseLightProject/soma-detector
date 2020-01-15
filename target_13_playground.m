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
component_properties = regionprops3(component_struct, 'Volume', 'Centroid', 'BoundingBox', 'Image') ;
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

for label = 1 : component_count ,
    %is_in_this_component_yxz = (label_stack_yxz == label) ;
    serial_voxel_indices = serial_voxel_indices_from_label{label} ;
    if ismember(serial_voxel_index_of_interest, serial_voxel_indices) ,
        label_of_interest = label ;
        break ;
    end
end

component_image = component_properties.Image{label_of_interest} ;
component_mip = max(component_image, [], 3) ;
component_centroid_xyz = component_centroid_xyz_from_label(label_of_interest, :) ;

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
     'XData', component_centroid_xyz(1), ...
     'YData', component_centroid_xyz(2), ...
     'Marker', 'o', ...
     'Color', [1 0 0]) ;
 
component_bounding_box = component_properties.BoundingBox(label_of_interest,:) ;
component_bounding_box_heckbertian_origin_ijk1 = component_bounding_box(1:3) ;
bounding_box_shape_ijk = component_bounding_box(4:6) ;
component_bouncing_box_origin_ijk1 = component_bounding_box_heckbertian_origin_ijk1 + 0.5 ;  % index of lowest-index voxel in the bounding box


 
     
     