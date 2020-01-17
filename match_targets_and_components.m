function component_index_from_target_index = ...
        match_targets_and_components(xyz_from_target_index, feature_struct_from_component_index, origin_at_zoom_level_xyz, spacing_at_zoom_level_xyz)
    
    target_count = size(xyz_from_target_index, 1) ;
    component_count = length(feature_struct_from_component_index) ;
    component_index_from_target_index = nan(target_count, 1) ;
    for target_index = 1 : target_count ,
        target_xyz = xyz_from_target_index(target_index, :) ;
        target_ijk1 = round( (target_xyz-origin_at_zoom_level_xyz) ./ spacing_at_zoom_level_xyz) + 1 ;
        for component_index = 1 : component_count ,
            feature_struct = feature_struct_from_component_index(component_index) ;
            target_ijk1_within_component_stack = target_ijk1 - feature_struct.component_stack_origin_ijk1 + 1 ;
            target_jik1_within_component_stack = target_ijk1_within_component_stack([2 1 3]) ;
            component_stack = feature_struct.component_stack ;
            if ismatrix(component_stack) ,
                component_stack_shape_jik = [size(component_stack) 1] ;
            else
                component_stack_shape_jik = size(component_stack) ;
            end
            if all(1<=target_jik1_within_component_stack & target_jik1_within_component_stack<=component_stack_shape_jik) ,
                is_target_in_component = component_stack(target_jik1_within_component_stack(1), ...
                                                         target_jik1_within_component_stack(2), ...
                                                         target_jik1_within_component_stack(3)) ;
                if is_target_in_component ,
                    component_index_from_target_index(target_index) = component_index ;
                    break  % move on to next target
%                 else
%                     % if target is near centroidoid, investigate
%                     centroid_xyz = feature_struct.centroid_xyz ;
%                     distance = sqrt(sum(centroid_xyz-target_xyz).^2) ;
%                     if distance < 20 ,
%                         keyboard
%                     end
                end
            end
        end
    end
end
