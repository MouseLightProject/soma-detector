function result = count_voxels_per_component(label_stack, component_count)
    result = zeros(1, component_count) ;
    for component_index = 1 : component_count ,
        is_in_this_component = (label_stack == component_index) ;
        voxel_count = sum(sum(sum(is_in_this_component))) ;
        result(compoent_index) = voxel_count ;
    end    
end
