function [soma_xyzs, names] = load_somata_xyzs_from_swc_files(swc_folder_path, soma_swc_file_name_template)
    soma_swc_file_path_template = fullfile(swc_folder_path, soma_swc_file_name_template) ;
    soma_file_names = simple_dir(soma_swc_file_path_template) ;

    soma_count = length(soma_file_names) ;
    soma_xyzs = zeros(soma_count, 3) ;
    names = repmat({''}, [soma_count 1]) ;
    for i = 1 : soma_count ,
        soma_file_name = soma_file_names{i} ;
        soma_file_path = fullfile(swc_folder_path, soma_file_name) ;
        [swc_array, name] = load_swc(soma_file_path) ;
        %node_count = size(swc_array,1) ;
        parent_node_id_from_node_id = swc_array(:,end) ;
        root_node_ids = find(parent_node_id_from_node_id==-1) ;
        root_count = size(root_node_ids) ;
        if root_count == 0 ,
            error('File %s seems to have zero root nodes', soma_file_name) ;
        elseif root_count>1 , 
            error('File %s seems to have more than one root node (it has %d root nodes)', soma_file_name, root_count) ;
        end
        root_node_id = root_node_ids ;
        soma_xyzs(i,:) = swc_array(root_node_id,3:5) ;
        names{i} = name ;
    end
end
