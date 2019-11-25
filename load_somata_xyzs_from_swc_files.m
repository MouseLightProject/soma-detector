function soma_xyzs = load_somata_xyzs_from_swc_files(swc_folder_path, soma_swc_file_name_template)
    soma_swc_file_path_template = fullfile(swc_folder_path, soma_swc_file_name_template) ;
    soma_file_names = simple_dir(soma_swc_file_path_template) ;

    soma_count = length(soma_file_names) ;
    soma_xyzs = zeros(soma_count, 3) ;
    for i = 1 : soma_count ,
        soma_file_name = soma_file_names{i} ;
        soma_file_path = fullfile(swc_folder_path, soma_file_name) ;
        swc_array = load_swc(soma_file_path) ;
        node_count = size(swc_array,1) ;
        if node_count == 0 ,
            error('File %s seems to have zero nodes', soma_file_name) ;
        elseif node_count>1 , 
            error('File %s seems to have more than one node (it has %d nodes)', soma_file_name, node_count) ;
        end
        soma_xyzs(i,:) = swc_array(1,3:5) ;
    end
end

    
    