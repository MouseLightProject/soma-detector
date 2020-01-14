function soma_xyzs = load_soma_targets()
    this_file_path = mfilename('fullpath') ;
    this_folder_path = fileparts(this_file_path) ;
    swc_folder_path = fullfile(this_folder_path, '2019-10-04-somata-gt-round-3') ;

    soma_swc_file_name_template = 'Neuron*.swc' ;
    soma_xyzs = load_somata_xyzs_from_swc_files(swc_folder_path, soma_swc_file_name_template) ;
end
