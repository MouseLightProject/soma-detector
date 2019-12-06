function [target_xyzs, is_target_a_soma] = load_targets()
    this_file_path = mfilename('fullpath') ;
    this_folder_path = fileparts(this_file_path) ;
    swc_folder_path = fullfile(this_folder_path, '2019-10-04_adam_AT') ;

    soma_swc_file_name_template = 'Neuron*.swc' ;
    soma_xyzs = load_somata_xyzs_from_swc_files(swc_folder_path, soma_swc_file_name_template) ;

    distractor_swc_file_name_template = '2019-10-04-auto-somata*.swc' ;
    distractor_xyzs = load_somata_xyzs_from_swc_files(swc_folder_path, distractor_swc_file_name_template) ;

    % Put them all into one array
    target_xyzs = vertcat(soma_xyzs, ...
                          distractor_xyzs) ;
    soma_count = size(soma_xyzs, 1) ;
    distractor_count = size(distractor_xyzs, 1) ;
    is_target_a_soma = vertcat(true(soma_count, 1) , ...
                               false(distractor_count, 1)) ;
end