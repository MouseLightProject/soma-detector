function [soma_xyzs, names] = load_traceable_soma_targets_from_tracers()
    this_file_path = mfilename('fullpath') ;
    this_folder_path = fileparts(this_file_path) ;
    swc_folder_path = fullfile(this_folder_path, 'ground-truth-traceable-somata-for-2019-10-04-from-monet-take-2') ;

    soma_swc_file_name_template = '*.swc' ;
    [soma_xyzs, names] = load_somata_xyzs_from_swc_files(swc_folder_path, soma_swc_file_name_template) ;
end
