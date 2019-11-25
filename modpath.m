this_file_path = mfilename('fullpath') ;
this_folder_path = fileparts(this_file_path) ;                

addpath('mouselight_toolbox');

tmt_path = fullfile(this_folder_path, 'tmt') ;
cd(tmt_path) ;
modpath ;
cd(this_folder_path) ;

