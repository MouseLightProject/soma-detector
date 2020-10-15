sample_date = '2020-09-15' ;
tag = 'v2' ;
output_folder_path = sprintf('/groups/mousebrainmicro/mousebrainmicro/cluster/Reconstructions/%s/soma-predictions-%s', sample_date, tag) ;
rendered_folder_path = sprintf('/nrs/mouselight/SAMPLES/%s', sample_date) ;

do_force_computation = false ;
do_use_bsub = false ;
do_actually_submit = true ;
stdout_file_path_template = fullfile(output_folder_path, 'find-somata-%d-%d-%d.out.txt') ;
bsub_options_template = ['-P mouselight -n8 -eo ' stdout_file_path_template ' -oo ' stdout_file_path_template ' -W 59 -J find-somata'] ;
%bsub_options = '-P mouselight -n8 -eo /dev/null -oo /dev/null -W 59 -J find-somata' ;

foreground_channel_index = 0 ;
background_channel_index = 1 ;
zoom_level = 4 ;  % The zoom level of the tiles we will analyze
pad_depth_in_um = 50 ; % um

% intensity_threshold = 40000 ;
% minimum_volume = 500 ;  % um^3
% maximum_volume = 15000 ;  % um^3
% maximum_sqrt_condition_number = 10 ;

intensity_threshold = 0.85 * 2^16 ;
minimum_volume = 550 ;  % um^3
%maximum_volume = 15000 ;  % um^3
maximum_volume = 25000 ;  % um^3
%maximum_sqrt_condition_number = 10 ;
maximum_sqrt_condition_number = 20 ;

find_somata
