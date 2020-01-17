sample_date = '2019-10-04' ;
rendered_folder_path = sprintf('/nrs/mouselight/SAMPLES/%s', sample_date) ;
roi_center_xyz = [73846.8, 18235.4, 36316.3] ;  % um
roi_corner_xyz = [74884.4, 19198.9, 36316.3] ; % um 
pad_depth_in_um = 50 ; % um
zoom_level = 4 ;
gfp_channel_index = 0 ;
background_channel_index = 1 ;

intensity_threshold = 0.75 * 2^16 ;
minimum_volume = 500 ;  % um^3
%maximum_volume = 15000 ;  % um^3
maximum_volume = 25000 ;  % um^3
maximum_sqrt_condition_number = 10 ;

debug_on_substack
