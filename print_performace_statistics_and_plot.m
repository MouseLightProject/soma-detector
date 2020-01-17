function [matching_candidate_index_from_target_index, matching_target_index_from_candidate_index] = ...
        print_performace_statistics_and_plot(what_are_these, ...
                                             xyz_from_target_index, ...
                                             feature_struct_from_candidate_index, ...
                                             heckbert_origin_xyz, ...
                                             stack_shape_xyz, ...
                                             spacing_at_zoom_level_xyz, ...
                                             substack_mip)
    
    is_there_a_mip = exist('substack_mip', 'var') ;
    if ~is_there_a_mip ,        
        substack_mip = [] ;
    end
    
    % Figure out what targets are in what candidates
    target_count = size(xyz_from_target_index, 1) ;
    candidate_count = length(feature_struct_from_candidate_index) ;
    origin_at_zoom_level_xyz = heckbert_origin_xyz + spacing_at_zoom_level_xyz/2 ;
    matching_candidate_index_from_target_index = ...
        match_targets_and_components(xyz_from_target_index, feature_struct_from_candidate_index, origin_at_zoom_level_xyz, spacing_at_zoom_level_xyz) ;
    matching_target_index_from_candidate_index = invert_partial_map_array(matching_candidate_index_from_target_index, candidate_count) ;
    is_there_a_matched_candidate_from_target_index = isfinite(matching_candidate_index_from_target_index) ;
    is_there_a_matched_target_from_candidate_index = isfinite(matching_target_index_from_candidate_index) ;

    %
    % End of new way
    %

    % Compute hit counts, precision, recall for the candidates
    fprintf('\n\n%s:\n', what_are_these) ;
    target_count
    candidate_count
    hit_count = sum(is_there_a_matched_candidate_from_target_index)
    assert( sum(is_there_a_matched_target_from_candidate_index) == hit_count ) ;
    miss_count = sum(~is_there_a_matched_candidate_from_target_index)
    chase_count = sum(~is_there_a_matched_target_from_candidate_index)

    precision = hit_count / candidate_count 
    recall = hit_count / target_count 


    % Plot the MIP image
    heckbert_far_corner_xyz = heckbert_origin_xyz + stack_shape_xyz ;
    f = figure('color', 'w', 'name', sprintf('targets-and-%s', what_are_these)) ;
    a = axes(f, 'YDir', 'reverse') ;
    if is_there_a_mip ,
        image(a, 'CData', substack_mip, ...
                 'XData', [padded_substack_origin_xyz(1) padded_substack_far_corner_xyz(1)], ...
                 'YData', [padded_substack_origin_xyz(2) padded_substack_far_corner_xyz(2)], ...
                 'CDataMapping', 'scaled') ;         
    end
    xlim([heckbert_origin_xyz(1) heckbert_far_corner_xyz(1)]) ;
    ylim([heckbert_origin_xyz(2) heckbert_far_corner_xyz(2)]) ;
    colormap(gray(256)) ;
    axis image    
    xlabel('x (um)') ;
    ylabel('y (um)') ;


    % plot each target, and each candidate, and draw a line between matches
    hold on ;
    for target_index = 1 : target_count ,
        target_xyz = xyz_from_target_index(target_index,:) ;
        marker_color = fif(is_there_a_matched_candidate_from_target_index(target_index), [0 0.5 1], [1 0 0]) ;    
        plot(target_xyz(1), target_xyz(2), 'Marker', '+', 'Color', marker_color) ;            
        %text(target_xyz(1)+5, target_xyz(2)+5, sprintf('t%d', target_index), 'Color', 0.5*[1 1 1]) ;            
    end
    for candidate_index = 1 : candidate_count ,
        candidate_xyz = feature_struct_from_candidate_index(candidate_index).centroidoid_xyz ;
        marker_color = fif(is_there_a_matched_target_from_candidate_index(candidate_index), [0 0.5 1], [1 0 0]) ;    
        plot(candidate_xyz(1), candidate_xyz(2), 'Marker', 'o', 'MarkerSize', 6, 'Color', marker_color) ;
        %text(candidate_xyz(1)-5, candidate_xyz(2)-5, sprintf('g%d', candidate_index), 'Color', 0.5*[1 1 1]) ;            
    end
    for target_index = 1 : target_count ,
        target_xyz = xyz_from_target_index(target_index,:) ;
        if is_there_a_matched_candidate_from_target_index(target_index) ,
            candidate_index = matching_candidate_index_from_target_index(target_index) ;
            candidate_xyz = feature_struct_from_candidate_index(candidate_index).centroidoid_xyz ;
            plot([target_xyz(1) candidate_xyz(1)], [target_xyz(2) candidate_xyz(2)], 'Color', [0 0.5 1]) ;    
        end
    end
    hold off ; 
    
    fprintf('\n\n') ;
end
