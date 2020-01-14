function [is_match, match_distance] = match_targets_and_guesses(distance_matrix, distance_threshold)
    % distance_matrix is target_count x guess_count
    [target_count, guess_count] = size(distance_matrix) ;
    working_distance_matrix = distance_matrix ;
    maximum_match_count = min(target_count, guess_count) ;
    is_match = false(target_count, guess_count) ;
    match_distance = inf(target_count, guess_count) ;
    for match_index = 1 : maximum_match_count ,
        [minimum_distance_left, target_index, guess_index] = find_min_matrix_element( working_distance_matrix ) ;
        if minimum_distance_left > distance_threshold ,
            break
        end
        is_match(target_index, guess_index) = true ;
        match_distance(target_index, guess_index) = minimum_distance_left ;
        % Setup for next iter
        working_distance_matrix(target_index,:) = inf ;
        working_distance_matrix(:,guess_index) = inf ;
    end    
end
