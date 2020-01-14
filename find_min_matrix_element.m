function [A_min, i_min, j_min] = find_min_matrix_element(A)
    [min_per_row, j_min_per_row] = min(A, [], 2) ;
    [A_min, i_min] = min(min_per_row) ;
    j_min = j_min_per_row(i_min) ;
end
