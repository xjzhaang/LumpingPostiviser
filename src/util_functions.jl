#################################################################################################################################################
### Sort algorithm to show new variables in reduced echelon form ###
#################################################################################################################################################
#To make intersected matrix into pseudo ref form (upper triangular)
function sort_matrix(matrix::Array{Int})
    new_matrix = zeros(Int, size(matrix)...)
    tuples = [(row_index, findfirst(x -> x != 0, matrix[row_index, :])) for row_index in 1:size(matrix)[1]]
    sort!(tuples, by = x -> x[2])
    for row_index in 1:size(tuples)[1]
        new_matrix[row_index, :] = matrix[tuples[row_index][1], :] 
    end
    return new_matrix
end

#################################################################################################################################################
### Function to find number of nonzero elements in array ###
#################################################################################################################################################
function find_nonzero(list)
    return sum([1 for x in list if x != 0])
end
