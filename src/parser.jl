using Nemo
include("poly_eval.jl")

function parse_matrix(txt)
    f = open(txt)
    lines = readlines(f)
    #first we create a zeros matrix of the correct dimension.
    matrix = zeros(length(lines), length(split(lines[1])))
    matrix = Array{Rational{Int64}}(matrix)

    for line_index in 1:length(lines)
        vector = split(lines[line_index])
        for v_index in 1:length(vector)
            max_deminominator = 1
            # first we check if the element is a fraction, if so we find its numerator and denominator and create a rational to represent it
            if occursin("/", vector[v_index])
                numerator = parse(Int,split(vector[v_index], "/")[1])
                denominator = parse(Int,split(vector[v_index], "/")[2])
                if denominator > max_deminominator
                    max_deminominator = denominator
                end
                rational = numerator // denominator
                matrix[line_index, v_index] = rational
            #if the element is an int, we parse the string to int and add it to the matrix
            else
                matrix[line_index, v_index] = parse(Int64, vector[v_index])
            end
        end
    end
    matrix = rational_to_int(matrix)
    return matrix
end



function parse_polynomial(txt)
    # parse the text file
    f = open(txt)
    lines = readlines(f)
    lines = julia_exponent(lines)

    variables_str = Array{String}([])
    for index in 0:length(lines) - 1
        push!(variables_str, "y$index")
    end
    
    #Create the PolynomiaL ring
    R, v = PolynomialRing(Nemo.QQ, variables_str)
    S = MatrixSpace(R, length(lines), length(lines))

    #Create a dictionary sending symbol to v[i]
    expr_dict = Dict(Symbol("y$num") => v[num + 1] for num in 0:(length(v) - 1))
    
    poly_system = Array{fmpq_mpoly}([])
    for line_index in 1:length(lines)
        if myeval(Meta.parse(lines[line_index]), expr_dict) == fmpq(0)
            push!(poly_system, R(0))
        else
            push!(poly_system, myeval(Meta.parse(lines[line_index]), expr_dict))
        end
    end
    return poly_system
end

#Transforms Array{Rational{Int}} matrix to Array{Int}
function rational_to_int(matrix::Array{Rational{Int}})
    for row in eachrow(matrix)
        denom = 1
        lcm_ = 1
        for i in row
            if typeof(i) == Rational{Int64}
                denom = denominator(i)
                lcm_ = lcm(lcm_, denom)
            end
        end
        row .= row * lcm_
    end
    matrix = Array{Int}(matrix)
    return matrix
end

function julia_exponent(lines)
    for index in 1:size(lines)[1]
        lines[index] = replace(lines[index], "**" => "^")
    end
    return lines
end


function parse_varnames(txt)
    f = open(txt)
    lines = readlines(f)
    
    variables_dict = Dict()
    for l in 0:length(lines) - 1
        variables_dict[l] = lines[l + 1]
    end
    return variables_dict
end
