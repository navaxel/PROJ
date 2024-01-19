mutable struct Graph
    name                ::String
    n                   ::Int
    s                   ::Int
    t                   ::Int
    S                   ::Int
    d1                  ::Int
    d2                  ::Int
    p                   ::Vector{Int}
    ph                  ::Vector{Int}
    d                   ::Matrix{Int}
    D                   ::Matrix{Float64}
end

function parse_file(file_name::String)::Graph
    # Define variables to store data
    n = 0
    s = 0
    t = 0
    S = 0
    d1 = 0
    d2 = 0
    p = []
    ph = []
    d = Matrix{Int}(undef, 0, 0)
    D = Matrix{Float64}(undef, 0, 0)

    # Open the file for reading
    file_path = "data/" * file_name
    file = open(file_path, "r")

    # Read each line from the file
    for line in eachline(file)
        # Split the line into tokens based on '=' or whitespace
        tokens = split(line, " = ")

        # Extract data based on the structure
        if tokens[1] == "n"
            n = parse(Int, tokens[2])
            d = zeros(Int, n, n)
            D = zeros(Float64, n, n)
        elseif tokens[1] == "s"
            s = parse(Int, tokens[2])
        elseif tokens[1] == "t"
            t = parse(Int, tokens[2])
        elseif tokens[1] == "S"
            S = parse(Int, tokens[2])
        elseif tokens[1] == "d1"
            d1 = parse(Int, tokens[2])
        elseif tokens[1] == "d2"
            d2 = parse(Int, tokens[2])
        elseif tokens[1] == "p"
            p = parse.(Int, split(tokens[2][2:end-1], ","))
        elseif tokens[1] == "ph"
            ph = parse.(Int, split(tokens[2][2:end-1], ","))
        elseif tokens[1] != "Mat"
            i, j, v, V = split(tokens[1][1:end-1], " ")
            i, j, v = parse.(Int, [i,j,v])
            V = parse(Float64, V)
            d[i,j] = v
            D[i,j] = V
        end
    end

    # Close the file when done
    close(file)

    return Graph(file_name, n, s, t, S, d1, d2, p, ph, d, D)                  
end 


function robust_path_eval(g::Graph, path::Vector{Int})

    path_edges = [(path[i], path[i+1]) for i = 1:length(path) - 1]

    worst_case_lengths = []
    for edge in path_edges
        push!(worst_case_lengths, g.d[edge[1], edge[2]])
    end

    sorted_indexes = sortperm(worst_case_lengths, rev=true)
    sorted_edges = path_edges[sorted_indexes]

    d1 = g.d1
    i = 1
    evaluation = 0

    while d1 >= 0 && i <= length(sorted_edges)
        edge = sorted_edges[i]
        delta = min(g.D[edge[1], edge[2]], d1)
        evaluation += g.d[edge[1], edge[2]] * (1 + delta)
        i += 1
        d1 -= delta
    end

    return evaluation
end


function robust_constraint_eval(g::Graph, path::Vector{Int})

    p_hats = []
    for node in path
        push!(p_hats, g.ph[node])
    end

    sorted_indexes = sortperm(p_hats, rev=true)
    sorted_nodes = path[sorted_indexes]

    d2 = g.d2
    i = 1
    evaluation = 0

    while d2 >= 0 && i <= length(sorted_nodes)
        node = sorted_nodes[i]
        delta = min(2, d2)
        evaluation += g.p[node] + delta * g.ph[i]
        i += 1
        d2 -= delta
    end

    return evaluation
end