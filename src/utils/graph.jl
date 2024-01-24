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
    d                   ::Matrix{Float64}
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


function has_path_with_removed_vertex(g::Graph, vertex_to_remove::Int)
    visited = falses(g.n)
    stack = [g.s]

    while !isempty(stack)
        current_vertex = pop!(stack)
        if current_vertex == g.t
            return true 
        end

        if !visited[current_vertex]
            visited[current_vertex] = true
            for neigh = 1:g.n
                if g.d[current_vertex, neigh] != 0 && neigh != vertex_to_remove
                    push!(stack, neigh)
                end
            end
        end
    end

    return false 
end