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

    p_hats = g.ph[path]

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

function static_path_eval(g::Graph, path::Vector{Int})
    distance = 0
    for i = 1:length(path)-1
        distance += g.d[path[i], path[i+1]]
    end
    return distance
end

function static_constraint_eval(g::Graph, path::Vector{Int})
    path_weight = 0
    for i = 1:length(path)
        path_weight += g.p[path[i]]
    end
    return path_weight
end

function display_results(g::Graph, path::Vector{Int})
    println(g.name)
    println("Path from ", g.s, " to ", g.t, " : ", path)
    println("Static path distance : ", static_path_eval(g, path))
    println("Static path weight (S = $(g.S)): ", static_constraint_eval(g, path))
    println("Robust path distance : ", robust_path_eval(g, path))
    println("Robust path weight (S = $(g.S)): ", robust_constraint_eval(g, path))
end