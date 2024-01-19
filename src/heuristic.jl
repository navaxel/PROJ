struct AntParameters
    m       ::Int
    alpha   ::Float64
    beta    ::Float64
    T       ::Int
    Q       ::Float64
    rho     ::Float64
    tau0    ::Float64
end

function ant_colony(g::Graph, param::AntParameters)
    n = g.n

    a = param.alpha
    b = param.beta

    best_dist = Inf
    best_path = []

    #tau initialization
    tau = Matrix{Float64}(undef, n, n)
    for i = 1:n
        for j = 1:n
            if g.d[i,j] != 0
                tau[i,j] = param.tau0
            end
        end
    end

    for t = 1:param.T
        for k = 1:param.m
            i = g.s
            path = [i]
            
            #Visited vertices
            visited = Dict()
            visited[i] = true

            #Add new vertice in path
            while i != g.t

                #Probability distribution
                distrib = Dict()
                for j = 1:n
                    d = g.d[i,j]
                    if d != 0 && !haskey(visited, j)
                        nu = 1 / d
                        distrib[j] = tau[i,j]^a * nu^b
                    end
                end
                
                #Random draw
                if !isempty(distrib)
                    p = rand() * sum(values(distrib))
                    cumulative_prob = 0
                    for (j,prob) in distrib
                        cumulative_prob += prob
                        if p <= cumulative_prob
                            push!(path, j)
                            visited[j] = true
                            i = j
                            break
                        end
                    end
                else
                    break
                end
            end

            #Path not completed
            if i != g.t 
                continue
            end

            #Worst path distance
            dist = robust_path_eval(g, path)

            #Worst path weight 
            path_weight = robust_constraint_eval(g, path)

            #println("t = $t, k = $k, path = $path, dist = $dist, path_weight = $path_weight")

            #Feasible path
            if path_weight <= g.S
                for idx = 1:length(path)-1
                    i = path[idx]
                    j = path[idx+1]
                    tau[i,j] += param.Q / dist
                end

                #Best solution improvement
                if dist < best_dist
                    best_dist = dist
                    best_path = path
                end
            end
        end
        tau .*= (1-param.rho)
    end

    return best_dist, best_path
end


function static_dijkstra(g::Graph)
    n = g.n
    visited = falses(n)
    dist = Inf * ones(Float64, n)
    dist[g.s] = 0
    predecessors = zeros(Int, n)

    for _ in 1:n
        # Trouver le sommet non visité avec la distance minimale
        current_vertex = argmin(dist .+ visited * Inf)
        visited[current_vertex] = true

        # Mettre à jour les distances des voisins non visités et enregistrer les prédécesseurs
        for neighbor in 1:n
            if !visited[neighbor] && g.d[current_vertex, neighbor] != 0
                new_distance = dist[current_vertex] + g.d[current_vertex, neighbor]
                if new_distance < dist[neighbor]
                    dist[neighbor] = new_distance
                    predecessors[neighbor] = current_vertex
                end
            end
        end
    end

    # Reconstruire le chemin à partir des prédécesseurs
    path = get_path(predecessors, g.s, g.t)

    return dist[g.t], path
end

function get_path(predecessors, start, target)
    current_vertex = target
    path = [current_vertex]

    while current_vertex != start
        current_vertex = predecessors[current_vertex]
        pushfirst!(path, current_vertex)
    end

    return path
end