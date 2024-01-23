function static_dijkstra(g::Graph)
    start_time = time()
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

    resolution_time = time() - start_time

    return dist[g.t], path, resolution_time
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

function robust_constraint_dijkstra(g::Graph, random=false::Bool, save=false::Bool)
    dist, path, resolution_time = static_dijkstra(g)
    g_copy = g
    while robust_constraint_eval(g_copy, path) > g.S
        g_copy = deepcopy(g_copy)

        path2 = path[2:end-1]
        i_to_delete = 0
        #Delete one vertice i from path with probability p[i] + ph[i]
        if random
            p_plus_ph = g.p[path2] + g.ph[path2]
            prob = rand() * sum(p_plus_ph)
            cumulative_prob = 0
            for i = 1:length(path2)
                cumulative_prob += p_plus_ph[i]
                if prob <= cumulative_prob
                    i_to_delete = path2[i]
                    break
                end
            end
        #Delete one vertive i from path with the heaviest weight p[i] + ph[i]
        else
            i_to_delete = path[argmax(g.p[path] + g.ph[path])]
        end
        
        for j = 1:g.n
            g_copy.d[i_to_delete,j] = 0
            g_copy.d[j,i_to_delete] = 0
        end
        dist, path, resolution_time2 = static_dijkstra(g_copy)
        resolution_time += resolution_time2
    end

    if save && length(path) > 1
        save_results("RobustDijkstra", g, path, resolution_time)
    end

    return robust_path_eval(g, path), path, resolution_time
end

