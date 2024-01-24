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

# function robust_constraint_dijkstra(g::Graph, random=false::Bool, save=false::Bool)
#     dist, path, resolution_time = static_dijkstra(g)
#     g_copy = g

#     deleted_vertices = []
#     threshold = div(g.n, 2)

#     i_to_delete = nothing
#     i_to_add = nothing

#     while robust_constraint_eval(g_copy, path) > g.S
#         g_copy = deepcopy(g_copy)

#         #Path without g.s and g.t
#         path2 = path[2:end-1]
        
        
#         #Delete one vertice i from path with probability p[i] + 2*ph[i]
#         # if random
#         #     p_plus_ph = g.p[path2] + g.ph[path2]
#         #     prob = rand() * sum(p_plus_ph)
#         #     cumulative_prob = 0
#         #     for i = 1:length(path2)
#         #         cumulative_prob += p_plus_ph[i]
#         #         if prob <= cumulative_prob
#         #             i_to_delete = path2[i]
#         #             break
#         #         end
#         #     end
#         # end

#         #Delete one vertive i from path with the heaviest weight p[i] + ph[i]
#         sorted_idx = sortperm(g.p[path2] + 2 * g.ph[path2], rev=true)
        
#         idx = 1
#         while idx <= length(sorted_idx)
#             i_to_delete = path2[sorted_idx[idx]]
#             if has_path_with_removed_vertex(g_copy, i_to_delete)
#                 break
#             end
#             idx += 1
#         end
        
#         # No more path -> add one vertex
#         if idx == length(sorted_idx)
#             i_to_add = splice!(deleted_vertices, rand(1:length(deleted_vertices)))
#             i_to_delete = nothing
        
#         # Still path -> delete one vertex and add one if too much vertex deleted
#         else
#             push!(deleted_vertices, i_to_delete)
#             if length(deleted_vertices) > threshold
#                 i_to_add = splice!(deleted_vertices, rand(1:threshold+1))
#             end
#         end
        
#         # i_to_delete = path[argmax(g.p[path] + g.ph[path])]
#         # println(i_to_delete)

#         if !isnothing(i_to_delete)
#             for j = 1:g.n
#                 g_copy.d[i_to_delete,j] = 0
#                 g_copy.d[j,i_to_delete] = 0
#             end
#         end
            
#         if !isnothing(i_to_add)
#             for j = 1:g.n
#                 g_copy.d[i_to_add, j] = g.d[i_to_add, j]
#                 g_copy.d[j, i_to_add] = g.d[j, i_to_add]
#             end
#         end

        
#         dist, path, resolution_time2 = static_dijkstra(g_copy)
        
#         resolution_time += resolution_time2
#     end

#     if save && length(path) > 1
#         save_results("RobustDijkstra", g, path, resolution_time)
#     end

#     return robust_path_eval(g, path), path, resolution_time
# end

function robust_dijkstra(g::Graph, save=false::Bool)
    dist, path, resolution_time = static_dijkstra(g)

    alpha = 1

    g_copy = deepcopy(g)

    #Changer le d_max ?????? (cas robuste ou non)
    d_max = maximum(g.d)
    p_max = maximum(g.p + 2*g.ph)

    norm = p_max / d_max

    p_dist = zeros(Float64, g.n, g.n)
    for i = 1:g.n
        for j = 1:g.n
            if g.d[i,j] != 0
                p_dist[i,j] = (g.p[i] + 2*g.ph[i] + g.p[j] + 2*g.ph[j])/2
            end
        end
    end

    while robust_constraint_eval(g, path) > g.S && alpha >= 0
        alpha -= 0.01
        g_copy.d = alpha * norm * g.d + (1-alpha) * p_dist
        dist, path, resolution_time2 = static_dijkstra(g_copy)

        # println("alpha : ", alpha)
        # println("path : ", path)
        # println("weight : ", robust_constraint_eval(g, path))
        
        resolution_time += resolution_time2
    end


    if save && length(path) > 1
        save_results("RobustDijkstra", g, path, resolution_time)
    end

    return robust_path_eval(g, path), path, resolution_time
end