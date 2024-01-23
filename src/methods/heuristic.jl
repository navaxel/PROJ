struct AntParameters
    m           ::Int
    alpha       ::Float64
    beta        ::Float64
    T           ::Int
    Q           ::Float64
    rho         ::Float64
    tau0        ::Float64
    init_path   ::Union{Nothing,Vector{Int}}
end

function ant_colony(g::Graph, param::AntParameters, save=false::Bool)
    start_time = time()

    n = g.n

    a = param.alpha
    b = param.beta

    best_dist = Inf
    best_path = []

    #tau initialization
    tau = zeros(Float64, n, n)
    for i = 1:n
        for j = 1:n
            if g.d[i,j] != 0
                tau[i,j] = param.tau0
            end
        end
    end
    
    #Path initialization
    if !isnothing(param.init_path)
        for i = 1:length(param.init_path)-1
            tau[param.init_path[i], param.init_path[i+1]] = 200 * param.tau0
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
            println("t = $t, k = $k, path_length = $(length(path)), dist = $dist, path_weight = $path_weight")

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

    resolution_time = time() - start_time

    if save && length(best_path) > 1
        save_results("AntColony", g, best_path, resolution_time)
    end

    return best_dist, best_path, resolution_time
end
