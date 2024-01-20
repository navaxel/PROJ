function display_results(g::Graph, path::Vector{Int}, resolution_time::Float64)
    println(g.name)
    println("Path from ", g.s, " to ", g.t, " : ", path)
    println("Static path distance : ", static_path_eval(g, path))
    println("Static path weight (S = $(g.S)): ", static_constraint_eval(g, path))
    println("Robust path distance : ", robust_path_eval(g, path))
    println("Robust path weight (S = $(g.S)): ", robust_constraint_eval(g, path))
    println("Resolution time : ", round(resolution_time, digits=2))
end

function save_results(method::String, g::Graph, path::Vector{Int}, resolution_time::Float64)
    file = open("results/$(g.name)", "a")
    obj = robust_path_eval(g, path)
    r_t = round(resolution_time, digits=2)
    write(file, "m : $method\np : $path\nt : $r_t\no : $obj\n\n")

    close(file)
end

function get_best_results(g::Graph)
    best_objectives = Dict{String, Float64}()
    best_paths = Dict{String, Vector{Int}}()
    best_times = Dict{String, Float64}()

    file = open("results/$(g.name)", "r")

    current_method = nothing
    current_path = nothing
    current_time = nothing

    for line in eachline(file)
        tokens = split(line, " : ")
        if tokens[1] == "m"
            current_method = tokens[2]

        elseif tokens[1] == "p"
            current_path = parse.(Int, split(tokens[2][2:end-1], ","))
        
        elseif tokens[1] == "t"
            current_time = parse(Float64, tokens[2])

        elseif tokens[1] == "o"
            obj = parse(Float64, tokens[2])
            if haskey(best_objectives, current_method)
                if best_objectives[current_method] > obj
                    best_objectives[current_method] = obj
                    best_paths[current_method] = current_path
                    best_times[current_method] = current_time
                end
            else
                best_objectives[current_method] = obj
                best_paths[current_method] = current_path
                best_times[current_method] = current_time

            end
        end
    end
    close(file)

    return best_objectives, best_paths, best_times
end