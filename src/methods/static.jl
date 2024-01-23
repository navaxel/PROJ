function static_resolution(g::Graph)
    start_time = time()

    n = g.n
    s = g.s
    t = g.t
    S = g.S

    #Model
    model = Model(CPLEX.Optimizer)
    set_optimizer_attribute(model, "CPX_PARAM_SCRIND", 0)

    #Variables
    @variable(model, x[1:n,1:n] >= 0, binary=true)

    #Constraint 0
    for i = 1:n
        for j = 1:n
            if g.d[i,j] == 0
                @constraint(model, x[i,j] == 0)
            end
        end
    end

    #Constraint 1 & 2
    for i = 1:n
        if i != s && i != t
            @constraint(model, sum(x[i,j] for j = 1:n if g.d[i,j] != 0) == sum(x[j,i] for j = 1:n if g.d[i,j] != 0))
        end
    end

    #Constraint 3
    @constraint(model, sum(x[i,s] for i = 1:n if g.d[i,s] != 0) == 0)

    #Constraint 4
    @constraint(model, sum(x[s,j] for j = 1:n if g.d[s,j] != 0) == 1)

    #Constraint 5
    @constraint(model, sum(x[i,t] for i = 1:n if g.d[i,t] != 0) == 1)

    #Constraint 6
    @constraint(model, sum(x[t,j] for j = 1:n if g.d[t,j] != 0) == 0)

    # COnstraint on the weights
    @constraint(model, sum(x[i,j] * (g.p[i] + g.p[j]) for i = 1:n for j = 1:n if g.d[i,j] != 0) + g.p[s] + g.p[t] <= 2 * S)

    #Objective
    @objective(model, Min, sum(g.d[i,j] * x[i,j] for i = 1:n for j = 1:n if g.d[i,j] != 0))

    optimize!(model)

    resolution_time = time() - start_time

    feasibleSolutionFound = primal_status(model) == MOI.FEASIBLE_POINT
    isOptimal = termination_status(model) == MOI.OPTIMAL

    if feasibleSolutionFound
        path = [s]
        i = s
        while i != t
            for j = 1:n
                if value(x[i,j]) == 1
                    push!(path, j)
                    i = j 
                    break
                end
            end
        end
        
        obj_value = JuMP.objective_value(model)

        if save && length(path) > 1
            save_results("Static", g, path, resolution_time)
        end

        return obj_value, path, resolution_time
    end
    
    return nothing, nothing, nothing
end