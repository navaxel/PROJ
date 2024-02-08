function warm_start(g::Graph)
    x_init = zeros(Int, g.n, g.n)

    obj_value, path, resolution_time = robust_dijkstra(g)

    for i = 1:length(path)-1
        x_init[path[i],path[i+1]] = 1
    end

    return x_init
end


function branch_and_cut_resolution(g::Graph, save=false::Bool, time_limit=nothing::Union{Nothing,Int}, warmstart=nothing::Union{Nothing,Bool})
    start_time = time()

    n = g.n
    s = g.s
    t = g.t
    S = g.S

    #Model
    model = Model(CPLEX.Optimizer)
    #set_optimizer_attribute(model, "CPX_PARAM_SCRIND", 0)
    if !isnothing(time_limit)
        set_time_limit_sec(model, time_limit)
    end

    #Variables
    @variable(model, x[1:n,1:n] >= 0, binary=true)
    @variable(model, eta >= 0)

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

    # Adding static constraint
    @constraint(model, sum(x[i,j] * (g.p[i] + g.p[j]) for i = 1:n for j = 1:n if g.d[i,j] != 0) + g.p[s] + g.p[t] <= 2 * g.S)
            
    #Objective
    @objective(model, Min, sum(g.d[i,j]*x[i,j] for i = 1:n for j = 1:n if g.d[i,j] != 0) + eta + sum(x[i,j] * 200 * i for i = 1:n for j = 1:n if g.d[i,j] != 0))
    #@objective(model, Min, sum(g.d[i,j]*x[i,j] for i = 1:n for j = 1:n if g.d[i,j] != 0) + eta)

    # Callbacks
    function callback_function(cb_data)

        status = callback_node_status(cb_data, model)
        if !(status == MOI.CALLBACK_NODE_STATUS_FRACTIONAL)
            
            x_val = zeros(Int, n, n)
            for i = 1:n
                for j = 1:n
                    if callback_value(cb_data, x[i,j]) >= 1 - 1e-6
                        x_val[i,j] = 1
                    end
                end
            end
            eta_val = callback_value(cb_data, eta)

            obj1, delta1 = sp_A(g, x_val, eta_val)
            obj2, delta2 = sp_B(g, x_val, eta_val)
            
            if obj1 >= 1e-6
                cstr1 = @build_constraint(sum(g.d[i,j]*x[i,j]*delta1[i,j] for i = 1:n for j = 1:n if g.d[i,j] != 0) <= eta)
                MOI.submit(model, MOI.LazyConstraint(cb_data), cstr1)    
            end

            if obj2 >= 1e-6
                cstr2 = @build_constraint(sum(x[i,j] * (g.p[i] + g.p[j] + delta2[i] * g.ph[i] + delta2[j] * g.ph[j]) for i = 1:n for j = 1:n if g.d[i,j] != 0) + 
                g.p[s] + delta2[s] * g.ph[s] + g.p[t] + delta2[t] * g.ph[t] <= 2 * g.S)
                MOI.submit(model, MOI.LazyConstraint(cb_data), cstr2)
            end

            # if obj1 <= 1e-6 && obj2 <= 1e-6
            #     println("should terminate")
            # end 
        end
    end 

    set_attribute(model, MOI.LazyConstraintCallback(), callback_function)

    # Warm start
    if warmstart
        x_init = warm_start(g)
        set_start_value.(x, x_init)
    end

    optimize!(model)

    feasibleSolutionFound = primal_status(model) == MOI.FEASIBLE_POINT
    isOptimal = termination_status(model) == MOI.OPTIMAL

    if feasibleSolutionFound
        x_opt = zeros(Int, n, n)
        for i = 1:n
            for j = 1:n
                if value(x[i,j]) >= 1 - 1e-6
                    x_opt[i,j] = 1
                end
            end
        end

        eta_opt = value(eta)
        opt_obj = JuMP.objective_value(model)

        # Recover optimal path
        path = [s]
        i = s
        while i != t
            for j = 1:n
                if x_opt[i,j] == 1
                    push!(path, j)
                    i = j 
                    break
                end
            end
        end
        obj_value = opt_obj

        resolution_time = time() - start_time

        if save && length(path) > 1 && robust_constraint_eval(g, path) <= g.S
            save_results("BranchCut", g, path, resolution_time)
        end
    
        return obj_value, path, resolution_time
    end
    return nothing, nothing, nothing
end
