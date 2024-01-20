function cutting_planes_resolution(g::Graph, save=false::Bool)
    start_time = time()
    n = g.n
    s = g.s
    t = g.t

    # Initial cutting planes
    u_1_prime = []
    u_2_prime = []

    separation_verified = false

    master_x = zeros(Int, n, n)
    master_obj = 0

    while !separation_verified
        master_obj, master_x, master_eta = master_problem(g, u_1_prime, u_2_prime)

        obj1, delta1 = sp_A(g, master_x, master_eta)
        if obj1 > 0
            push!(u_1_prime, delta1)
        end

        obj2, delta2 = sp_A(g, master_x, master_eta)
        if obj2 > 0
            push!(u_2_prime, delta2)
        end

        if obj1 <= 0 && obj2 <= 0
            separation_verified = true
        end
    end 

    # Recover optimal path
    path = [s]
    i = s
    while i != t
        for j = 1:n
            if master_x[i,j] == 1
                push!(path, j)
                i = j 
                break
            end
        end
    end
    obj_value = master_obj

    resolution_time = time() - start_time

    if save
        save_results("CuttingPlanes", g, path, resolution_time)
    end
    
    return obj_value, path, resolution_time
end


function master_problem(g::Graph, u_1_prime::Vector{Any}, u_2_prime::Vector{Any})

    n = g.n
    s = g.s
    t = g.t
    S = g.S

    #Model
    master_model = Model(CPLEX.Optimizer)
    set_optimizer_attribute(master_model, "CPX_PARAM_SCRIND", 0)

    #Variables
    @variable(master_model, x[1:n,1:n] >= 0, binary=true)
    @variable(master_model, eta >= 0)

    #Constraint 0
    for i = 1:n
        for j = 1:n
            if g.d[i,j] == 0
                @constraint(master_model, x[i,j] == 0)
            end
        end
    end

    #Constraint 1 & 2
    for i = 1:n
        if i != s && i != t
            @constraint(master_model, sum(x[i,j] for j = 1:n if g.d[i,j] != 0) == sum(x[j,i] for j = 1:n if g.d[i,j] != 0))
        end
    end

    #Constraint 3
    @constraint(master_model, sum(x[i,s] for i = 1:n if g.d[i,s] != 0) == 0)

    #Constraint 4
    @constraint(master_model, sum(x[s,j] for j = 1:n if g.d[s,j] != 0) == 1)

    #Constraint 5
    @constraint(master_model, sum(x[i,t] for i = 1:n if g.d[i,t] != 0) == 1)

    #Constraint 6
    @constraint(master_model, sum(x[t,j] for j = 1:n if g.d[t,j] != 0) == 0)

    #Constraint on U_1^{*}
    for delta_1 in u_1_prime
        @constraint(master_model, sum(g.d[i,j]*x[i,j]*delta_1[i,j] for i = 1:n for j = 1:n if g.d[i,j] != 0) <= eta)
    end

    #Constraint on U_2^{*}
    for delta_2 in u_2_prime
        @constraint(master_model, sum(x[i,j] * (g.p[i] + g.p[j] + delta_2[i] * g.ph[i] + delta_2[j] * g.ph[j]) for i = 1:n for j = 1:n if g.d[i,j] != 0) + 
        g.p[s] + delta_2[s] * g.ph[s] + g.p[t] + delta_2[t] * g.ph[t] <= 2 * S)
    end

    #Objective
    @objective(master_model, Min, sum(g.d[i,j]*x[i,j] for i = 1:n for j = 1:n if g.d[i,j] != 0) + eta)

    optimize!(master_model)

    feasibleSolutionFound = primal_status(master_model) == MOI.FEASIBLE_POINT
    isOptimal = termination_status(master_model) == MOI.OPTIMAL

    if feasibleSolutionFound
        x_opt = zeros(Int, n, n)
        for i = 1:n
            for j = 1:n
                x_opt[i,j] = value(x[i,j])
            end
        end

        eta_opt = value(eta)
        master_obj = JuMP.objective_value(master_model)
    
        return master_obj, x_opt, eta_opt
    end
    
    return nothing, nothing, nothing
end


function sp_A(g::Graph, x::Matrix{Int}, eta::Float64)
    n = g.n
    s = g.s
    t = g.t

    model = Model(CPLEX.Optimizer)
    set_optimizer_attribute(model, "CPX_PARAM_SCRIND", 0)


    #Variables
    @variable(model, delta[1:n,1:n] >= 0)

    #Constraint on the sum of delta's coordinates
    @constraint(model, sum(delta[i,j] for i = 1:n for j = 1:n if g.d[i,j] != 0) <= g.d1)

    #Constraint on each delta[i,j] 
    for i = 1:n
        for j = 1:n
            if g.d[i,j] != 0
                @constraint(model, delta[i,j] <= g.D[i,j])
            else
                @constraint(model, delta[i,j] == 0)
            end
        end
    end

    #Objective
    @objective(model, Max, sum(g.d[i,j]*x[i,j]*delta[i,j] for i = 1:n for j = 1:n if g.d[i,j] != 0) - eta)

    optimize!(model)

    feasibleSolutionFound = primal_status(model) == MOI.FEASIBLE_POINT
    isOptimal = termination_status(model) == MOI.OPTIMAL

    if feasibleSolutionFound
        delta_opt = zeros(Float64, n, n)
        for i = 1:n
            for j = 1:n
                delta_opt[i,j] = value(delta[i,j])
            end
        end
        
        obj_value = JuMP.objective_value(model)

        return obj_value, delta_opt
    end
    
    return nothing, nothing
end


function sp_B(g::Graph, x::Matrix{Int}, eta::Float64)
    n = g.n
    s = g.s
    t = g.t
    S = g.S

    model = Model(CPLEX.Optimizer)
    set_optimizer_attribute(model, "CPX_PARAM_SCRIND", 0)

    #Variables
    @variable(model, delta[1:n] >= 0)

    #Constraint on the sum of delta's coordinates
    @constraint(model, sum(delta[i] for i = 1:n) <= g.d2)

    #Constraint on each delta[i] 
    for i = 1:n
        @constraint(model, delta[i] <= 2)
    end

    #Objective
    @objective(model, Max, sum(x[i,j] * (g.p[i] + g.p[j] + delta[i] * g.ph[i] + delta[j] * g.ph[j]) for i = 1:n for j = 1:n if g.d[i,j] != 0) + 
    g.p[s] + delta[s] * g.ph[s] + g.p[t] + delta[t] * g.ph[t] - 2 * S)

    optimize!(model)

    feasibleSolutionFound = primal_status(model) == MOI.FEASIBLE_POINT
    isOptimal = termination_status(model) == MOI.OPTIMAL

    if feasibleSolutionFound
        delta_opt = zeros(Float64, n)
        for i = 1:n
            delta_opt[i] = value(delta[i])
        end
        
        obj_value = JuMP.objective_value(model)

        return obj_value, delta_opt
    end
    
    return nothing, nothing
end

