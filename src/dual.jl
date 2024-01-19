function dual_resolution(g::Graph)
    n = g.n
    s = g.s
    t = g.t

    model = Model(CPLEX.Optimizer)

    @variable(model, x[1:n,1:n] >= 0, binary=true)
    @variable(model, alpha >= 0)
    @variable(model, beta[1:n, 1:n] >= 0)
    @variable(model, gamma >= 0)
    @variable(model, sigma[1:n] >= 0)

    #Constraint 0
    for i = 1:n
        for j = 1:n
            if g.d[i,j] == 0
                @constraint(model, x[i,j] == 0)
                @constraint(model, beta[i,j] == 0)
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

    #Constraint 7
    @constraint(model, sum(x[i,j]*(g.p[i] + g.p[j]) for i = 1:n for j = 1:n if g.d[i,j] != 0) + 2*sum(sigma[i] for i = 1:n) + gamma*g.d2 <= 2*g.S - g.p[s] - g.p[t])

    #Constraint 13
    for i = 1:n
        if i != s && i != t
            @constraint(model, g.ph[i] * (sum(x[i,j] for j = 1:n if g.d[i,j] != 0) + sum(x[j,i] for j = 1:n if g.d[i,j] != 0)) <= gamma + sigma[i])
        end
    end

    #Constraint 14
    @constraint(model, g.ph[s] * (1 + sum(x[s,j] for j = 1:n if g.d[s,j] != 0)) <= gamma + sigma[s])

    #Constraint 15
    @constraint(model, g.ph[t] * (1 + sum(x[i,t] for i = 1:n if g.d[i,t] != 0)) <= gamma + sigma[t])

    #Constraint 16
    for i = 1:n
        for j = 1:n
            if g.d[i,j] != 0
                @constraint(model, g.d[i,j]*x[i,j] <= alpha + beta[i,j])
            end
        end
    end

    #Objective
    @objective(model, Min, sum(g.d[i,j]*x[i,j] for i = 1:n for j = 1:n if g.d[i,j] != 0) + alpha*g.d1 + sum(beta[i,j]*g.D[i,j] for i = 1:n for j = 1:n if g.d[i,j] != 0))

    optimize!(model)

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
        return obj_value, path
        # println(g.name)
        # println("Objective value : ", JuMP.objective_value(model))
        # println("Path from ", g.s, " to ", g.t, " : ", path)
    end
    
    return nothing, nothing, nothing
    # for i = 1:n
    #     for j = 1:n
    #         if value(x[i,j]) == 1
    #             println((i,j))
    #         end
    #     end
    # end
end