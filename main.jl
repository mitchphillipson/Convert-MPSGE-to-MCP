using MPSGE
using JuMP
using PATHSolver

PZ_NUMS = [1 3; 2 2]

function mpsge_model()

    M = MPSGEModel()

    @parameter(M, tax,0)

    @sector(M, X[i=1:2])

    @commodities(M, begin
        PX[i=1:2]
        PY[i=1:2]
        PZ[i=1:2]
    end)

    @consumer(M, RA)

    @production(M, X[i=1:2], [s=2, t=0, nest=>s=1], begin
        @output(PY[i], 14, t)
        @input(PX[lamia=1:2], 5, s, taxes = [Tax(RA, tax)])
        @input(PZ[j=1:2], PZ_NUMS[i,j], nest)
    end)

    @demand(M, RA, begin
        @final_demand(PY[j=1:2], 14)
        @endowment(PX[j=1:2], 10)
        @endowment(PZ[j=1:2], sum(PZ_NUMS[:,j]))
    end)

    return M
end

M_mps = mpsge_model()

solve!(M_mps, cumulative_iteration_limit=0)

function mcp_model()
    M = Model(PATHSolver.Optimizer)

    @variables(M, begin
        X[i=1:2], (start = 1, lower_bound = 0)
        PX[i=1:2], (start = 1, lower_bound = 0)
        PY[i=1:2], (start = 1, lower_bound = 0)
        PZ[i=1:2], (start = 1, lower_bound = 0)
        RA, (start = 28)

        tax, (start = 0)
    end)

    fix(tax, 0, force=true)

    # X - Zero Profit

    # @input(PZ[j=1:2], 2, nest)
    @expression(M, ucf_X_nest[i=1:2], prod(PZ[j]^(PZ_NUMS[i,j]/sum(PZ_NUMS[i,:])) for j=1:2))

    @expression(M, px_price[i=1:2], PX[i]*(1+tax))

    @expression(M, ucf_X_s[i=1:2], (sum(PZ_NUMS[i,:])/(10+sum(PZ_NUMS[i,:]))*ucf_X_nest[i]^(1/(1-2)) + sum(5/(10+sum(PZ_NUMS[i,:]))*(px_price[j])^(1/(1-2)) for j=1:2))^(1-2))
    
    @constraint(M, x_zp[i=1:2], 14*ucf_X_s[i] -  14*PY[i] ⟂ X[i])

    
    # PY - Market Clearance
    @constraint(M, py_mc[i=1:2], 14*X[i] - 14/28*RA/PY[i] ⟂ PY[i])

    # PX - Market Clearance
    @constraint(M, px_mc[i=1:2], sum(5/14 * 14 * -(ucf_X_s[j]/px_price[i])^2*X[j] for j=1:2) + 10 ⟂ PX[i])

    

    # PZ - Market Clearance
    @constraint(M, pz_mc[i=1:2], sum(PZ_NUMS[j,i]/14 * 14 * -(ucf_X_s[j]/ucf_X_nest[j])^2*(ucf_X_nest[j]/PZ[i])^1*X[j] for j=1:2) + sum(PZ_NUMS[:,i]) ⟂ PZ[i])

    # Income Balance
    @constraint(M, income_balance, RA - (sum(PX[j]*10 for j=1:2) + sum(PZ[j]*sum(PZ_NUMS[:,j]) for j=1:2)) + sum(5* tax*PX[j]*X[i]*(ucf_X_s[i]/px_price[j])^2 for i=1:2,j=1:2)  ⟂ RA)

    return M
end

M_mcp = mcp_model()

JuMP.set_attribute(M_mcp, "cumulative_iteration_limit", 0)
optimize!(M_mcp)


JuMP.set_attribute(M_mcp, "cumulative_iteration_limit", 10_000)


set_value!(M_mps[:tax], .1)
fix(M_mcp[:tax], .1, force=true)

solve!(M_mps)

optimize!(M_mcp)

value.(M_mcp[:X])
value.(M_mps[:X])

value.(M_mcp[:PX])
value.(M_mps[:PX])