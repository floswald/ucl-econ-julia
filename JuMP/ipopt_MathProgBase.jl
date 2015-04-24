## This file implements Example H071 using the MathProgBase Interface.
## This is the example code from the MathProgBase test file from https://github.com/JuliaOpt/MathProgBase.jl/blob/master/test/nlp.jl , slightly modified to remove the tests.
## This interface is solver independent, but will use Ipopt as the solver here.

using MathProgBase
import MathProgBase.MathProgSolverInterface   ##In Julia, importing means that the prefix 'MathProgSolverInterface.' must be included before calling any functions within this module
using Ipopt

tic()  #Start the clock

# Here the type represents the complete instance, but it
# could also store instance data.
# The abstract type AbstractNLPEvaluator is used by solvers for accessing the objective function f and constraints g.
# Solvers may query the value, gradients, Hessian-vector products, and the Hessian of the Lagrangian.
type HS071 <: MathProgSolverInterface.AbstractNLPEvaluator
end

# hs071
# min x1 * x4 * (x1 + x2 + x3) + x3
# st  x1 * x2 * x3 * x4 >= 25
#     x1^2 + x2^2 + x3^2 + x4^2 = 40
#     1 <= x1, x2, x3, x4 <= 5
# Start at (1,5,5,1)
# End at (1.000..., 4.743..., 3.821..., 1.379...)

## Initialize must be called before any other methods. The vector requested_features lists features requested by the solver.
function MathProgSolverInterface.initialize(d::HS071, requested_features::Vector{Symbol})
    for feat in requested_features
        if !(feat in [:Grad, :Jac, :Hess])
            error("Unsupported feature $feat")
        end
    end
end

## Returns which features are available for the model HS071
MathProgSolverInterface.features_available(d::HS071) = [:Grad, :Jac, :Hess]

################################
# Define the Objective Function
################################
MathProgSolverInterface.eval_f(d::HS071, x) = x[1] * x[4] * (x[1] + x[2] + x[3]) + x[3]

###################
# Constraint set
##################
function MathProgSolverInterface.eval_g(d::HS071, g, x)
    g[1] = x[1]   * x[2]   * x[3]   * x[4]
    g[2] = x[1]^2 + x[2]^2 + x[3]^2 + x[4]^2
end

#####################################
# Gradient of the Objective Function
#####################################
function MathProgSolverInterface.eval_grad_f(d::HS071, grad_f, x)
    grad_f[1] = x[1] * x[4] + x[4] * (x[1] + x[2] + x[3])
    grad_f[2] = x[1] * x[4]
    grad_f[3] = x[1] * x[4] + 1
    grad_f[4] = x[1] * (x[1] + x[2] + x[3])
end

###################################################################################################################################
##These functions collect the sparsity structures of the Jacobian matrix and the Hessian of the Lagrangian (row first then columns)
MathProgSolverInterface.jac_structure(d::HS071) = [1,1,1,1,2,2,2,2],[1,2,3,4,1,2,3,4]
# lower triangle only
MathProgSolverInterface.hesslag_structure(d::HS071) = [1,2,2,3,3,3,4,4,4,4],[1,1,2,1,2,3,1,2,3,4]

#########################################################################################################################
##OR you could define the values of the non-zero elements. The function below does this for the Jacobian matrix and saves
## them to the matrix J. There is no need to create the matrix J.
####################################################################
function MathProgSolverInterface.eval_jac_g(d::HS071, J, x)
    # Constraint (row) 1
    J[1] = x[2]*x[3]*x[4]  # 1,1
    J[2] = x[1]*x[3]*x[4]  # 1,2
    J[3] = x[1]*x[2]*x[4]  # 1,3
    J[4] = x[1]*x[2]*x[3]  # 1,4
    # Constraint (row) 2
    J[5] = 2*x[1]  # 2,1
    J[6] = 2*x[2]  # 2,2
    J[7] = 2*x[3]  # 2,3
    J[8] = 2*x[4]  # 2,4
end

####################################################################################################
# Similarly for the Hessian -  σ = weight on objective function, μ = multipliers on the constraints
####################################################################################################
function MathProgSolverInterface.eval_hesslag(d::HS071, H, x, σ, μ)
    # Again, only lower left triangle (since it is symmetric)
    # Objective
    H[1] = σ * (2*x[4])               # 1,1
    H[2] = σ * (  x[4])               # 2,1
    H[3] = 0                          # 2,2
    H[4] = σ * (  x[4])               # 3,1
    H[5] = 0                          # 3,2
    H[6] = 0                          # 3,3
    H[7] = σ* (2*x[1] + x[2] + x[3])  # 4,1
    H[8] = σ * (  x[1])               # 4,2
    H[9] = σ * (  x[1])               # 4,3
    H[10] = 0                         # 4,4

    # First constraint
    H[2] += μ[1] * (x[3] * x[4])  # 2,1
    H[4] += μ[1] * (x[2] * x[4])  # 3,1
    H[5] += μ[1] * (x[1] * x[4])  # 3,2
    H[7] += μ[1] * (x[2] * x[3])  # 4,1
    H[8] += μ[1] * (x[1] * x[3])  # 4,2
    H[9] += μ[1] * (x[1] * x[2])  # 4,3

    # Second constraint
    H[1]  += μ[2] * 2  # 1,1
    H[3]  += μ[2] * 2  # 2,2
    H[6]  += μ[2] * 2  # 3,3
    H[10] += μ[2] * 2  # 4,4

end


# Now put this all together to define and solve the model.
# First, define a model object and the solver to use
    m = MathProgSolverInterface.model(MathProgBase.defaultNLPsolver )
# The vector of lower bounds on the variables
    l = [1,1,1,1]
# The vector of upper bounds on the variables
    u = [5,5,5,5]
# Vector of lower bounds on the constraints
    lb = [25, 40]
# Vector of constraint upper bounds
    ub = [Inf, 40]
# Load the nonlinear programming problem into the model - syntax for loadnonlinearproblem!() is
# loadnonlinearproblem!(m::AbstractMathProgModel, numVar, numConstr, l, u, lb, ub, sense, d::AbstractNLPEvaluator)
    MathProgSolverInterface.loadnonlinearproblem!(m, 4, 2, l, u, lb, ub, :Min, HS071())
## Set the initial values for the optimisation
    MathProgSolverInterface.setwarmstart!(m,[1,5,5,1])

## Solve the optimisation problem
    MathProgSolverInterface.optimize!(m)
## Let's see the termination status after solving
    stat = MathProgSolverInterface.status(m)

# Now let's see the solution.
    x = MathProgSolverInterface.getsolution(m)

toc()

println(x)

