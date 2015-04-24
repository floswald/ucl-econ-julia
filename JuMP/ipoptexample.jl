#Ipopt.jl example - This is taken from the documentation available on:
# http://ipoptjl.readthedocs.org/en/latest/ipopt.html#

# Install Package if you don't have it
Pkg.add("Ipopt")

tic()
# Load the package
using Ipopt

# HS071
# min x1 * x4 * (x1 + x2 + x3) + x3
# st  x1 * x2 * x3 * x4 >= 25
#     x1^2 + x2^2 + x3^2 + x4^2 = 40
#     1 <= x1, x2, x3, x4 <= 5
# Start at (1,5,5,1)
# End at (1.000..., 4.743..., 3.821..., 1.379...)

## Set the number of arguments to solve for and bounds for these
## Aside for new Julia users: the '.0' after each number is important - it tells Julia
## that the vectors x_L, etc are of Float type.
## Note that Julia is a functional language, with methods (e.g. a function) working for the types
## for which it is defined (e.g. integer, Float). So, a function defined for an integer will not
## work if you enter a Float variable, and vice-versa.
n = 4
x_L = [1.0, 1.0, 1.0, 1.0]
x_U = [5.0, 5.0, 5.0, 5.0]


# Set the number of constraints and bounds on these - Note that equality
# constraints can be accommodated for by setting the upper bound = lower bound
m = 2
g_L = [25.0, 40.0]
g_U = [2.0e19, 40.0]

# Set the objective function - Returns the value of the objective function at the current solution x
# Note that the '::Vector{Float64}' below defines the type of x. This is a Julia-specific thing, and
# need not be defined. The function would be created even without including the type for x.
function eval_f(x::Vector{Float64})
  return x[1] * x[4] * (x[1] + x[2] + x[3]) + x[3]
end

# Set the constraint set - Sets the value of the constraint functions g at the current solution x
# Note that there is no need to create a separate vector g in which to save the constraints.
# The values of g are set 'in-place'
function eval_g(x::Vector{Float64}, g::Vector{Float64})
  g[1] = x[1]   * x[2]   * x[3]   * x[4]
  g[2] = x[1]^2 + x[2]^2 + x[3]^2 + x[4]^2
end

# Gradient of the objective function - Sets the value of the gradient of the objective function at the current solution x
# Note that there is no need to create a vector grad_f to fill in. The values of grad_f are set 'in-place'.
function eval_grad_f(x::Vector{Float64}, grad_f::Vector{Float64})
  grad_f[1] = x[1] * x[4] + x[4] * (x[1] + x[2] + x[3])
  grad_f[2] = x[1] * x[4]
  grad_f[3] = x[1] * x[4] + 1
  grad_f[4] = x[1] * (x[1] + x[2] + x[3])
end


# Gradient of the constraint set (Jacobian matrix). There are two ways of setting this:
# (1) Under the ':Structure' mode, indicate the rows and columns with non-zero values (i.e. tell Ipopt
# about the sparsity structure of the matrix)
# (2) Undet the ':Values' mode, tell it the actual values.
# Each non-zero value is indexed, starting from 1.
###########################################
# General specification is as follows:
###########################################
# function eval_jac_g(
#  x::Vector{Float64},         # Current solution
#  mode,                       # Either :Structure or :Values
#  rows::Vector{Int32},        # Sparsity structure - row indices
#  cols::Vector{Int32},        # Sparsity structure - column indices
#  values::Vector{Float64})    # The values of the Hessian

#  if mode == :Structure
    # rows[...] = ...
    # ...
    # cols[...] = ...
#  else
    # values[...] = ...
#  end
# end
#############################################
function eval_jac_g(x, mode, rows, cols, values)
  if mode == :Structure
    # Constraint (row) 1
    rows[1] = 1; cols[1] = 1
    rows[2] = 1; cols[2] = 2
    rows[3] = 1; cols[3] = 3
    rows[4] = 1; cols[4] = 4
    # Constraint (row) 2
    rows[5] = 2; cols[5] = 1
    rows[6] = 2; cols[6] = 2
    rows[7] = 2; cols[7] = 3
    rows[8] = 2; cols[8] = 4
  else
    # Constraint (row) 1
    values[1] = x[2]*x[3]*x[4]  # 1,1
    values[2] = x[1]*x[3]*x[4]  # 1,2
    values[3] = x[1]*x[2]*x[4]  # 1,3
    values[4] = x[1]*x[2]*x[3]  # 1,4
    # Constraint (row) 2
    values[5] = 2*x[1]  # 2,1
    values[6] = 2*x[2]  # 2,2
    values[7] = 2*x[3]  # 2,3
    values[8] = 2*x[4]  # 2,4
  end
end

# Hessian of the Lagrangian - this is an optional object. In terms of the Julia syntax, it is
# similar to the Jacobian.
# The general definition and types for the objects in the function below are:
# function eval_h(
# x::Vector{Float64},         # Current solution
# mode,                       # Either :Structure or :Values
# rows::Vector{Int32},        # Sparsity structure - row indices
# cols::Vector{Int32},        # Sparsity structure - column indices
# obj_factor::Float64,        # Lagrangian multiplier for objective
# lambda::Vector{Float64},    # Multipliers for each constraint
# values::Vector{Float64})    # The values of the Hessian
# if....else...end            # the function text goes here
# end                         # close the function

function eval_h(x, mode, rows, cols, obj_factor, lambda, values)
  if mode == :Structure
    # Symmetric matrix, fill the lower left triangle only
    idx = 1
    for row = 1:4
      for col = 1:row
        rows[idx] = row
        cols[idx] = col
        idx += 1
      end
    end
  else
    # Again, only lower left triangle
    # Objective
    values[1] = obj_factor * (2*x[4])  # 1,1
    values[2] = obj_factor * (  x[4])  # 2,1
    values[3] = 0                      # 2,2
    values[4] = obj_factor * (  x[4])  # 3,1
    values[5] = 0                      # 3,2
    values[6] = 0                      # 3,3
    values[7] = obj_factor * (2*x[1] + x[2] + x[3])  # 4,1
    values[8] = obj_factor * (  x[1])  # 4,2
    values[9] = obj_factor * (  x[1])  # 4,3
    values[10] = 0                     # 4,4

    # First constraint
    values[2] += lambda[1] * (x[3] * x[4])  # 2,1
    values[4] += lambda[1] * (x[2] * x[4])  # 3,1
    values[5] += lambda[1] * (x[1] * x[4])  # 3,2
    values[7] += lambda[1] * (x[2] * x[3])  # 4,1
    values[8] += lambda[1] * (x[1] * x[3])  # 4,2
    values[9] += lambda[1] * (x[1] * x[2])  # 4,3

    # Second constraint
    values[1]  += lambda[2] * 2  # 1,1
    values[3]  += lambda[2] * 2  # 2,2
    values[6]  += lambda[2] * 2  # 3,3
    values[10] += lambda[2] * 2  # 4,4
  end
end

# Creates and returns an IpoptProblem with the given options. Raises error if something goes wrong during construction.
# If you do not provide a callback for the Hessian, you must set the Hessian approximation option:
# addOption(prob, "hessian_approximation", "limited-memory")
prob = createProblem(n, x_L, x_U, m, g_L, g_U, 8, 10,
                     eval_f, eval_g, eval_grad_f, eval_jac_g, eval_h)

# Set starting solution
prob.x = [1.0, 5.0, 5.0, 1.0]

# Solve
status = solveProblem(prob)

# Display the solution
println(Ipopt.ApplicationReturnStatus[status])
println(prob.x)
println(prob.obj_val)

toc()
