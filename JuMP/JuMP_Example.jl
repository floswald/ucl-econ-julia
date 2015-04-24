# JuMP example for the Julia Working Group

# Reproduce the example HS071.
# HS071
# min x1 * x4 * (x1 + x2 + x3) + x3
# st  x1 * x2 * x3 * x4 >= 25
#     x1^2 + x2^2 + x3^2 + x4^2 = 40
#     1 <= x1, x2, x3, x4 <= 5
# Start at (1,5,5,1)
# End at (1.000..., 4.743..., 3.821..., 1.379...)

# First install the package, and a solver
#Pkg.add("JuMP")
#Pkg.add("Ipopt") #solver

# Load these packages
using JuMP
using Ipopt
#using NLopt

tic()  #start clock

#Create the model object
m = Model()
#m = Model(solver = NLoptSolver(algorithm=:LD_SLSQP))
#m=Model(solver=IpoptSolver(tol=1e-8, max_iter=5000)) # can use the argument solver to pass solver options

#Defining Variables - These are defined using the @defVar macro. The first argument of
# @defVar(m,x) will always be the model to associate the variable with. You can also
# impose the bounds on the variable here; and start values (using the start option)
@defVar(m, 1<=x[1:4]<=5)
# Set start value
setValue(x[1],1)
setValue(x[2],5)
setValue(x[3],5)
setValue(x[4],1)

# One can also create variables manually one by one as:
# x = Variable(m::Model, lower::Number, upper::Number, category::Symbol, name::String)
# Bounds can also be set using setLower(x::Variable, lower) and setUpper(x::Variable, lower)


# Objective function - I will show the syntax for a nonlinear objective function. For the
# syntax for linear objective functions, see the documentation. Essentially one needs to
# use setObjective() or @setObjective() rather than @setNLObjective()
# For nonlinear objectives, these are defined using the macro @setNLObjective()

@setNLObjective(m::Model, :Min, x[1]*x[4]*(x[1]+x[2]+x[3])+x[3] )

# Next, set constraints.
# Linear constraints are added through the macro @addConstraint()
# Nonlinear constraints are added through the @addNLConstraint() macro
@addNLConstraint(m::Model, x[1]*x[2]*x[3]*x[4]>=25 )
@addNLConstraint(m::Model, x[1]^2 + x[2]^2 + x[3]^2 + x[4]^2 == 40)

# Finally, solve the model
solve(m)

# Let's look at the solutions
println("x=", getValue(x))

toc()
