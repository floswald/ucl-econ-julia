##  SCRIPT NAME:   solveValueFunction.jl
##  DATE:          January 2015
#
#   This script solves the household problem by value function iteration.
#   The function obtains the value function for each time period and the policy
#   function (i.e. optimal next-period asset choice). From there I can work out
#   the optimal consumption level. The solution approach taken is by backwards
#   recursion. The optimization each period is carried out using NLopt.
#
#
##  ALEXANDROS THELOUDIS, UCL

##  -------------------------------------------------------------------------------------------
cd("C:\\Dropbox\\JuliaMeeting3")

#   Load packages:
using Grid
using NLopt

#   Load functions:
include("objectiveFunction.jl")
include("utilityH.jl")
include("utilityW.jl")
include("equalityConstraints1.jl")
include("mUtilH_cH.jl")
include("mUtilW_cW.jl")
include("mUtilH_k.jl")
include("mUtilW_k.jl")


##  ---------------------------------------------------------------------------------------------
#   Economy constants:
const r                  = 0.01            # Deterministic interest rate
const beta               = 0.98            # Discount factor
const cminH              = .1              # Min consumption husband (normalised /1000)
const cminW              = .1              # Min consumption wife (normalised /1000)
const kmin               = .1              # Min (inputs to) public consumption (normalised /1000)
const maxAnnWorkingHrsH  = 2.0             # 40 hrs/week times 50 weeks/year (normalised /1000)
const maxAnnWorkingHrsW  = 2.0             # 40 hrs/week times 50 weeks/year (normalised /1000)

#   Pareto weights related constants:
const paretoH            = 0.50            # Pareto weight husband
const paretoW            = 1 - paretoH     # Pareto weight wife

#   Preference-related constants (husband):
const gammaCH            = 1.5             # CRRA on private consumption
const gammaKH            = 1.5             # CRRA on public consumption
const alphaCH            = 0.5             # Weight of private consumption

#   Preference-related constants (wife):
const gammaCW            = 1.5             # CRRA on private consumption
const gammaKW            = 1.5             # CRRA on public consumption
const alphaCW            = 0.5             # Weight of private consumption

#   Solution constant:
const lengthA1           = 50              # Length of asset grid tomorrow


##  ---------------------------------------------------------------------------------------------
#   Current and future state:
A         = 10          # current assets
wH        = 5           # husband's wage today
wW        = 5           # wife's wage today
priceC    = 1.5         # price of private consumption

#   Calculate min and max assets tmr:
# ..lowest assets: cannot die with debts
lbA1      = 0
# ..highest assets: everything saved
ubA1      = (A + wH*maxAnnWorkingHrsH + wW*maxAnnWorkingHrsW - priceC*(cminH + cminW) -kmin)*(1+r)
# ..generate an asset grid for tmr assets (to be used in the EVal interpolation)
xcoord   = range(lbA1, (ubA1 - lbA1)/(lengthA1-1), lengthA1)

#   Expected Value function tomorrow (tmr is terminal period):
EVal1     = zeros(lengthA1)


#  ------------------------------------------------------------------------------------
#  Declare the relationship (linear) between next period's assets and Eval1:
VA!i     = CoordInterpGrid(xcoord, EVal1, BCnil, InterpLinear)


#  ------------------------------------------------------------------------------------
#  Define function to be optimized:

# Set properties of NL optimization problem:
valueFunctionOpt = Opt(:LD_SLSQP, 4)          # algorithm and dimensionality (:LN_COBYLA  :LD_SLSQP)
xtol_rel!(valueFunctionOpt,1e-6)   # tolerance

# Set bounds of optimization:
low_bound = [lbA1 ; cminH ; cminW ; kmin]     # lower bounds
up_bound  = [ubA1 ; Inf   ; Inf   ; Inf ]     # upper bounds
lower_bounds!(valueFunctionOpt,low_bound)     # lower bounds
upper_bounds!(valueFunctionOpt,up_bound)      # upper bounds

# Set objective function, min/max mode, and constraints:
max_objective!(valueFunctionOpt, (x,g) -> objectiveFunction(x, g, A, wH, wW, VA!i))
equality_constraint!(valueFunctionOpt, (x,g) -> equalityConstraints1(x, g, A, wH, wW, priceC), 1e-6)

# Set starting value and initialize optimization:
start_point = [ (ubA1 + lbA1)/2 ;
                (ubA1 - lbA1)/3 ;
                (ubA1 - lbA1)/3 ;
                (ubA1 - lbA1)/3 ]
startpt_valueFunctionOpt     = start_point             # starting point
(fval,x_optimal,return_flag) = optimize(valueFunctionOpt,[startpt_valueFunctionOpt])

# Store results of optimization:
policyA1 = maximum(x_optimal[1])    # policy for next period's assets
policyCH = maximum(x_optimal[2])    # policy for husband's consumption
policyCW = maximum(x_optimal[3])    # policy for wife's consumption
policyK  = maximum(x_optimal[4])    # policy for public consumption
Val      = fval                     # value function today

# Print results:
println("NLopt finished solving. It found an fval=",Val," with CH=",policyCH,", CW=",policyCW,", K=",policyK," and next period's assets are ",policyA1)
