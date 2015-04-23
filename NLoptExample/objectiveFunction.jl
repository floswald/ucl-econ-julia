##  SCRIPT NAME:   objectiveFunction.jl
##  DATE:          January 2015
#
#   This script returns the objective function in each time period, that is:
#   μH uH(cH,K) + μW uW(cW,K) + β EV(a')
#   where the way it is written I assume the choice variables are a', cH, cW,
#   and K (i.e. next period's assets, private consumption, and public consumption).
#   It is the felicity the family enjoys if current state variables are: a
#   (assets), wH, and wW (wages), and the family chooses next period assets a',
#   private consumption cH and cW, and public consumption. The continuation value
#   is inputted in the function.
#
#
##  ALEXANDROS THELOUDIS, UCL

##  ---------------------------------------------------------------------------------------------

function objectiveFunction(choice_var::Vector, grad::Vector, A0, wH0, wW0, interpol_A!_EVal!)

  # Declare choice variables:
  A!        = choice_var[1]
  consH     = choice_var[2]
  consW     = choice_var[3]
  kons      = choice_var[4]

  # Calculate gradient:
  if length(grad) > 0
    (value,derivative) = valgrad(interpol_A!_EVal!, A!)

    grad[1] = beta * derivative
    grad[2] = paretoH * mUtilH_cH(consH, kons)
    grad[3] = paretoW * mUtilW_cW(consW, kons)
    grad[4] = paretoH * mUtilH_k(consH, kons) + paretoW * mUtilW_k(consW, kons)
  end

  # Calculate objective (with interpolated future value function):
  objective = paretoH * utilityH(consH, kons) + paretoW * utilityW(consW, kons) + beta * interpol_A!_EVal![A!]


  return objective

end


