##  SCRIPT NAME:   equalityConstraints1.jl
##  DATE:          March 2015
#
#   This script returns the equality constraints subject to which the
#   family maximises its objective function in each period. If a constraint
#   is of the following form:
#   x + y = a
#   then this must be returned by this script as:
#   result[i] = x + y - a
#   which is set to 0 by the optimiser. The family problem is subject to the
#   following constraints:
#   (a) sequential budget constraint
#
#
##  ALEXANDROS THELOUDIS, UCL

##  ---------------------------------------------------------------------------------------------

function equalityConstraints1(choice_var::Vector, grad::Vector, A0, wH0, wW0, p0)

  # Declare choice variables:
  A!    = choice_var[1]
  consH = choice_var[2]
  consW = choice_var[3]
  kons  = choice_var[4]

  # Calculate gradient:
  if length(grad) > 0
    grad[1] = -1/(1+r)
    grad[2] = -p0
    grad[3] = -p0
    grad[4] = -1
  end

  # Write equality constraints and return results:
  # (a) Sequential budget constraint:
  A0 + wH0*maxAnnWorkingHrsH + wW0*maxAnnWorkingHrsW - p0*(consH + consW) - kons - A!/(1+r)


end
