##  SCRIPT NAME:   mUtilH_cH.jl
##  DATE:          March 2015
#
#   This script calculates the husband's marginal utility with
#   respect to his private consumption.
#
#
##  ALEXANDROS THELOUDIS, UCL

##  ---------------------------------------------------------------------------------------------

function mUtilH_cH(cH, K)

  marginal_utility = alphaCH * (cH^(-gammaCH))


  return marginal_utility

end
