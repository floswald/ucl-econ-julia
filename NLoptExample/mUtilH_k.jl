##  SCRIPT NAME:   mUtilH_k.jl
##  DATE:          March 2015
#
#   This script calculates the husband's marginal utility with
#   respect to the family's public consumption.
#
#
##  ALEXANDROS THELOUDIS, UCL

##  ---------------------------------------------------------------------------------------------

function mUtilH_k(cH, K)

  marginal_utility = (1-alphaCH) * (K^(-gammaKH))


  return marginal_utility

end
