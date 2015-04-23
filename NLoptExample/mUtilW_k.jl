##  SCRIPT NAME:   mUtilW_k.jl
##  DATE:          March 2015
#
#   This script calculates the wife's marginal utility with
#   respect to the family's public consumption.
#
#
##  ALEXANDROS THELOUDIS, UCL

##  ---------------------------------------------------------------------------------------------

function mUtilW_k(cW, K)

  marginal_utility = (1-alphaCW) * (K^(-gammaKW))


  return marginal_utility

end
