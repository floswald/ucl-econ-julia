##  SCRIPT NAME:   mUtilW_cW.jl
##  DATE:          March 2015
#
#   This script calculates the wife's marginal utility with
#   respect to her private consumption.
#
#
##  ALEXANDROS THELOUDIS, UCL

##  ---------------------------------------------------------------------------------------------

function mUtilW_cW(cW, K)

  marginal_utility = alphaCW * (cW^(-gammaCW))


  return marginal_utility

end
