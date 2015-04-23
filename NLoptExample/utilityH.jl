##  SCRIPT NAME:   utilityH.jl
##  DATE:          March 2015
#
#   This script calculates the felicity the husband enjoys for every value
#   of his private consumption and the family level public consumption.
#
#
##  ALEXANDROS THELOUDIS, UCL

##  ---------------------------------------------------------------------------------------------

function utilityH(cH, K)


  utility = alphaCH*(cH^(1-gammaCH))/(1-gammaCH) + (1-alphaCH)*(K^(1-gammaKH))/(1-gammaKH)


  return utility

end
