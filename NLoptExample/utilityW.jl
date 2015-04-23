##  SCRIPT NAME:   utilityW.jl
##  DATE:          March 2015
#
#   This script calculates the felicity the wife enjoys for every value
#   of her private consumption and the family level public consumption.
#
#
##  ALEXANDROS THELOUDIS, UCL

##  ---------------------------------------------------------------------------------------------

function utilityW(cW, K)


  utility = alphaCW*(cW^(1-gammaCW))/(1-gammaCW) + (1-alphaCW)*(K^(1-gammaKW))/(1-gammaKW)


  return utility

end
