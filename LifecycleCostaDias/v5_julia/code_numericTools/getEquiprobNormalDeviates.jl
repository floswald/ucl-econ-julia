##  SCRIPT NAME:   getEquiprobNormalDeviates.jl
##  DATE:          Jnauary 2015
#
#   This script generates N equiprobable line segments in the stationary
#   univariate normal distribution of income Y:
#          Y ~ N [ alpha/1-rho  ;  sigma_epsilon^2/(1-rho)^2 ]
#   where the mean of the above distribution is declared below as 'Ymean'
#   and the standard deviation as 'Ysd'. It delivers two arrays: 'Bounaries'
#   which holds the N+1 bounds of the N equiprobable segments, and
#   'ExpValues' which holds the expected value in each segment.
#
#
##  ALEXANDROS THELOUDIS, UCL
#
#  Needs: A constant called "normBnd"

##  ---------------------------------------------------------------------------------------------

function getEquiprobNormalDeviates(Ymean, Ysd, N)

  #  Initialise output:
  Boundaries = zeros(N+1)
  ExpValues  = zeros(N)


  #  --------------------------------------------------------------------------------------------
  #  Calculate applicable support of univariate normal distribution (first and
  #  last point in the range that I divide equiprobably). Note: in principle
  #  this is +-infinity but I am ruling some extremely unlikely values out:
  Boundaries[1]   = Ymean - normBnd * Ysd
  Boundaries[N+1] = Ymean + normBnd * Ysd


  #  --------------------------------------------------------------------------------------------
  #  Sequentially calculate the remaining bounds of the equiprobable regions:
  for ixi = 2:1:N
    Boundaries[ixi] = quantile(TruncatedNormal(Ymean,Ysd,Boundaries[1],Boundaries[N+1]), (ixi-1)/N)
  end


  #  --------------------------------------------------------------------------------------------
  #  Calculate the expected value of a random variable X ~ N[Ymean,Ysd^2] which
  #  is bounded by the equiprobable line segment above. Note: this is essentially
  #  the mean of a univariate 2-sided truncated normal distribution; read more
  #  on http://en.wikipedia.org/wiki/Truncated_normal_distribution
  for ixi = 1:1:N
    ExpValues[ixi] = Ymean + Ysd * N * (pdf(Normal(),(Boundaries[ixi]-Ymean)/Ysd) - pdf(Normal(),(Boundaries[ixi+1]-Ymean)/Ysd))
  end

  return (Boundaries, ExpValues)

end

