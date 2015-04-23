
function checkInputs()
#Check various inputs

 println("Let's check inputs...!")


# Check that the hard coded income and PDF have the correct dimension
if (isUncertainty == 1) && (uncertaintyMethod == 0)
    if (size(hcIncome, 1)!= numPointsY) || (size(hcIncPDF, 1)!= numPointsY)
        error("hcIncome and must have dimension numPointsY * 1 and numPointsY * 1 respectively")
    end
end

#
# Check if uncertainty is off that numPointsY is set equal to 1
if (numPointsY!=1) && (isUncertainty == 0)
    error("There is no uncertainty but numPointsY was not set to 1.")
end

#
# Check that the value of rho is not greater than 1. Warn if close to 1
if (rho>0.999) || (rho <-0.999)
    if (rho>1)  || (rho < 1)     
    error("rho is greater than 1. This code solves a stationary income process. rho greater than 1 implies a non-stationary process")
    else
    warning("rho is greater than 0.99. This code solves a stationary income process. As rho gets closer to 1 the process becomes nonstationary - possibility of numerical instability")        
    end
end

#
# Check that the standard deviation is not too small
if (sigma<1e-10) && (isUncertainty == 1)
    if (sigma<=0) 
        error("sigma is less than or equal to zero")
    else
        warning("sigma is very small and close to zero - possibility of numerical instability. Consider turning uncertainty off")
    end
end



# Check a number of inputs that need to be either 0 or 1
if (linearise != 0) && (linearise !=1)
    error("linearise should be either 0 or 1")
end

if (uncertaintyMethod != 0) && (uncertaintyMethod !=1)
    error("uncertaintyMethod should be either 0 or 1")
end

if (borrowingAllowed != 0) && (borrowingAllowed !=1)
    error("borrowingAllowed should be either 0 or 1")
end

if (isUncertainty != 0) && (isUncertainty !=1)
    error("isUncertainty should be either 0 or 1")
end

println("Inputs are fine :)")

end
