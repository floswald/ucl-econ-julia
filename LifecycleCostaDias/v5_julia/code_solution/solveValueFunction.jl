function solveValueFunction()

#This function obtains the value function for each time period and
#the policy function (i.e. optimal next-period asset choice) for each time
#period. From there we can work the optimal consumption level.

#The approach taken is by backwards recursion. The optimisation each period
#is carried out using 'fminbnd'. This is an in-built optimiser in Matlab.
#The optimisation routine it uses is known as the 'golden search method'

#global T r tol minCons
#global numPointsA numPointsY Agrid Ygrid incTransitionMrx 
#global Agrid1 EV1

#% ------------------------------------------------------------------------ 
# GENERATE MATRICES TO STORE NUMERICAL APPROXIMATIONS AND INITIATE AS NAN

# Matrices to hold the policy, value and marginal utility functions 
V        = Array(Float64,T+1, numPointsA, numPointsY);
policyA1 = Array(Float64,T,   numPointsA, numPointsY);
policyC  = Array(Float64,T,   numPointsA, numPointsY);        

#Matrices to hold expected value and marginal utility functions 
EV  = Array(Float64,T+1, numPointsA, numPointsY);


#% ------------------------------------------------------------------------ 
#Set the terminal value function and expected value function to 0

EV[T + 1, :,:]  = .0;          # continuation value at T-1
V[T + 1,:,:]    = .0; 
#% ------------------------------------------------------------------------ 
# SOLVE RECURSIVELY THE CONSUMER'S PROBLEM, STARTING AT TIME T-1 AND MOVING
# BACKWARDS TO ZERO, ONE PERIOD AT A TIME


for ixt=T:-1:1                                # Loop from time T-1 to 1
    Agrid1 = vec(Agrid[ixt + 1, :]);          # The grid on assets tomorrow, vector
    
    for ixA = 1:1:numPointsA                  # points on asset grid

        
        # STEP 1. solve problem at grid points in assets and income
        # ---------------------------------------------------------
        A    = Agrid[ixt, ixA];            # assets today
        lbA1 = Agrid[ixt + 1, 1];          # lower bound: assets tomorrow

        for ixY = 1:1:numPointsY               # points on income grid
                        
            # Value of income and information for optimisation            
            Y    = Ygrid[ixt, ixY];            # income today            
            ubA1 = (A + Y - minCons)*(1+r);    # upper bound: assets tomorrow
            EV1  = vec(EV[ixt + 1,:, ixY]);    # relevant section of EV matrix (in assets tomorrow), vector
            
            intfV  = InterpIrregular(Agrid1,EV1, BCnil, InterpLinear); # Define the interpolation

            # Compute solution 
            if (ubA1 - lbA1 < minCons)                       # if liquidity constrained
                negV = objectivefunc(lbA1, A, Y,Agrid1,EV1,intfV); 
                policyA1[ixt,ixA,ixY] = lbA1;
            else                                             # if interior solution

                function obj(A1)
                    return objectivefunc(A1, A, Y, Agrid1,EV1,intfV)
                end

                Res = optimize(obj,lbA1,ubA1);          #Find the solution...

                policyA1[ixt,ixA,ixY]=Res.minimum;
                negV=Res.f_minimum;

                if isnan(negV)
                    println("Sendo error en ixt=$ixt, ixA=$ixA, ixY=$ixY")
                end

            end # if (ubA1 - lbA1 < minCons)         

            # Store solution and its value
            policyC[ixt, ixA, ixY] = A + Y - policyA1[ixt, ixA, ixY]/(1+r);
            V[ixt, ixA, ixY]       = -negV;          
        end #ixY


        # STEP 2. integrate out income today conditional on income
        # yesterday to get EV and EdU
        # --------------------------------------------------------
        realisedV = squeeze(V[ixt, ixA, :],1);
        for ixY = 1:1:numPointsY   
            d=incTransitionMrx[ixY,:]*realisedV'
            EV[ixt, ixA, ixY]  = d[1,1];
        end #ixY
    end #ixA

    println("Passed period $ixt of $T")
    
end #ixt

return policyA1, policyC, V, EV 

end #function

