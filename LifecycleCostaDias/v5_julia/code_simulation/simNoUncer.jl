function simNoUncer(policyA1,EV,startingA)

# This function takes the policy functions and value functions in an environment
# where there is no uncertainty, along with starting assets and returns 
# simulated paths of income, consumption assets and value

#% ------------------------------------------------------------------------ 
# Declare global we need this file have access to
#global T r 
#global Agrid  Ygrid numSims interpMethod 


#% ------------------------------------------------------------------------
# Initialise arrays that will hold the paths of income consumption, value
# and assets

# Arguments for output
y = Array(Float64,T, numSims);            # income
c = Array(Float64,T, numSims);            # consumption
v = Array(Float64,T, numSims);            # value
a = Array(Float64,T + 1,numSims);         # this is the path at the start of each period, so we include the 'start' of death

#% ------------------------------------------------------------------------
# Obtain paths using the initial condition and the policy and value
# functions
#-------------------------------------------------------------------------%

for s = 1:1:numSims                   # loop through individuals
    a[1, s] = startingA;   
    for t = 1:1:T                     # loop through time periods for a particular individual
		intfV  = InterpIrregular(vec(Agrid[t, :]),vec(EV[t, :]), BCnil, InterpLinear);       # Define the interpolation func for Value
		intfP  = InterpIrregular(vec(Agrid[t, :]),vec(policyA1[t, :]), BCnil, InterpLinear); # Define the interpolation func for Polic

        y[t  , s]   = Ygrid[t];     
        v[t  , s]   =interp1extra(vec(Agrid[t, :]),vec(EV[t, :]),a[t, s],intfV);
        a[t+1, s]   =interp1extra(vec(Agrid[t, :]),vec(policyA1[t, :]),a[t, s],intfP); 
        c[t  , s] = a[t, s]  + y[t, s] - (a[t+1, s]/(1+r));
    end   #t      
end # s
  
return  y, c, a, v
 
end
