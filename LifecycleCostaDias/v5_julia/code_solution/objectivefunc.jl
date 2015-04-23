
function objectivefunc(A1, A0, Y,Agrid1,EV1,intfV)

#-------------------------------------------------------------------------------%
# This function returns the following quantity:
# u(c) +  b V( A1)
# where c is calculated from today's assets and tomorrow's assets
#-------------------------------------------------------------------------------%

# ------------------------------------------------------------------------ 
# Declare global we need this file have access to
#global beta r interpMethod            % structural model parameters
#global Agrid1 EV1                     % tomorrow's asset grid and tomorrow's expected value function

# ------------------------------------------------------------------------ 
#Get tomorrow's consumption (cons), the value of left over assets (VA1) and
#total value (u(c) + b * VA1
cons = A0  + Y - (A1)/(1+r);

# In order to obatin the future value function, you need to interpolate on the grid
# of assets. We don't need that for Income as it was "integrated out" at the end
# of the solution loop

    #intfV  = InterpIrregular(Agrid1,EV1, BCnil, InterpLinear); #Already done outside this function
    VA1=interp1extra(Agrid1,EV1,A1,intfV);

value = utility(cons) + beta * VA1;

# ------------------------------------------------------------------------ 
#The optimisation routine that we will use searches for the minimum of the
#function. We want the maximum. So we multiply out function here by -1 so
#that the optimiser will fill the minimum of the negative of our function,
#i.e. the maximum of our functino 

return value = - value;

end

