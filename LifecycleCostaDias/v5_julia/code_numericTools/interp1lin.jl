function interp1extra(Xgrid,Ygrid,val,interpObj)
# I just don't know why the Grid guys didn't add linear extrapolation!
# interpObj: You should call the interpolation function first, this
#            program just perform the linear extrapolation if you are
#            out of the top boundary of Xgrid
# Example:
#   intfV  = InterpIrregular(Agrid1,EV1, BCnil, InterpLinear);
#   VA1    = interp1extra(Agrid1,EV1,A1,intfV); # instead of intfV[A1]


    if val<Xgrid[1] # Below the minimum point
    	Vval=Ygrid[1];
    elseif val>Xgrid[1] && val<Xgrid[length(Xgrid)]  # Interior Solution
		#interpObj  = InterpIrregular(Xgrid,Ygrid, BCnil, InterpLinear); #Already done outside this function
        Vval=interpObj[val]
    else  # Above the maximum point (linear extrapolation: li=x*b+c)
        x2=Xgrid[length(Xgrid)]
        x1 =Xgrid[length(Xgrid)-1]

        l2 =Ygrid[length(Xgrid)]
        l1 =Ygrid[length(Xgrid)-1]

        b=(l2-l1)/(x2-x1);
        c=l2-b*x2;

    	Vval=val*b+c;
    end

    return Vval

end