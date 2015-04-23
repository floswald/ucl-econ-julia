function interp2D( x1Vec0, x2Vec0, yMat, x1Val, x2Val )
#UNTITLED Summary of this function goes here
#   Detailed explanation goes here

x1Vec=vec(x1Vec0)
x2Vec=vec(x2Vec0)

#Check that both x1Vec and x2Vec are n_i * 1 vectors
sizex1Vec = size(x1Vec);
sizex2Vec = size(x2Vec);
sizeyMat = size(yMat);


#println("Enter $sizex1Vec and $sizex2Vec and $sizeyMat")

# if (sizex1Vec[1]<=1) || (sizex2Vec[1]<=1) || (sizex1Vec[2]!=1) || (sizex2Vec[2]!=1)
#     error("Both x1Vec and x2Vec need to be n_i by 1 vectors. The length of the vectors may differ")
# end

#Check whether x1Vec and x2Vec are monotonic
# isx1VecMonotone = all(diff[x1Vec]>=0);
# isx2VecMonotone = all(diff[x2Vec]>=0);

# if (isx1VecMonotone*isx2VecMonotone == 0)
#     error("One of x1Vec or x2Vec is non-monotonic");
# end

n1 = sizex1Vec[1];

function nonNeg(args)
    return args>=0
end

#Find where x1val fits into x1Vec and pick out the vectors in yMat
#conditional on these values
if (x1Val<=x1Vec[1])
    x1LowerIndex = 1; 
    x1UpperIndex = 2;        
elseif (x1Val>=x1Vec[n1])
    x1UpperIndex = n1;
    x1LowerIndex = n1 - 1;
else #x1Val is between the lower and the upper components of x1Vec
    #Get a vect
    x1VecLessX1Val = x1Vec - x1Val;

    # Get the smallest non-negative 
    posI=find(nonNeg, x1VecLessX1Val)
    (a,b)=findmin(x1VecLessX1Val[posI],1)
    x1UpperIndex=posI[b];
    x1LowerIndex = x1UpperIndex - 1.;
end

newYMatAfterInterpOnX2 = Array(Float64,2, 1);

yMatCondOnx1Lower = yMat[x1LowerIndex, :];
yMatCondOnx1Upper = yMat[x1UpperIndex, :];

intfL  = InterpIrregular(vec(x2Vec),vec(yMatCondOnx1Lower), BCnil, InterpLinear);
intfH  = InterpIrregular(vec(x2Vec),vec(yMatCondOnx1Upper), BCnil, InterpLinear);

# The last bit of interpolation!

newYMatAfterInterpOnX2[1]  =interp1extra(vec(x2Vec),vec(yMatCondOnx1Lower),x2Val,intfL);
newYMatAfterInterpOnX2[2]  =interp1extra(vec(x2Vec),vec(yMatCondOnx1Upper),x2Val,intfH);

newX1VecShort=[x1Vec[x1LowerIndex] x1Vec[x1UpperIndex] ]'

intfF  = InterpIrregular(vec(newX1VecShort),vec(newYMatAfterInterpOnX2), BCnil, InterpLinear); 
yVal   = interp1extra(vec(newX1VecShort),vec(newYMatAfterInterpOnX2),x1Val,intfF); #  intfF[x1Val];

return yVal
    
end

