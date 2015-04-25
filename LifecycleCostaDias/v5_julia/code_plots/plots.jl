
function plots()

actualDir=pwd() 
cd(joinpath(actualDir,"output","images"))



xi=Agrid[whichyear, plotNode1:plotNodeLast]

# Policy Fucntion (Consumption) -----------------------------------
y1=policyC[whichyear, plotNode1:plotNodeLast, 1]
y2=policyC[whichyear, plotNode1:plotNodeLast, numPointsY]

bsample = DataFrame(A1=vec(xi),LowInc=vec(y1),HighInc=vec(y2))
bsample[:mtindex] = 1:size(bsample,1)     # Add an identifier
datamia=stack(bsample,[2:3])              # Reshape the data so we can use colors!

plot1=plot(datamia,x="A1",y="value",color="variable",Geom.line,Guide.title("Policy Function A1"))
draw(PNG("policyA1.png", 24cm, 12cm), plot1)




# Value Function (A1) --------------------------------------------

y1=val[whichyear, plotNode1:plotNodeLast, 1]
y2=val[whichyear, plotNode1:plotNodeLast, numPointsY]

bsample = DataFrame(Val=vec(xi),LowInc=vec(y1),HighInc=vec(y2))
bsample[:mtindex] = 1:size(bsample,1)     # Add an identifier
datamia=stack(bsample,[2:3])              # Reshape the data so we can use colors!

plot2=plot(datamia,x="Val",y="value",color="variable",Geom.line,Guide.title("Expected Value Function"))
draw(PNG("valueFunc.png", 24cm, 12cm), plot2)

 
cd(actualDir)

end
