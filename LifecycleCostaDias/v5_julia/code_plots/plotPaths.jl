function plotPaths( )


bsample = DataFrame(t=vec(1:T+1),Sim1=vec(apath[:,1]),Sim2=vec(apath[:,2]))
bsample[:mtindex] = 1:size(bsample,1)     # Add an identifier
datamia=stack(bsample,[2:3])              # Reshape the data so we can use colors!

plot1=plot(datamia,x="t",y="value",color="variable",Geom.line,Guide.title("Time path of assets"))
draw(PNG("output\\images\\pathAssets.png", 24cm, 12cm), plot1)

# -------------------------------------------------------------------------

bsample = DataFrame(t=vec(1:T),C1=vec(cpath[:,1]),Y1=vec(ypath[:,1]))
bsample[:mtindex] = 1:size(bsample,1)     # Add an identifier
datamia=stack(bsample,[2:3])              # Reshape the data so we can use colors!

plot2=plot(datamia,x="t",y="value",color="variable",Geom.line,Guide.title("Time path of income and consumption Individual 1"))
draw(PNG("output\\images\\pathIndiv1.png", 24cm, 12cm), plot2)

# -------------------------------------------------------------------------

bsample = DataFrame(t=vec(1:T),C2=vec(cpath[:,2]),Y2=vec(ypath[:,2]))
bsample[:mtindex] = 1:size(bsample,1)     # Add an identifier
datamia=stack(bsample,[2:3])              # Reshape the data so we can use colors!

plot3=plot(datamia,x="t",y="value",color="variable",Geom.line,Guide.title("Time path of income and consumption Individual 2"))
draw(PNG("output\\images\\pathIndiv2.png", 24cm, 12cm), plot3)

end


