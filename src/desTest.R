#Simple test
con = list()
con$lowerBounds=c(-5,-5)
con$upperBounds=c(5,5)
con$maxIterations=4000
rosen = function(x) (1-x[1])^2+100*(x[2]-x[1]^2)^2

result = des(rosen,2,control=con)

