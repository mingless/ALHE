#simple test of des
michal <- function(xx, m=10)
{
  ii <- c(1:length(xx))
  sum <- sum(sin(xx) * (sin(ii*xx^2/pi))^(2*m))

  y <- -sum
  return(y)
}

#create the grid for contour plot
x <- seq(0,pi,length.out = 51)
y <- seq(0,pi,length.out = 51)
z <- outer(x,y,Vectorize(function(x,y)michal(c(x,y))))

#plot colours
my.cols = c("#9E0142","#D53E4F","#F46D43","#FDAE61","#FEE08B","#FFFFBF",
            "#E6F598","#ABDDA4","#66C2A5","#3288BD","#5E4FA2")

#plot the function
contour(x,y,z,col=my.cols)

#make a list with control parameteres for the des function
#feel free to experiment with these to see how different parameters
#change the convergence
con <- list(lowerBounds=c(0,0), upperBounds=c(pi,pi), maxIterations=100,
            noiseCoeff=1e-1, horizonSize=12, migrationCoeff=0.6)

#create the output function to let us see the populations
totalpop <- numeric()
out <- function(pop,scores)
{
  totalpop <<- rbind(totalpop,pop)
}

#optimize the function
result <- des(michal,2,outputFcn=out,control=con)
print(result)

#plot all the generated individuals
points(totalpop[,1],totalpop[,2],cex=.4,pch=19)
