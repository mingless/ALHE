#let's find a minimum of levy function (spoiler: f(1,1) = 0)
#when starting in a small square (x and y both between -9 and -10)
#to check effects of the noise coefficient
levy <- function(xx)
{
  d <- length(xx)
  w <- 1 + (xx - 1)/4

  term1 <- (sin(pi*w[1]))^2
  term3 <- (w[d]-1)^2 * (1+1*(sin(2*pi*w[d]))^2)

  wi <- w[1:(d-1)]
  sum <- sum((wi-1)^2 * (1+10*(sin(pi*wi+1))^2))

  y <- term1 + sum + term3
  return(y)
}

#generate a grid for the contour plot
x <- seq(-15,15,length.out = 71)
y <- seq(-15,15,length.out = 71)
z <- outer(x,y,Vectorize(function(x,y)levy(c(x,y))))

#make a list with control parameteres for the des function
con <- list(lowerBounds=c(-10,-10), upperBounds=c(-9,-9), maxIterations=4000,
           noiseCoeff=1e-1, horizonSize=6, migrationCoeff=0.9)

#create an output function to help us plot all the individuals
totalpop <- numeric()
out <- function(pop,scores)
{
  totalpop <<- rbind(totalpop,pop)
}

#try to optimize levy using the given parameters
result <- des(levy,2,outputFcn=out,control=con)
print(result)

#plot colours
my.cols = c("#9E0142","#D53E4F","#F46D43","#FDAE61","#FEE08B","#FFFFBF",
            "#E6F598","#ABDDA4","#66C2A5","#3288BD","#5E4FA2")

#plots!
contour(x,y,z,levels=as.integer(exp(-5:9)),col=my.cols)
points(totalpop[,1],totalpop[,2],cex=.4,pch=19)

#for noise = 1e-1 the function stops in a local minimum, as we can
#see on the plot

#let's change noise to 1e0
con$noiseCoeff = 1e0

#and try optimizing again
result <- des(levy,2,outputFcn=out,control=con)
print(result)

#more plots!
contour(x,y,z,levels=as.integer(exp(-5:9)),col=my.cols)
points(totalpop[,1],totalpop[,2],cex=.4,pch=19)

#function now finds the correct minimum, although it's not very precise
