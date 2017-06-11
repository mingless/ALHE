library("GA")
library("cec2017")

desUniformInit <- function(fitnessFcn,nvars,control)
{
  population = numeric()
  for(i in 1:nvars){
    population = cbind(population,runif(control$popSize,control$lowerBounds[i],control$upperBounds[i]))
  }
  return(population)
}


#Main des function.
#Minimizess the given fitnessFcn, which should be defined to take
# an nvars-sized vector as its argument.

#initFcn - function used to initialize the population. It is called with
# parameters (fitnessFcn,nvars,control)
#The function checks for following parameters in the control list:
#scalingCoeff - scaling coefficient. Default = 1/sqrt(2)
#popSize - population size. Default = 100*nvars
#eliteSize - the number of individuals that take part in defining the
# reference point for the next population
#migrationCoeff - coefficient in range <0.,1.> determining how fast
# the center of population moves. default = 0.5
#horizonSize - amount of the past populations that take part in creating
# the offspring. default = 10
#noiseCoeff - coefficient specifying the amount of noise added to each new individual
# default = 1
#lowerBounds - lower bounds of the initial population if using the default initFcn.
# default = rep(-1000,nvars)
#upperBounds - upper bounds --||--. Default = rep(1000,nvars)
#maxIterations - number of iterations after which optimizations stops. Currently
# the only stop condition. default = 1000*nvars
#
#TODO: Better param checks. Currently assumes that the user knows what he's doing.
#TODO: Other stop conditions
#
des <- function(fitnessFcn,nvars,initFcn=desUniformInit,control=list())
{
  if(length(control$scalingCoeff)==0)
    control$scalingCoeff = 0.70710678118
  if(length(control$popSize)==0)
    control$popSize = nvars*4
  if(length(control$eliteSize)==0)
    control$eliteSize = nvars*30
  if(length(control$migrationCoeff)==0)
    control$migrationCoeff=0.5
  if(length(control$horizonSize)==0)
    control$horizonSize = 10
  if(length(control$noiseCoeff)==0)
    control$noiseCoeff = 1
  if(length(control$lowerBounds)==0)
    control$lowerBounds = rep(-1000,nvars)
  if(length(control$upperBounds)==0)
    control$upperBounds = rep(1000,nvars)
  if(length(control$maxIterations)==0)
    control$maxIterations = 1000*nvars
  #current iteration and work index
  t = 1
  i = 1
  #to save memory past scores and populations are saved as a psedo queue
  #of size up to con$horizonSize
  Del = list()
  score = list()
  pop = list()
  bestScore = Inf
  Del[[1]] = rep(0,nvars)
  del = numeric()
  for(it in 1:nvars)
    del[it] = (sqrt(nvars) + (nvars/sqrt(nvars)))/2
  pop[[1]] = initFcn(fitnessFcn,nvars,control)
  repeat
  {
    #apply the fitness function to each row of current population
    score[[i]] = apply(pop[[i]],1,fitnessFcn)
    #mean individual from a complete population
    meanInd = colMeans(pop[[i]])
    #get the indices of sorted scores and then reorder them and the population
    ord = order(score[[i]])
    score[[i]]=score[[i]][ord]
    pop[[i]]=pop[[i]][ord,]
    #mean individual from the elite
    eliteMeanInd = colMeans(pop[[i]][1:control$eliteSize,])

    Del[[i+1]] = (1-control$migrationCoeff)*Del[[i]]+control$migrationCoeff*(eliteMeanInd-meanInd)


    #update best
    if(score[[i]][1]<bestScore)
    {
      bestScore = score[[i]][1]
      bestInd = pop[[i]][1,]
    }

    #break condition
    if(t >= control$maxIterations)
      break

    #vectors of uniformly distributed random integer variables
    h = sample(1:i,control$popSize,replace=TRUE)
    j = sample(1:control$eliteSize,control$popSize,replace=TRUE)
    k = sample(1:control$eliteSize,control$popSize,replace=TRUE)

    #new pop init
    pop[[i+1]] = matrix(0,control$popSize,nvars)
    #calculating the new individuals
    for(it in 1:control$popSize)
    {
      d <- control$scalingCoeff*(pop[[h[it]]][j[it],]-pop[[h[it]]][k[it],]) +
        Del[[i+1]]*del*rnorm(1,0,1)
      pop[[i+1]][it,] = eliteMeanInd + d + control$noiseCoeff*rnorm(nvars,0,1)
    }
    t = t+1

    #work index doesn't increase if we've reached horizonSize iterations
    #instead we remove old populations
    if(t>=control$horizonSize)
    {
      pop = pop[-1] #remove the oldest population from history
      Del = Del[-1]
      score = score[-1]
    }
    else
    {
      i = i+1;
    }
  }
  #prepare the result structure
  result = list()
  result$x = bestInd
  result$score = bestScore
  result$it = t
  result$flag = 1
  return(result)
}


desfun <- function(x) cec2017::cec2017(1,x);
gafun <- function(x) {
  horizon <- as.integer(x[3])
  elite <- as.integer(x[4])
  con=list(lowerBounds=rep(-100,10),upperBounds=rep(100,10),migrationCoeff=x[1],
           noiseCoeff=x[2],eliteSize=elite,horizonSize=horizon,maxIterations=100000)
  best <- Inf
  ressum <- 0
  for(i in 1:1)
  {
    result <- des(desfun,10,initFcn=desUniformInit,con)$score
    print(result)
    ressum <- ressum + result
    if(best > result) best <- result
  }
  resmean <- ressum/21
  return( -(resmean + (best-resmean)/4)) #negative because ga maximizes the function

}

res <- ga(type="real-valued",fitness=gafun,min=c(0,0,2,2),max=c(1,5,100,40),monitor=TRUE,maxiter=1)

