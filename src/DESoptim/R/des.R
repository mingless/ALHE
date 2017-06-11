#'Uniformly distributed initial population
#'@description
#'Returns initial populations for the des algorithm using
#'uniform distribution.
#'
#'@param fitnessFcn Des algorithm fitness function. Not used in this
#' init, but might be useful in a custom function
#'@param nvars Amount of variables in the optimized problem
#'@param control A control list from the des algorithm. Important
#' elements are control$upperBounds and control$lowerBounds,
#' specyfying the initial search space, and control$popSize defining
#' the amount of generated individuals.
#'
#'@return A numeric matrix with dimensions (control$popSize x nvars),
#' with each row containing subsequent variables of each individual
desUniformInit <- function(fitnessFcn,nvars,control)
{
  population = numeric()
  for(i in 1:nvars){
    population = cbind(population,runif(control$popSize,control$lowerBounds[i],control$upperBounds[i]))
  }
  return(population)
}


#'Main des function.
#'@description
#'Minimizes the given function.
#'
#'@param fitnessFcn Optimized function. Should be defined to take
#' an nvars-sized vector as it's argument.
#'@param nvars Amount of variables in the optimized function.
#'@param initFcn Function used to initialize the population. It is called with
#' parameters (fitnessFcn,nvars,control)
#'@param outputFcn function that takes the current population and it's score
#' as parameters. Can be used to plot the progress. Default = NULL
#'@param control List containing the algorithm settings. Check details for more information.
#'
#'@details
#'The function is not safe against using wrong parameters. Make sure
#'to use reasonable values, or the results might unpredictable.
#'
#'Most of the algorithm parameters is defined in the control list:
#'\itemize{
#'\item scalingCoeff - scaling coefficient. Default = 1/sqrt(2)
#'\item popSize - population size. Default = 4*nvars
#'\item eliteSize - the number of individuals that take part in defining the
#' reference point for the next population. Default = popSize/2
#'\item migrationCoeff - coefficient in range <0.,1.> determining how fast
#' the center of population changes. default = 0.5
#'\item horizonSize - amount of the past populations that take part in creating
#' the offspring. default = 10
#'\item noiseCoeff - coefficient specifying the amount of noise added to each new individual
#' default = 1e-4
#'\item lowerBounds - lower bounds of the initial population if using the default initFcn.
#' default = rep(-100,nvars)
#'\item upperBounds - upper bounds --||--. Default = rep(100,nvars)
#'\item maxIterations - number of iterations after which optimizations stops. Currently
#' the only stop condition. default = 1000*nvars
#'\item stopScore - if not NULL, stops further evaluations at given score.
#' Useful for benchmark problems, where the global minimum is already known.
#' Default - NULL
#'}
#'
#'@return A list containing the best individual, his score, number of
#' iterations passed and an exit flag (currently not used)
#'
#'@examples
#'#rosenbrock function with global minimum f(1,1)=0
#'rosen <- function(x) (1-x[1])^2+100*(x[2]-x[1]^2)^2
#'con <- list(upperBounds=c(10,10),lowerBounds=c(-10,-10),maxIterations=4000)
#'des(rosen,2,control=con)
des <- function(fitnessFcn,nvars,initFcn=desUniformInit,outputFcn=NULL,control=list())
{
  if(length(control$scalingCoeff)==0)
    control$scalingCoeff = 0.70710678118
  if(length(control$popSize)==0)
    control$popSize = nvars*4
  if(length(control$eliteSize)==0)
    control$eliteSize = control$popSize / 2
  if(length(control$migrationCoeff)==0)
    control$migrationCoeff=0.5
  if(length(control$horizonSize)==0)
    control$horizonSize = 10
  if(length(control$noiseCoeff)==0)
    control$noiseCoeff = 1e-4
  if(length(control$lowerBounds)==0)
    control$lowerBounds = rep(-100,nvars)
  if(length(control$upperBounds)==0)
    control$upperBounds = rep(100,nvars)
  if(length(control$maxIterations)==0)
    control$maxIterations = 1000*nvars
  #current iteration and work index
  t = 1
  i = 1
  #to save memory past scores and populations are saved as a pseudo queue
  #of size up to con$horizonSize

  #not really optimal most likely
  Del = list()
  score = list()
  pop = list()
  bestScore = Inf
  Del[[1]] = rep(0,nvars)
  del = numeric()
  for(it in 1:nvars)
    del[it] = sqrt(nvars)
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

    if(length(outputFcn)!=0)
      outputFcn(pop[[i]],score[[i]])

    #break condition
    if(t >= control$maxIterations)
      break
    if(length(control$stopScore)!=0)
      if(bestScore <= control$stopScore)
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

