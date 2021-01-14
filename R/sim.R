get_parents <- function(weights){
  # sample from the parents based on weights
  # should fail if weights are not all >= 0
  # returns a vector of indexes into the genotype matrix representing parents
  N = length(weights)
  parents = sample.int(n=N, size=N, replace=TRUE, prob=weights)
  return(parents)
}

reproduce <- function(genos, parents1, parents2){
  # produces genotypes for one generation of reproduction
  # takes as input genotypes from the parental generation
  # and parents for each offspring

  # infer the number of loci and offspring
  L = dim(genos)[[2]] # number of loci (columns)
  N = dim(genos)[[1]] # number of individuals (rows)
  stopifnot(length(parents1) == length(parents2))
  stopifnot(N == length(parents1))

  # construct matrices of the parental genotypes
  parents_geno1 = genos[parents1, ]
  parents_geno2 = genos[parents2, ]

  # offspring are a combination of the two parents
  # binomial sampling from allele freq of each parent
  offspring <- rbinom(n=N*L, size=1, p=parents_geno1 / 2) +
      rbinom(n = N*L, size=1, p=parents_geno2 / 2)
  # convert offspring to a matrix
  offspring <- matrix(data=offspring, nrow=N, ncol=L)
  return(offspring)
}

get_initial_genotypes <- function(N, L, weights){
  # create an initial set of genotypes
  # starts from HW proportions at AF=50% and then has one generation of reproduction
  genos0 = sample(c(0,1,2), size=N * L, replace=TRUE, prob=c(1,2,1))
  genos0 = matrix(data=genos0, nrow=N, ncol=L)
  parents1 <- get_parents(weights)
  parents2 <- get_parents(weights)
  return(reproduce(genos0, parents1, parents2))
}

get_tempF <- function(P1, P2){
  # computes Nei & Tajima 1981 version of temporal F
  # takes the (unweighted) mean F across loci
  # P1 and P2 are vectors of allele frequencies from two different points in time
  stopifnot(length(P1) == length(P2))
  PP = (P1+P2)/2
  top = (P1-P2)^2
  bottom = PP - P1*P2
  tempF = mean(top/bottom)
  return(tempF)
}

get_inbreeding_Ne <- function(parents1, parents2){
  # calculate inbreeding Ne based on the number of offspring from each parent
  N = length(parents1)
  bigparents = c(parents1,parents2)
  # count number of offspring for each parent with >0 offspring
  K = as.data.frame(table(bigparents))$Freq
  sumk = sum(K)
  sumk2 = sum(K^2)
  inbreeding_Ne = (sumk-1)/((sumk2/sumk)-1)
  return(inbreeding_Ne)
}

get_expected_Ne <- function(weights){
  # calculate expected inbreeding Ne based on parental weights
  N = length(weights)
  A = mean(weights)
  B = var(weights)*(N-1)/N
  CVsq = B/A^2
  expected_Ne = N/(1+CVsq)
  return(expected_Ne)
}

get_expected_F <- function(weights, G){
  # calculate expected F based on parental weights
  # G is the number of generations
  # assumes the set weights are constant across generations
  expected_Ne = get_expected_Ne(weights)
  expected_F = 1 - (1 - 1/(2*expected_Ne))^(G)
  return(expected_F)
}

sim <- function(N, L, G, weights, Nrep=1){
  # simulate reproduction and genotypes
  # allows specification of a vector of parental weights
  # the weights vercor is permuted and reused each generation
  # reports observed and expected inbreeding Ne each generation
  # reports observed and expected temporal F each generation vs first generation
  # allows replication with Nrep
  # N = number of (diploid) individuals
  # L = number of loci
  # G = number of generations
  # weights = a vector or parental weights
  # Nrep = performs this many replicate simulations
  # output has have one row per generation per replicate
  output <- matrix(data = NA, nrow=G*Nrep, ncol= 6)
  index = 1
  for (rep in 1:Nrep){
    shuffle_weights = sample(weights) # permute weights
    parental_genotypes = get_initial_genotypes(N, L, shuffle_weights)
    initial_frequencies = colMeans(parental_genotypes)/2
    for (g in 1:G){
      shuffle_weights = sample(weights)
      parents1 = get_parents(shuffle_weights)
      parents2 = get_parents(shuffle_weights)
      offspring_genotypes = reproduce(parental_genotypes, parents1, parents2)
      offspring_freqs = colMeans(offspring_genotypes)/2
      expected_Ne = get_expected_Ne(shuffle_weights)
      expected_F = get_expected_F(shuffle_weights, G=g)
      realized_Ne = get_inbreeding_Ne(parents1, parents2)
      realized_F = get_tempF(initial_frequencies, offspring_freqs)

      # prepare for next generation and fill out results matrix
      output[index,1] = rep
      output[index,2] = g
      output[index,3] = expected_Ne
      output[index,4] = realized_Ne
      output[index,5] = expected_F
      output[index,6] = realized_F
      index = index+1
      parental_genotypes = offspring_genotypes
    }
  }
  dimnames(output) <- list(NULL, c('rep', 'g', 'expected_Ne','realized_Ne', 'expected_F',
     'realized_F' ))
  return(output)
}

get_gen_means <- function(res){
  # take mean at each generation over replicates
  mean_res = aggregate(res[, c('expected_Ne','realized_Ne', 'expected_F',
      'realized_F')], list(res[,'g']), mean)
  return(mean_res)
}

plot_tempF <- function(mean_res, title=NULL){
  # plot observed vs expected temporal F
  plot(mean_res$Group.1, mean_res$realized_F, type="p", lwd=3, col="red",
       xlab="Generation", ylab="Temporal F")
  lines(mean_res$Group.1, mean_res$expected_F,type="l",lty=3,lwd=3,col="blue")
  legend(0.15*max(mean_res$Group.1), 0.8*max(mean_res$realized_F),
         legend=c("Observed","Expected"),
         col=c("red", "blue"), lty=1:2, cex=1.5)
  title(title)
}

compare_Ne <- function(res)  {
  I = 1/res[,"realized_Ne"]
  HMNe = 1/mean(I)
  ENe = res[1,"expected_Ne"]
  Obsexp = HMNe/ENe
  col1 = c("ObsNe","ExpNe","Obs/Exp")
  col2 = c(round(HMNe, 3), round(ENe, 3), round(Obsexp, 3))
  ResultsNe = noquote(cbind(col1,col2))
  dimnames(ResultsNe) = NULL
  return(ResultsNe)
}


## choose parental weights; some possible weighting schemes:
# average over the replicates and plot temporal F

## Wright-Fisher equal probability
res_WF = sim(N=100, L=100, G=20, weights=rep(1,100), Nrep = 50)
print(compare_Ne(res_WF))
plot_tempF(get_gen_means(res_WF), title = 'WF')

## linear weights
res_seq = sim(N=100, L=100, G=20, weights=1:100, Nrep = 50)
print(compare_Ne(res_seq))
plot_tempF(get_gen_means(res_seq), title = '1:N')

## inverse weights
res_Inv = sim(N=100, L=100, G=20, weights=1/1:100, Nrep = 50)
print(compare_Ne(res_Inv))
plot_tempF(get_gen_means(res_Inv), title = '1 / 1:N')
