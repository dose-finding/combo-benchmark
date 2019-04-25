# This code can be used to reproduce the results in Table 4 (Section 4.2) of the paper
# ``A benchmark for dose finding studies with unknown ordering''
# by Mozgunov, Paoletti, Jaki (2019)

# The code computes the benchmark for partial ordering for binary toxicity endpoints and continuous efficacy endpoint
# with 4 doses of drug A and 2 doses of drug B, and the monotonicity assumption satisfied within each agent.
 
# Defining the orderings
orderings<-mat.or.vec(14,8)
orderings[1,]<-c(1,2,3,4,5,6,7,8)
orderings[2,]<-c(1,2,3,5,4,6,7,8)
orderings[3,]<-c(1,2,5,3,4,6,7,8)
orderings[4,]<-c(1,5,2,3,4,6,7,8)
orderings[5,]<-c(1,2,3,5,6,4,7,8)
orderings[6,]<-c(1,2,5,3,6,4,7,8)
orderings[7,]<-c(1,5,2,3,6,4,7,8)
orderings[8,]<-c(1,2,5,6,3,4,7,8)
orderings[9,]<-c(1,5,2,6,3,4,7,8)
orderings[10,]<-c(1,2,3,5,6,7,4,8)
orderings[11,]<-c(1,2,5,3,6,7,4,8)
orderings[12,]<-c(1,5,2,3,6,7,4,8)
orderings[13,]<-c(1,2,5,6,3,7,4,8)
orderings[14,]<-c(1,5,2,6,3,7,4,8)

#Defining the weights as defined in Equation (6)
powers<-mat.or.vec(nrow(orderings),ncol(orderings))
powers[which(orderings==1)]<-0
powers[which(orderings==2)]<-1
powers[which(orderings==3)]<-2
powers[which(orderings==4)]<-3
powers[which(orderings==5)]<-3
powers[which(orderings==6)]<-2
powers[which(orderings==7)]<-1
powers[which(orderings==8)]<-0
powers

# Defining the function to order a sequence $x$ with respect to order $z$
reorder<-function(x,z){
  y<-mat.or.vec(length(x),1)
  for (i in 1:length(x)){
    y[z[i]]<-sort(x)[i]
  }
  return(y)
}

# Defining the function for the normal log likelihood
normal.log.likelihood<-function(x,mu,sigma){
  y1<-sum(x-mu)^2
  y2<-(-length(x)/2)*log(2*pi)+(-length(x)/2)*log(n)
  y<-y2 + (-1)/(2*sigma)*y1
  return(y)
}


# Setting up the parameters of the trial
n<-72 # number of patients
nsims<-100 # number of simulations (takes 10 seconds on average)
normalisation<-1/7000 # normalisation constant for log normal likelihood
sigma<-1 # s.d. of the continuous efficacy response



# Choose scenario

# Scenario 1
true.tox<-c(0.01,0.10,0.40,0.50,0.05,0.15,0.45,0.55)
true.eff<-c(0.5,0.0,-1.5,-2.5,-1.5,-2.0,-3.5,-4.5)

# Scenario 2
# true.tox<-c(0.01,0.05,0.15,0.45,0.45,0.50,0.60,0.90)
# true.eff<-c(0.0,-0.5,-3.5,-5.5,-1.0,-1.5,-4.5,-6.5)

# Scenario 3
# true.tox<-c(0.01,0.15,0.40,0.50,0.05,0.20,0.45,0.55)
# true.eff<-c(0.0,-2.0,-2.0,-2.0,0.0,-2.0,-2.0,-2.0)




# Creating the matrices
result.tox<-mat.or.vec(n,length(true.tox))
result.eff<-mat.or.vec(n,length(true.eff))
selection<-mat.or.vec(nsims,length(true.tox))
distance.all.1<-distance.all.2<-mat.or.vec(nrow(orderings),1)
like<-mat.or.vec(length(true.eff),1)


# Running the PO-Benchmark

set.seed(100)
for (z in 1:nsims){
  # Generating toxicity and efficacy profiles and corresponding responses
  for (i in 1:n){
    eff.profile<-runif(1)
    tox.profile<-runif(1)
    result.tox[i,]<-as.numeric(tox.profile<true.tox)
    result.eff[i,]<-sqrt(sigma)*qnorm(eff.profile)+true.eff
  }
  tox.est<-colSums(result.tox) # storing the number of toxicity outcomes
  
  # Compute the probability of the number of toxicities to be generated under each particular ordering 
  # as defined in Equation (4)
  for (u in 1:nrow(orderings)){
    y1<-prod((dbinom(tox.est,n,reorder(true.tox,orderings[u,])))^(1/(1+powers[u,])))
    distance.all.1[u]<-y1
  }
  probabilities.1<-distance.all.1/sum(distance.all.1)

  # Compute the probability of efficacy outcomes to be generated under each particular ordering 
  # as defined in Equation (8)
  for (u in 1:nrow(orderings)){
    for (q in 1:length(like)){
      like[q]<-normal.log.likelihood(result.eff[,q],mu=-reorder(-true.eff,orderings[u,])[q],sigma=sigma)
    }
    y2<-sum(like)
    distance.all.2[u]<-y2
  }
  probabilities.2<-exp(distance.all.2*normalisation)
  SUM<-sum(probabilities.2)  
  probabilities.2<-probabilities.2/SUM

  # Computing the probabilities of all possible combinations of orderings (as defined after Equation (8))
  all.probabilities<-probabilities.1  %*% t(probabilities.2)
  all.probabilities<-all.probabilities/sum(all.probabilities)

  # Finding the target combination under each combination of orderings
  for (u in 1:nrow(orderings)){
    for (v in 1:nrow(orderings)){
      tox.est.new<-reorder(tox.est/n,orderings[u,])
      eff.est.new<-(-reorder(-colMeans(result.eff),orderings[v,]))
      eff.est.new[which(tox.est.new>0.30 | eff.est.new>0)]<-0
      target.dose<-which(eff.est.new==min(eff.est.new))
      selection[z,target.dose]<-selection[z,target.dose]+all.probabilities[u,v]*(1/(length(target.dose)))
    }
  }

} 
new.benchmark<-colMeans(selection) # PCS by the novel benchmark


#Conducting the original benchmark
nsims<-10000
selection<-mat.or.vec(nsims,length(true.tox))
for (z in 1:nsims){
  result.tox<-mat.or.vec(n,length(true.tox))
  for (i in 1:n){
    eff.profile<-runif(1)
    tox.profile<-runif(1)
    result.tox[i,]<-as.numeric(tox.profile<true.tox)
    result.eff[i,]<-sqrt(sigma)*qnorm(eff.profile)+true.eff
  }
  
  tox.est<-colMeans(result.tox)
  eff.est<-colMeans(result.eff)
  eff.est[which(tox.est>0.30 | eff.est>0)]<-0
  target.dose<-which(eff.est==min(eff.est))
  selection[z,target.dose]<-1/length(target.dose)
}
original.benchmark<-colMeans(selection)

# Results

# The original benchmark
round(new.benchmark,2)

# The proposed PO-Benchmark
round(original.benchmark,2)


