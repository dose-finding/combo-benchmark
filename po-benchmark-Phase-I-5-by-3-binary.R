# This code can be used to reproduce the results in Table 2 (Section 4.1) of the paper
# ``A benchmark for dose finding studies with unknown ordering''
# by Mozgunov, Paoletti, Jaki (2019)

# The code computes the benchmark for partial ordering for binary toxicity endpoints
# with 3 doses of drug A and 5 doses of drug B, and the monotonicity assumption
# satisfied within each agent.
 
# (!) The code requires the orderings to be specified (run file compute-orderings-5-by-3.R)
 
#Defining the weights as defined in Equation (6)
weight<-mat.or.vec(nrow(all.orderings),ncol(all.orderings))
weight[which(all.orderings==1)]<-0
weight[which(all.orderings==2)]<-2
weight[which(all.orderings==3)]<-4
weight[which(all.orderings==4)]<-6
weight[which(all.orderings==5)]<-8
weight[which(all.orderings==6)]<-4
weight[which(all.orderings==7)]<-4
weight[which(all.orderings==8)]<-4
weight[which(all.orderings==9)]<-4
weight[which(all.orderings==10)]<-4
weight[which(all.orderings==11)]<-8
weight[which(all.orderings==12)]<-6
weight[which(all.orderings==13)]<-4
weight[which(all.orderings==14)]<-2
weight[which(all.orderings==15)]<-0


# Defining the function to order a sequence $x$ with respect to order $z$
reorder<-function(x,z){ # function to reorder scenarios with respect to all possible orderings
  y<-mat.or.vec(length(x),1)
  for (i in 1:length(x)){
    y[z[i]]<-sort(x)[i]
  }
  return(y)
}

# Simulation Setting
n<-60 # Total number of Patients
target<-0.30 # Target Toxicity Level 
nsims<-20  # Number of simulations (average time for 20 simulations and 6006 orderigns is 1.5 minutes)
# Note: The results are accurate enough even for small number of simulations

# Choose scenario
# true.tox<-c(0.05,0.10,0.15,0.30,0.45,0.10,0.15,0.30,0.45,0.55,0.15,0.30,0.45,0.50,0.60) # Scenario 1
# true.tox<-c(0.15,0.30,0.45,0.50,0.60,0.30,0.45,0.50,0.60,0.75,0.45,0.55,0.60,0.70,0.80) # Scenario 2
# true.tox<-c(0.02,0.07,0.10,0.15,0.30,0.07,0.10,0.15,0.30,0.45,0.10,0.15,0.30,0.45,0.55) # Scenario 3
# true.tox<-c(0.30,0.45,0.60,0.70,0.80,0.45,0.55,0.65,0.75,0.85,0.50,0.60,0.70,0.80,0.90) # Scenario 4
# true.tox<-c(0.01,0.02,0.08,0.10,0.11,0.03,0.05,0.10,0.13,0.15,0.07,0.09,0.12,0.15,0.30) # Scenario 5
true.tox<-c(0.05,0.08,0.10,0.13,0.15,0.09,0.12,0.15,0.30,0.45,0.15,0.30,0.45,0.50,0.60) # Scenario 6
# true.tox<-c(0.07,0.10,0.12,0.15,0.30,0.15,0.30,0.45,0.52,0.60,0.30,0.50,0.60,0.65,0.75) # Scenario 7
# true.tox<-c(0.02,0.10,0.15,0.50,0.60,0.05,0.12,0.30,0.55,0.70,0.08,0.15,0.45,0.60,0.80) # Scenario 8
# true.tox<-c(0.005,0.01,0.02,0.04,0.07,0.02,0.05,0.08,0.12,0.15,0.15,0.30,0.45,0.55,0.65) # Scenario 9
# true.tox<-c(0.05,0.10,0.15,0.30,0.45,0.45,0.50,0.60,0.65,0.70,0.70,0.75,0.80,0.85,0.90) # Scenario 10

# Create matrices to store the output
result.tox<-mat.or.vec(n,length(true.tox))
selection<-mat.or.vec(nsims,length(true.tox))
distance.all<-mat.or.vec(nrow(all.orderings),1)

# Run the novel benchmark
set.seed(100)
for (z in 1:nsims){
  # generating complete information knowing the ordering
  for (i in 1:n){ # for each patient
    tox.profile<-runif(1) # generate toxicity profile
    result.tox[i,]<-as.numeric(tox.profile<true.tox) # transform the toxicity profile for each dose to 0/1
  }
  
  tox.est<-colSums(result.tox) # store number of toxicities at each combination
  
  
  # Compute the probability of the number of toxicities to be generated under each particular ordering 
  # as defined in Equation (4)
  for (u in 1:nrow(all.orderings)){
    distance.all[u]<-prod((dbinom(tox.est,n,reorder(true.tox,all.orderings[u,])))^(1/(1+weight[u,]))) #weight function as in Equation (6)
  }
  probabilities<-distance.all/sum(distance.all) # Normalising as in Equation (5)
  
  # Find the selected combination under each ordering
  for (u in 1:nrow(all.orderings)){
    tox.est.new<-reorder(tox.est/n,all.orderings[u,])
    target.dose<-which(abs(tox.est.new-target)==min(abs(tox.est.new-target)))
    selection[z,target.dose]<-selection[z,target.dose]+probabilities[u]*(1/(length(target.dose)))
  }
  
  # cat(z,"out of",nsims,"\n")
}
new.benchmark<-colMeans(selection)

# Run the original benchmark
selection<-mat.or.vec(nsims,length(true.tox))
for (z in 1:nsims){
  result.tox<-mat.or.vec(n,length(true.tox))
  for (i in 1:n){
    tox.profile<-runif(1)
    result.tox[i,]<-as.numeric(tox.profile<true.tox)
  }
  
  tox.est<-colMeans(result.tox)
  selection[z,which(abs(tox.est-target)==min(abs(tox.est-target)))]<-1/length(which(abs(tox.est-target)==min(abs(tox.est-target))))
}
original.benchmark<-colMeans(selection)


#PCS for the original benchmark
round(sum(original.benchmark[which(abs(true.tox-target)==min(abs(true.tox-target)))]),2)

# PCS for the novel benchmark
round(sum(new.benchmark[which(abs(true.tox-target)==min(abs(true.tox-target)))]),2)


