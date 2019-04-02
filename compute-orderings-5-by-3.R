nsims<-1000000
y.store<-mat.or.vec(nsims,15)
orders.store<-mat.or.vec(nsims,15)
y<-mat.or.vec(15,1)
y[1]<-0
y[15]<-1

set.seed(100)
for (i in 1:nsims){
y[2]<-runif(1,min=y[1])
y[3]<-runif(1,min=y[2])
y[4]<-runif(1,min=y[3])
y[5]<-runif(1,min=y[4])
y[6]<-runif(1,min=y[1])
y[7]<-runif(1,min=max(y[6],y[2]))
y[8]<-runif(1,min=max(y[7],y[3]))
y[9]<-runif(1,min=max(y[4],y[8]))
y[10]<-runif(1,min=max(y[9],y[5]))
y[11]<-runif(1,min=y[6])
y[12]<-runif(1,min=max(y[11],y[7]))
y[13]<-runif(1,min=max(y[12],y[8]))
y[14]<-runif(1,min=max(y[9],y[13]))
y.store[i,]<-y
}


for (i in 1:nsims){
  orders.store[i,]<-order(y.store[i,]) 
}

all.orderings<-unique(orders.store)
nrow(all.orderings)
