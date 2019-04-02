# install.packages("Deducer")
library("Deducer")
# install.packages("Deducer", dependencies=TRUE)
############################################################
#########          3 by 3 combination trial        #########
#########      Number of orderings computation     #########
############################################################
X<-perm(2:8)
R<-nrow(X)
C<-ncol(X)

for (i in 1:R){
  y<-X[i,]
  if(y[1]!=0){
    y1<-which(y==2) > which(y==3)
    y2<-which(y==4) > which(y==5)
    y3<-which(y==5) > which(y==6)
    y4<-which(y==7) > which(y==8)
    y5<-which(y==2) > which(y==5)
    y6<-which(y==3) > which(y==6)
    y7<-which(y==4) > which(y==7)
    y8<-which(y==5) > which(y==8)
    if(y1|y2|y3|y4|y5|y6|y7|y8){
      X[i,]<-0
    }
  }
  # cat(i,"out of",R,"\n")
}
all.orderings<-X[X[,1]!=0,]
all.orderings<-cbind(rep(1,42),all.orderings,rep(9,42))
all.orderings
nrow(all.orderings)