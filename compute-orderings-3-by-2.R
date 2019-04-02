library("Deducer")

############################################################
#########         3 by 2 combination trial         #########
#########      Number of orderings computation     #########
############################################################
X<-perm(1:6)
R<-nrow(X)
C<-ncol(X)

for (i in 1:R){
  y<-X[i,]
  if(y[1]!=0){
    y1<-which(y==1) > which(y==2)
    y2<-which(y==2) > which(y==3)
    y3<-which(y==1) > which(y==4)
    y4<-which(y==4) > which(y==5)
    y5<-which(y==5) > which(y==6)
    y6<-which(y==2) > which(y==5)
    y7<-which(y==3) > which(y==6)
    if(y1|y2|y3|y4|y5|y6|y7){
      X[i,]<-0
    }
  }
  # cat(i,"out of",R,"\n")
}
all.orderings<-X[X[,1]!=0,]



