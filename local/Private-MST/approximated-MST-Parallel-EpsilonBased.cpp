#include <Rcpp.h>

//[[Rcpp::depends(BH)]]
// [[Rcpp::plugins(cpp11)]]
#include <iostream>
#include <string>
#include <sstream>
#include <stdlib.h>
#include <limits.h>
#include <random>

#include <boost/graph/adjacency_list.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/utility.hpp>
#include <approximated-MST-header-V4.h>

using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::NumericMatrix PrivateMST(int order, Rcpp::DataFrame Elist, double privacydegree){

  Rcpp::NumericMatrix Edgelist = DFtoNM( Elist); // transform the edgelist (dataframe type) to a
  //NumericMatrix type edgelist  for C++ to understand it

  int numberofnode = order-1;
  int numberofedge =  Edgelist.nrow();

  Graph g(order) ;

  for (int i=0;i<numberofedge;++i) {// filling the weights and creating the edges

    NodeID u = Edgelist(i,0)-1; // setting the nodes that the edge descriptor will use
    NodeID v = Edgelist(i,1)-1;
    EdgeID edge; // initializing an edge descriptor

    bool ok;
    boost::tie(edge, ok) = boost::add_edge(u,v, g); // boost::add_edge return std::pair<EdgeID,bool>,tie do it for us

    if(Edgelist(i,2)!=0){ // if the weight is not 0 fill the weight
      if (ok){  // also if the edge creation worked
        g[edge].weight = Edgelist(i,2);
      }
    }
  }

  std::cout<< "Remplissage du graph OK" <<std::endl;
  std::cout<<"-------------------------"<<std::endl;


  Rcpp::NumericMatrix mstexpmechanism=approximated_MST(g,(privacydegree*numberofedge)/(2*(numberofnode)));
  //call the c++ function to approximate the MST using prim+exponential mechanism,
  //the privacy parameter here is the parameter that the exponential distribution
  //will use therefore it is epsilon'/2*delta_u (remeber that we use epsilon'=epsilon/(Nmber of nodes-1))
  return(mstexpmechanism);
}

/*** R
#install.packages('igraph')
#install.packages('foreach')
#install.packages('itertools')
#source("https://bioconductor.org/biocLite.R")

library('igraph')
  library('RBGL')
  library('foreach')
  library('itertools')
  library('abind')
  library('doMC')
  library(matrixStats)
  registerDoMC(cores=20)


  dense=function(x,y){
    return(abs(x-y))
  }

simlap=function(n=1,loc=0,scale=1){
  u=runif(n)
  return(loc-scale*(log(2*u)*(u<1/2)-log(2*(1-u))*(u>1/2)))
}
sup_q=function(x,q){
  return(as.numeric(x>q))
}

n=1000
niter=100
lengthepsilon=3
lengthprob=9
epsiloniter=c(0.1,0.5,1)
probiter=seq(0.1,0.9,length=lengthprob)
alpha=0.05
q=qnorm(1-alpha/2)

for (epsilon in epsiloniter){
  compare_tensor<- foreach(prob=probiter,.combine = abind) %:%
foreach ( j=1:niter , .combine='rbind')%dopar%{

#epsilon=1
  graph = erdos.renyi.game(n, prob, type="gnp",directed = FALSE, loops = FALSE)
  E(graph)$weight=rnorm(ecount(graph),0,1)
  graphlap=graph
  E(graphlap)$weight=E(graph)$weight + simlap(ecount(graphlap),scale=1/epsilon)


  g=igraph.to.graphNEL(graphlap) # turning to graph class from boost
  gTrue=igraph.to.graphNEL(graph)
  MST=mstree.kruskal(g) #Computing mst kruskal using boost
  MSTTrue=mstree.kruskal(gTrue)
  Truecost=sum(MSTTrue$weights)


  edgelist=data.frame("from"=as.numeric(MST$edgeList[1,]),"to"=as.numeric(MST$edgeList[2,]))
  Dataframe=igraph::as_data_frame(graph,what="edges")
  Dataframe2=rbind(Dataframe,data.frame("from"=c(Dataframe[,2]),"to"=c(Dataframe[,1]),"weight"=c(Dataframe[,3])))
  Privatemst=PrivateMST(n,Dataframe,epsilon)

  Laplace=sum(merge(Dataframe2,edgelist,by.x=c('from','to'),by.y=c('from','to'))$weight)
  Exponential=sum(Privatemst[,3])
  real=Truecost

  c(Laplace,Exponential,real)
}

indices=1:(lengthprob*3)

Mean=colMeans(compare_tensor)
Sd=colSds(compare_tensor)
meanLaplace=Mean[indices%%3==1]
meanExponential=Mean[indices%%3==2]
meanReal=Mean[indices%%3==0]

sdLaplace=Sd[indices%%3==1]
sdExponential=Sd[indices%%3==2]
sdReal=Sd[indices%%3==0]

titrepdf=paste("Courbe-MST-",epsilon,"-Epsilon-Based.pdf")
pdf(titrepdf, height=10,width=10)

#on trace le graphique
titre = paste("comparison of the two Private MST techniques for ",n ,"nodes \n and Epsilon: " , epsilon)
plot(epsiloniter,meanLaplace,type="o",cex=2,col="grey",ylab="weight of the MST",xlab="probability for Erdos-Renyi simulation",ylim=c(0,max(meanLaplace)+5),main=titre,pch=16)
lines(epsiloniter,meanExponential,col="red",pch=18,type="o",cex=2)

lines(epsiloniter,meanExponential+sdExponential*q/sqrt(niter),lty=2,col="red")
lines(epsiloniter,meanExponential-sdExponential*q/sqrt(niter),lty=2,col="red")

lines(epsiloniter,meanLaplace+sdLaplace*q/sqrt(niter),lty=2,col="grey")
lines(epsiloniter,meanLaplace-sdLaplace*q/sqrt(niter),lty=2,col="grey")

legend("topright",c("Laplace","PAMST"), cex=2, col=c("grey","red"), lty=c(1,1),pch=c(16,18), bty="n")
dev.off()
}


  */
