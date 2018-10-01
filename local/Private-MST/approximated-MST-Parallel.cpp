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
library('foreach')
library('itertools')
library('abind')
library('doMC')
library('matrixStats')
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
lengthepsilon=20
epsiloniter=seq(0.1,1,length=lengthepsilon)
prob=0.1
alpha=0.05
q=qnorm(1-alpha/2)

compare_tensor<- foreach(epsilon=epsiloniter,.combine = abind) %:%
  foreach ( j=1:niter , .combine='rbind')%dopar%{

    #epsilon=1
    graph = erdos.renyi.game(n, prob, type="gnp",directed = FALSE, loops = FALSE)
    E(graph)$weight=runif(ecount(graph),0,10)
    graphlap=graph
    E(graphlap)$weight=E(graphlap)$weight + simlap(ecount(graphlap),scale=1/epsilon)


    MST=mst(graphlap,E(graphlap)$weight)
    MSTTrue=mst(graph,E(graph)$weight)
    Truecost=sum(E(MSTTrue)$weight)

    edgelist=igraph::as_data_frame(MST,what="edges")
    Dataframe=igraph::as_data_frame(graph,what="edges")
    Dataframe2=rbind(Dataframe,data.frame("from"=c(Dataframe[,2]),"to"=c(Dataframe[,1]),"weight"=c(Dataframe[,3])))
    Privatemst=PrivateMST(n,Dataframe,epsilon)

    Laplace=sum(merge(Dataframe2,edgelist,by.x=c('from','to'),by.y=c('from','to'))$weight.x)-Truecost
    Privatemst=PrivateMST(n,Dataframe,epsilon)
    Exponential=sum(Privatemst[,3])-Truecost

    c(Laplace,Exponential)
  }

indices=1:(lengthepsilon*2)
Mean=colMeans(compare_tensor)
Sd=colSds(compare_tensor)
meanLaplace=Mean[indices%%2==1]
meanExponential=Mean[indices%%2==0]

sdLaplace=Sd[indices%%2==1]
sdExponential=Sd[indices%%2==0]


titrepdf=paste("Courbe-MST-Erdos",prob,"Unif(0,10).pdf")
pdf(titrepdf, height=10,width=10)

#on trace le graphique

# ATTENTION IL N4Y A PAS DE TITRE POUR CE GRAPHIQUE

write.csv(cbind(epsiloniter,meanLaplace,meanExponential,sdExponential,sdLaplace),"PAMSTComparaisonLaplace.csv")
titre = paste("Comparison of the two Private MST techniques for a graph of size O(", n*n/10 ,") \n simulated by Erdos-Renyi method of probability: " , prob)
plot(epsiloniter,meanLaplace,type="o", cex.lab=1.7,cex.axis=1.7,cex=2,col="blue",ylab="Expected approximation error",xlab="Epsilon value",ylim=c(min(meanExponential)-5,max(meanLaplace)+5),pch=17)
lines(epsiloniter,meanExponential,col="red",pch=16,type="o",cex=2)

lines(epsiloniter,meanExponential+sdExponential*q/sqrt(niter),lty=2,col="grey")
lines(epsiloniter,meanExponential-sdExponential*q/sqrt(niter),lty=2,col="grey")

lines(epsiloniter,meanLaplace+sdLaplace*q/sqrt(niter),lty=2,col="grey")
lines(epsiloniter,meanLaplace-sdLaplace*q/sqrt(niter),lty=2,col="grey")

legend("topright",c("Laplace","PAMST"), cex=2, col=c("blue","red"), lty=c(1,1),pch=c(17,16), bty="n")
dev.off()



*/
