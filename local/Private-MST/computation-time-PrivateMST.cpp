#include <Rcpp.h>

//[[Rcpp::depends(BH)]]
// [[Rcpp::plugins(cpp11)]]
#include <iostream>
#include <string>
#include <sstream>
#include <stdlib.h>
#include <limits.h>
#include <random>

#include <boost/graph/copy.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/utility.hpp>
#include <boost/graph/kruskal_min_spanning_tree.hpp>
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
library('matrixStats')

registerDoMC()
options(cores=10)


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


niter=10
epsilon=1
prob=0.9
alpha=0.05
q=qnorm(1-alpha/2)
#seq=c(5000)
n=8000
#CourbeComputationTime<-foreach(n=seq,.combine='cbind')%do%{

  ComputationTime<-foreach( i=1:niter,.combine='rbind')%dopar%{
    graph = erdos.renyi.game(n, prob, type="gnp",directed = FALSE, loops = FALSE)
    E(graph)$weight=rnorm(ecount(graph),0,1)
    Dataframe=igraph::as_data_frame(graph,what="edges")
    start.time <- Sys.time()
    PrivateMST(n,Dataframe,epsilon)
    end.time <- Sys.time()
    time.taken1 <- end.time - start.time

    km=igraph.to.graphNEL(graph)
    start.time <- Sys.time()
    mstree.kruskal(km)
    end.time <- Sys.time()
    time.taken2 <- end.time - start.time

    start.time <- Sys.time()
    mst(graph)
    end.time <- Sys.time()
    time.taken3 <- end.time - start.time
    c(n,time.taken1,time.taken2,time.taken3)
  }
    print(n)
    print(ComputationTime)
    CourbeComputationTime=colMeans(ComputationTime)
    print(CourbeComputationTime)
#}

write.csv(CourbeComputationTime,"Valeurcomputationtime12.csv")
#pdf("CourbeComputationTime5.pdf", height=10,width=10)
#plot(seq,CourbeComputationTime[,1],col="red",cex.main=1.5, cex=2,cex.lab=1.5,cex.axis=1.5,cex=2,main="Comparison of the computationnal for time PAMST, Kruskal and Prim", xlab="Graph size",ylab="Time in seconds:",ylim=c(-1,max(CourbeComputationTime[,2])),type="o",pch=18)
#lines(seq,CourbeComputationTime[,2],type="o",pch=16,col='grey',cex=2)
#lines(seq,CourbeComputationTime[,3],type="o",pch=17,col="green",cex=2)
#legend("topleft",c("Kruskal","PAMST","Prim"), cex=2, col=c("grey","red","green"), lty=c(1,1,1),pch=c(16,18,17), bty="n")

*/