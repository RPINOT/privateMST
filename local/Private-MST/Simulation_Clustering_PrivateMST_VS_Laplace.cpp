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
##################################################################
# SIMULATION ET CLUSTERING Rapport de Stage
################################################################## 


###
# Simulation
###############################################################

f1=function(x){
  n=length(x)
  return(-x*x + rnorm(n,0,15) +200)
}

f2=function(x){
  x=x
  n=length(x)
  return(x*x + rnorm(n,0,15))
}


n=100

x1=runif(n,-10,10)
  y1=f1(x1)/10
plot(x1,y1,ylim=c(-10,30),xlim=c(-10,20))
  x2=runif(n,-10,10)
  y2=f2(x2)/10
points(x2+10,y2)
  
  DB1=cbind(c(x1,x2+10),c(y1,y2))
  
  
  
  
###
# DBClust
###############################################################
  
#########################################
#  Some Useful Packages
#########################################
  
  
#install.packages('igraph')
#install.packages('foreach')
#install.packages('stats')
  
  library('stats')
    library('igraph')
    library('foreach')
    
#########################################
#  Simple Functions
#########################################
    
# "Disp" function checks if the igraph object has at least one edge. 
# If it does, it return the maximum weight of the edges if it doesn't it returns 0
    
    Disp=function(ClusterGraph){ # ClusterGraph is an igraph object representing a non-oriented simple weighted graph
      
      if(ecount(ClusterGraph)==0){ #check if the graph has at least one edge.
        return(0) #if not, return 0
      }else{
        return(max(E(ClusterGraph)$weight)) #if it does, return the maximum of the edges weights 
      }
    }
    
################
    
# "Sep" function checks if the dataframe is not empty. 
# It returns 1 if it is, and the minimum value of the 3rd columns if it isn't.
    
    Sep=function(Incidentremoved){ # Incidentremoved is a dataframe with 3 columns
      
      Incidentremoved=matrix(as.matrix(Incidentremoved),ncol=3) # cast the Incidentremoved into a matrix 
      if(dim(Incidentremoved)[1]==0){ # check if the matrix is empty 
        return(1) # if it is return 1
      }else{
        return(min(Incidentremoved[,3])) # else return the minimum value of the 3rd columns
      }
      
    }
    
################
    
# "ValidityIndex" returns the validity index of a list composed of a igraph object and a matrix.
# the computation is regarding the above functions (Sep(matrix)-Disp(graph))/max(Sep(matrix),Disp(graph))*NumberOfNodeInTheGraph
    
    ValidityIndex=function(Cluster){ # Cluster is a list containing, in this order, a igraph object and a matrix (which is an edgelist) )
      
      Disptemp=Disp(Cluster[[1]]) # Compute Disp(graph)
      Septemp=Sep(Cluster[[2]]) # Compute Sep(matrix)
      
      return(((Septemp-Disptemp)/max(Septemp,Disptemp))*vcount(Cluster[[1]])) # Return the Validity index 
    }
    
################
    
# "getDBCVI" returns the sum of the "ValidityIndex" of every element of a list (each one of those element is a list(graph,matrix) as seen above) 
# this sum is normalized by "Graphorder" representing the total amount of node in the graphs. "Clusternumber" is the size of the list and is here for simplicity.
    
    getDBCVI=function(Forestpartition,Graphorder,Clusternumber){ # Forestpartition is a list(list(graph,matrix)) ; Graphorder is a normalization parameter (integer) ; Clusternumber is the size of the first list in Forestpatition
      
      DBCVI<-foreach(j=1:Clusternumber, .combine = '+')%do%{ # .combine ='+' is going to sum the on-going values of Validity index
        ValidityIndex(Forestpartition[[j]]) # for every Graph listed, compute its ValidityIndex.s
      }
      return(DBCVI/Graphorder) # return the value of the sum normalized by the total number of nodes in the graphs
    }
    
    
#########################################
#  Perform and Cut 
######################################### 
    
#'performCut' takes a list(list(graph,matrix)), a number indicating the location of a particular graph in the list and an edge belonging to the graph.
# It returns an update of the list(list(graph,matrix)) where the particular graph is demerged into two graphs.
# The matrix for each graph is updated according to the incidence of the edges it contains to one or the other graph. 
    
    performCut=function(Forestpartition,cut,partitionnumber=1,currentnumberofcluster=1){ # Forestpartition is a list(list(graph,matrix)); cut is a 1x3 (vertex,vertex,weight) vector;  partitionnumber is an int corresponding to the position of the graph we want to demerge; currentnumberofnode is the size of Forestpartition and is here for convinience
      
      ClusterGraph=Forestpartition[[partitionnumber]][[1]] # storing the chosen igraph object in a variable ClusterGraph
      
      edge=paste(cut[1],"|",cut[2], sep="") # the edge we want to cut has to be written like that to be understoof by igraph
      ClusterGraphtemp=igraph::delete_edges(ClusterGraph,edge) # deleting the edge with the igraph function
      ctemp=igraph::clusters(ClusterGraphtemp) # getting the connexes component from the updated graph (there is two of them because the graph is actually a tree)
      
      indicecluster1=names(which(ctemp$membership==1)) #getting the names of the nodes in each cluster
      indicecluster2=names(which(ctemp$membership==2)) #getting the names of the nodes in each cluster
      
      Ref=rbind(Forestpartition[[partitionnumber]][[2]],cut) # adding the deleted edge to the incident deleted edges stored in the matrix 
      
      incidentcutscluster1=Ref[Ref[,1] %in% indicecluster1 | Ref[,2] %in% indicecluster1,] # splitting the matrix into two matrix according 
      incidentcutscluster2=Ref[Ref[,1] %in% as.numeric(indicecluster2) | Ref[,2] %in% indicecluster2,] # to the incidence of each edge to one of the two clusters
      
      Forestpartition[[partitionnumber]]=list(igraph::induced_subgraph(ClusterGraph, indicecluster1),incidentcutscluster1) # update the old value of list(graph,matrix) to one of the new ones
      Forestpartition[[currentnumberofcluster+1]]=list(igraph::induced_subgraph(ClusterGraph, indicecluster2),incidentcutscluster2) # add a new element to the list which is the second couple list(graph,matrix) from the splitting
      
      return(Forestpartition) # return the updated list(list(graph,matrix))
    }
    
#########################################
#  Body of the Algorithm
#########################################
    
#'ClusteringAlgorithm' is a function that takes a graph ( in fact it is a tree ) and return a 
# partition subtrees that form clusters according to a smart splitting explained in the associated paper.
    
    
    
    ClusteringAlgorithm=function(MSTigraph){ # MSTigraph is a igraph object representing a tree
      
      V(MSTigraph)$name = V(MSTigraph) # giving names to the nodes so that a graph can be splitted while keeping track of the initiale names
      
      Forestpartition=list(list(MSTigraph,matrix(c(0,0,0),ncol=3))) # Initializing the partition 
      
      splitDBCVI = -1 #initializing the current score
      splitDBCVItemp = -1
      k=1
      N=vcount(MSTigraph) # geting the size of the graph 
      
      while(splitDBCVI < 1){ #while the value has not been chosen to break the loop
        
        OptimalCut=rep(NA,4) # at each round initialize the Optimal cut to (NA,NA,NA)
        for( i in 1:k){ #for every tree in the partition
          ClusterEdgelist=as.matrix(igraph::as_data_frame(Forestpartition[[i]][[1]])) #store the tree as an edgelist
          
          if(dim(ClusterEdgelist)[1]>0){ # if the list is not empty 
            for(numline in 1:(dim(ClusterEdgelist)[1]) ){ # line in the edgelist  
              cut=ClusterEdgelist[numline,1:3] # choose the line to represent the edge to be cut 
              newCluster=performCut(Forestpartition,cut,i,k) # get a new partition to split the current tree according to the currrent cut
              newDBCVI = getDBCVI(newCluster,N,k+1) # get the score of the new partition
              
              if(is.na(newDBCVI)){
                print("newDBCVI is NA")
              }
              else if(is.na(splitDBCVItemp)){
                print("split is NA")
              }
              else if(newDBCVI>=splitDBCVItemp){ # if this score is better than te current better score 
                OptimalCut=c(cut,i) # set the current cut to be the optimal cut (adding i is the indication of the tree the cut belonged to)
                splitDBCVItemp=newDBCVI # set the current score to become the better score
              }
            }
          }
        }
        
        if(!is.na(OptimalCut[1])){ # if the optimal cut has been modified at least once
          
          Forestpartition=performCut(Forestpartition,OptimalCut[1:3],OptimalCut[4],k) #update the real partition to be the one stemming from the splitting according to the optimale cut
          splitDBCVI = splitDBCVItemp # update current score to be the score of this update 
          k=k+1 # update the number of partitions in Forestpartition
        }else{
          splitDBCVI=1 # else set te current score to be a score that can not be obtain in a real setting to break the loop 
        }
      }
      
      return(list(k,Forestpartition)) # return the partition
    }
    
    
    
    
    
###
# Comparison of Private MSDR and None Private MSDR
###

### Functions
    simlap=function(n=1,loc=0,scale=1){
      u=runif(n)
      return(loc-scale*(log(2*u)*(u<1/2)-log(2*(1-u))*(u>1/2)))
    }
    
    dist=function(A,B){
      return((A-B)^2)
    }
    
    
### None Private DB1
    
    A=matrix(data = 0,nrow = n*2,ncol=n*2)
      A=outer(DB1[,1],DB1[,1],FUN=dist)+outer(DB1[,2],DB1[,2],FUN=dist)
      diag(A)=0
    
    G=graph_from_adjacency_matrix(A,weighted = TRUE,mode="undirected")
      MSTigraph=igraph::mst(G)
      
      P=ClusteringAlgorithm(MSTigraph)
      
      
      k=P[[1]]
    Classif=P[[2]]
    
    plot(DB1[,1],DB1[,2],ylim=c(-10,30),xlim=c(-10,20),col='black')
      for(i in 1:k){
        C=V(Classif[[i]][[1]])$name
        points(DB1[C,1],DB1[C,2],col=i+1,main="1")
      }
#dev.off()
### Private With epsilon = 1 DB1
      
      n_edges=ecount(G)
        E(G)$weight=E(G)$weight + simlap(n_edges,scale=1)
        E(G)$weight=E(G)$weight-rep(min(0,min(E(G)$weight)),length(E(G)$weight))+1
      
      MSTigraph=igraph::mst(G)
        P=ClusteringAlgorithm(MSTigraph)
        
        k=P[[1]]
      Classif=P[[2]]
      
      plot(DB1[,1],DB1[,2],ylim=c(-10,30),xlim=c(-10,20),col='black',main="2")
        for(i in 1:k){
          C=V(Classif[[i]][[1]])$name
          points(DB1[C,1],DB1[C,2],col=i)
        }
#dev.off()
### Private With epsilon = 0.1 DB1
        
        G=graph_from_adjacency_matrix(A,weighted = TRUE,mode="undirected")
          
          n_edges=ecount(G)
          E(G)$weight=E(G)$weight + simlap(n_edges,scale=10)
          E(G)$weight=E(G)$weight-rep(min(0,min(E(G)$weight)),length(E(G)$weight))+1
        
        MSTigraph=igraph::mst(G)
          P=ClusteringAlgorithm(MSTigraph)
          
          
          k=P[[1]]
        Classif=P[[2]]
        
        plot(DB1[,1],DB1[,2],ylim=c(-10,30),xlim=c(-10,20),col='black',main="3")
          for(i in 1:k){
            C=V(Classif[[i]][[1]])$name
            points(DB1[C,1],DB1[C,2],col=i)
          }
#dev.off()
### Private MST With epsilon = 1 DB1
        
        G=graph_from_adjacency_matrix(A,weighted = TRUE,mode="undirected")
        Dataframe=as_data_frame(G,what = "edges")
        MSTigraph=graph_from_data_frame(PrivateMST(2*n,Dataframe,1), directed = FALSE)
        
        P=ClusteringAlgorithm(MSTigraph)
        
        k=P[[1]]
        Classif=P[[2]]
        
        plot(DB1[,1],DB1[,2],ylim=c(-10,30),xlim=c(-10,20),col='black',main="4")
        for(i in 1:k){
          C=V(Classif[[i]][[1]])$name
          points(DB1[C,1],DB1[C,2],col=i+1,main="1")
        }
        
        ##### Private MST epsilon=0.1
        
        G=graph_from_adjacency_matrix(A,weighted = TRUE,mode="undirected")
        Dataframe=as_data_frame(G,what = "edges")
        MSTigraph=graph_from_data_frame(PrivateMST(2*n,Dataframe,0.1), directed = FALSE)
        P=ClusteringAlgorithm(MSTigraph)
        
        k=P[[1]]
        Classif=P[[2]]
        
        plot(DB1[,1],DB1[,2],ylim=c(-10,30),xlim=c(-10,20),col='black',main="5")
        for(i in 1:k){
          C=V(Classif[[i]][[1]])$name
          points(DB1[C,1],DB1[C,2],col=i+1,main="1")
        }
        
dev.off()
          
*/
