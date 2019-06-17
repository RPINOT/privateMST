#include <Rcpp.h>

//[[Rcpp::depends(BH)]]
//[[Rcpp::plugins(cpp11)]]

#include <iostream>
#include <algorithm>
#include <vector>

#include <string>
#include <sstream>
#include <stdlib.h>
#include <limits.h>
#include <random>
#include <math.h>
#include <list>


#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/utility.hpp>
#include <boost/graph/incremental_components.hpp>
#include <boost/pending/disjoint_sets.hpp>
#include <boost/foreach.hpp>
#include <boost/graph/graph_utility.hpp>
#include <boost/graph/adjacency_list.hpp>

using namespace boost;
using namespace Rcpp;
using namespace std;


// function to transform the Data frame that R will
// provide into a NumericMatrix for c++ to understand it
Rcpp::NumericMatrix DFtoNM( Rcpp::DataFrame x) {
  int nRows=x.nrows();
  Rcpp::NumericMatrix y(nRows,x.size());
  for (int i=0; i<x.size();i++) {
    y( Rcpp::_,i) = Rcpp::NumericVector(x[i]); // y( Rcpp::_,i) is a manner of choosing only the columns of y according to i
  }
  return y;
}


struct Nodeclass{ // define the characteristics of a node (here it has only a name)
  std::string name;
};

struct Edgeclass{ // define the characteristics of an edge (here it has only a weight, given that its name will be denoted by the two nodes forming the edge)
  float weight;
};

typedef adjacency_list <vecS, vecS, undirectedS, Nodeclass, Edgeclass> Graph; // basic construction of a graph using boost::adjacency_list according to the characteristics we gave above
typedef graph_traits<Graph>::vertex_descriptor Vertex; // vertex and edges descriptors are basics from  boost::graph that will help us call the node by their name and the edges by a couple of nodes forming it
typedef graph_traits<Graph>::vertices_size_type VertexIndex;
typedef graph_traits<Graph>::edge_descriptor Edge;




Rcpp::NumericMatrix approximated_MST_with_PAMST(Graph g,double const& privacy_parameter){ // The main function taking the graph g and a privacy_parameter
  // that is in fact the parameter at which the exponential distribution have to be parametrized such that the mechanism is an exponential mechanism. It outputs a spanning tree topology computed by PAMST algorithm (the one from the paper).

  unsigned long int N_Vertex= num_vertices(g) ;
  unsigned long int counter =0;
  Rcpp::NumericMatrix S_E (N_Vertex-1, 3);

  std::unordered_set<Vertex> S_V;
  S_V.reserve(N_Vertex);


  Edge e_current;
  Vertex selected_node;

  std::default_random_engine generator{static_cast<long unsigned int>(time(0))};
  std::uniform_real_distribution<double> distribution(0.0,1.0);
  double candidate =-1;
  double current_weight;

  Vertex start = 0;
  S_V.insert(start);


  boost::graph_traits<Graph>::out_edge_iterator e,e_end;
  selected_node=start;

  std::vector<Edge> edges;

  while(counter < (N_Vertex-1) ){ // do it until every node has been labeled but one

    current_weight=std::numeric_limits<float>::infinity(); // redefine current weight to be infinity
    auto S_Vstart = S_V.begin();
    auto S_Vfinish = S_V.end();

    for ( auto current_node = S_Vstart; current_node != S_Vfinish; ++current_node){ // for every node in the set of labeled nodes

      edges.clear();

      boost::tie(e, e_end) = out_edges(*current_node, g);
      for ( ; e != e_end; ++e ) {
        if( S_V.find(target(*e, g)) == S_V.end() ) { // verification that the edge is not between to nodes of S_V
          edges.push_back(*e);
        }
      }

      for (auto e : edges) {
          candidate= g[e].weight - log(1- distribution(generator))/privacy_parameter; // computing the noisy weight of the candidate
          if(candidate <= current_weight){ // if the weight is better than the current better one, change places
            e_current = e;
            current_weight = candidate;
          }
      }


    }

    S_E(counter,0)= (double) source(e_current,g);
    S_E(counter,1)=(double) target(e_current,g); // Forcing the double convertion
    S_E(counter,2)=(double) g[e_current].weight;// add the edge with the min noisy value to the selected edges

    selected_node=source(e_current, g); // select the node forming the selected edge that is not already in S_V and insert it in the set of node labeled
    if( S_V.find(selected_node) == S_V.end() ) {
      S_V.insert(selected_node);
    }else{
      selected_node=target(e_current, g);
      S_V.insert(selected_node);
    }
    ++counter;
  }

  return(S_E);
}


 Rcpp::NumericMatrix approximated_MST_with_Kruskal(Graph const& g,double const& privacy_parameter){ // The main function taking th graph g and a privacy_parameter
   // that is in fact the parameter at which the exponential distribution have to be parametrized such that the mechanism is an exponential mechanism. The usePAMST boolean indicates if the Prim like implementation (the one from the paper) should be used or not. If usePAMST==False, a Kruskal-like implementation is used insted (the one from mitrovic et al. 2017).

   unsigned long int N_Vertex= num_vertices(g) ;
   unsigned long int counter =0;
   Rcpp::NumericMatrix S_E (N_Vertex-1, 3);
   Edge e_current;

   std::default_random_engine generator{static_cast<long unsigned int>(time(0))};
   std::uniform_real_distribution<double> distribution(0.0,1.0);
   double candidate =-1;
   double current_weight;

   auto alledges = boost::edges(g);
   boost::graph_traits<Graph>::out_edge_iterator e,e_end;

   std::vector<VertexIndex> rank(num_vertices(g));
   std::vector<Vertex> parent(num_vertices(g));

   typedef VertexIndex* Rank;
   typedef Vertex* Parent;

   disjoint_sets<Rank, Parent> ds(&rank[0], &parent[0]);
   initialize_incremental_components(g, ds);

   while(counter < (N_Vertex-1)){// do it until every node has been labeled but one

     double current_weight=std::numeric_limits<float>::infinity(); // redefine current weight to be infinity

     for(auto current_edge=alledges.first; current_edge != alledges.second; ++current_edge){
       auto u =ds.find_set(source(*current_edge,g));
       auto v =ds.find_set(target(*current_edge,g));

       if( u!=v || counter==0 ){ // verification that the edge is not between to nodes of S_V
         candidate= g[*current_edge].weight - log(1- distribution(generator))/privacy_parameter;// computing the noisy weight of the candidate

          if(candidate <= current_weight){// if the weight is better than the current better one, change places
            e_current=*current_edge;
            current_weight=candidate;
          }
        }
        if(candidate==-1){
            std::cout<<"Error ?!"<<std::endl;
        }

      }

    S_E(counter,0)= (double) source(e_current,g);
    S_E(counter,1)=(double) target(e_current,g); // Forcing the double convertion
    S_E(counter,2)=(double) g[e_current].weight;// add the edge with the min noisy value to the selected edges

    ds.union_set(source(e_current, g),target(e_current, g));
    ++counter;
  }
   return(S_E);
 }

//' Compute approximate MST from an arbitrary graph using PAMST method
//' @param order Number of nodes in the graph
//' @param Elist edges list as a data.frame
//' @param privacydegree degree of privacy
//' @return The approximate MST using the PAMST algorithm
//' @export
//' @examples
//' n <- 60
//' prob <- 0.1
//' ## Generate random Erdos-Renyi graph
//' graph <- erdos.renyi.game(n, prob, type="gnp",directed = FALSE, loops = FALSE)
//' ## Assign random weights to the edges, using an uniform probability distribution
//' E(graph)$weight <- runif(ecount(graph),0,10)
//' eps <- 0.6
//' Dataframe=igraph::as_data_frame(graph,what="edges")
//' PMST <- PrivateSpanningTree(n,Dataframe,eps,TRUE)
//' print(sum( PMST[,3] ))
//' print(sum(E(mst(graph))$weight))
//'
//' ## plot the resulting MST
//' ## Be careful the ids are not allowed to be 0 !
//' approxMST <- igraph::graph_from_edgelist(cbind(PMST[,1]+1,PMST[,2]+1),directed = FALSE)
//' mylayout <- layout.auto(graph)
//' par(mfrow=c(1,2))
//' plot(graph, layout=mylayout, vertex.size=5, vertex.label=NA)
//' plot(approxMST, layout=mylayout, vertex.size=5, vertex.label=NA)
//'
// [[Rcpp::export]]
Rcpp::NumericMatrix PrivateSpanningTree(int order, Rcpp::DataFrame Elist, double privacydegree, bool usePAMST){
  Rcpp::NumericMatrix Edgelist = DFtoNM( Elist); // transform the edgelist (dataframe type) to a
  //NumericMatrix type edgelist  for C++ to understand it
  int numberofnode = order-1;
  int numberofedge =  Edgelist.nrow();
  Graph g(order);

  for (int i=0;i<numberofedge;++i) {// filling the weights and creating the edges
    Vertex u = Edgelist(i,0)-1; // setting the nodes that the edge descriptor will use
    Vertex v = Edgelist(i,1)-1;
    Edge e; // initializing an edge descriptor
    bool ok;
    boost::tie(e, ok) = boost::add_edge(u,v, g); // boost::add_edge return std::pair<Edge,bool>,tie do it for us
    if(Edgelist(i,2)!=0){ // if the weight is not 0 fill the weight
      if (ok){  // also if the edge creation worked
        g[e].weight = Edgelist(i,2);
      }
    }
  }
  Rcpp::NumericMatrix mstexpmechanism;

  if(usePAMST){ // if UsePAMST is true, return the private spanning tree computed by PAMST algorithm
  mstexpmechanism=approximated_MST_with_PAMST(g,privacydegree/(4*(numberofnode)));
  }
  else{ // otherwise return the private spanning tree computed by Kruskal algorithm
  mstexpmechanism=approximated_MST_with_Kruskal(g,privacydegree/(4*(numberofnode)));
  }
  //call the c++ function to approximate the MST either using prim+exponential mechanism or kruskal+exponential,
  //the privacy parameter here is the parameter that the exponential distribution
  //will use therefore it is epsilon'/2*delta_u (remeber that we use epsilon'=epsilon/(graph oreder -1), and here delta_u=2)
return(mstexpmechanism);
}
