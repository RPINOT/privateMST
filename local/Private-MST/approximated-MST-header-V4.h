
Rcpp::NumericMatrix DFtoNM( Rcpp::DataFrame x) { // function to transform the Data frame that R will provide into a NumericMatrix for c++ to understand it
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

typedef boost::adjacency_list<boost::listS,boost::vecS,boost::undirectedS,Nodeclass,Edgeclass> Graph; // basic construction of a graph using boost::adjacency_list according to the characteristics we gave above
typedef Graph::vertex_descriptor NodeID; // vertex and edges descriptors are basics from  boost::graph that will help us call the node by their name and the edges by a couple of nodes forming it
typedef Graph::edge_descriptor EdgeID;

Rcpp::NumericMatrix approximated_MST(Graph const& g,double const& privacy_parameter ){ // The main function taking th graph g and a privacy_parameter
  // that is in fact the parameter at which the exponential distribution have to be parametrized such that the mechanism is an exponential mechanism

  unsigned long int N_Vertex= num_vertices(g) ;
  unsigned long int counter =0;

  Rcpp::NumericMatrix S_E (N_Vertex-1, 3);

  NodeID start = 0;
  std::unordered_set<NodeID> S_V ;
  S_V.insert(start);

  boost::graph_traits<Graph>::out_edge_iterator e,e_end;
  NodeID selected_node=start;
  EdgeID e_current;

  std::default_random_engine generator;
  std::exponential_distribution<double> rexp(privacy_parameter);

  std::vector<EdgeID> CurrentEdgeVector(N_Vertex);
  std::vector<double> Currentweight(N_Vertex);
  double candidate =-1;

  while(counter < (N_Vertex-1) ){// do it until every node has been labeled but one
    double current_weight=std::numeric_limits<float>::infinity(); // redefine current weight to be infinity
    std::fill(CurrentEdgeVector.begin(), CurrentEdgeVector.end(), e_current);
    std::fill(Currentweight.begin(), Currentweight.end(), current_weight);

  for ( auto current_node = S_V.begin(); current_node != S_V.end(); ++current_node){ // for every node in the set of labeled nodes
  for (boost::tie(e, e_end) = out_edges(*current_node, g); e != e_end; ++e){ // for every edge that is incident to this node

    if( S_V.find(target(*e,g)) == S_V.end()  ){ // verification that the edge is not between to nodes of S_V
      candidate= g[*e].weight - rexp(generator); // computing the noisy weight of the candidate
      if(candidate <= current_weight){// if the weight is better than the current better one, change places
        CurrentEdgeVector[*current_node]=*e;
        Currentweight[*current_node]=candidate;
        current_weight=candidate;
      }
    }

  }
}

for ( auto current_node = S_V.begin(); current_node != S_V.end(); ++current_node){
  if(Currentweight[*current_node]<=candidate){
    e_current=CurrentEdgeVector[*current_node];
  }
}

S_E(counter,0)= (double) source(e_current,g);
S_E(counter,1)=(double) target(e_current,g); // on force la conversion en int pour que R comprenne ce qui se passe
S_E(counter,2)=(double) g[e_current].weight;// add the edge with the min noisy value to the selected edges

selected_node=source(e_current, g); // select the node forming the selected edge that is not already in S_V and insert it in the set of node labeled
if( S_V.find(selected_node) == S_V.end()) {
  S_V.insert(selected_node);
}else{
  selected_node=target(e_current, g);
  S_V.insert(selected_node);
}

++counter;
  }
 //std::cout<<"remplissage du MST ok"<<std::endl;
  return(S_E);
}
