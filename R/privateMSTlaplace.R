simlap <- function(n=1,loc=0,scale=1){
  u <- runif(n)
  return(loc-scale*(log(2*u)*(u<1/2)-log(2*(1-u))*(u>1/2)))
}


#' Compute a private MST using the Laplace mechanism
#'
#' @param graph igraph graph object from which a private Minimum Spanning Tree should be computed
#' @param epsilon privacy parameter of the algorithm
#' @return Approximate MST using Laplace mechanism
#' @references
#' \insertRef{Sealfon_2016}{privateMST}
#' @importFrom igraph mst
#' @importFrom igraph ecount
#' @export
#' @examples
#' n <- 70
#' prob <- 0.1
#' ## Generate random Erdos-Renyi graph
#' graph <- erdos.renyi.game(n, prob, type="gnp",directed = FALSE, loops = FALSE)
#' ## Assign random weights to the edges, using an uniform probability distribution
#' E(graph)$weight <- runif(ecount(graph),0,10)
#' eps <- 0.6
#' approxMSTlaplace <- laplaceMST(graph, epsilon = eps)
#' print(sum(E(approxMSTlaplace)$weight))
#' print(sum(E(mst(graph))$weight))
#'
#' ## plot the resulting MST
#' mylayout <- layout.auto(graph)
#' par(mfrow=c(1,2))
#' plot(graph, layout=mylayout, vertex.size=5, vertex.label=NA)
#' plot(approxMSTlaplace, layout=mylayout, vertex.size=5, vertex.label=NA)
#'
laplaceMST <- function(graph, epsilon) {
  aux <- graph
  E(aux)$weight <- E(aux)$weight + simlap(ecount(aux),scale=1/epsilon)
  return(mst(graph=aux, weights = E(aux)$weight, algorithm = "prim"))
}
