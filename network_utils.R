construct_network <- function(edges, type="laplacian")
{
  nodes <- sort(unique(c(edges[,1], edges[,2])))
  size <- length(nodes)
  #Construct interaction matrix
  indices_s <- cbind(match(edges[,1], nodes), match(edges[,2], nodes))
  indices_g <- cbind(match(edges[,2], nodes), match(edges[,1], nodes))
  
  indices <- rbind(indices_s, indices_g)
  mat <- Matrix(0, nrow = size, ncol= size, sparse = TRUE)
  mat[indices] <- 1
  diag(mat) <- 0
  
  if (type == "laplacian")
  { 
    g <- graph.adjacency(mat, mode="undirected")
    L <- graph.laplacian(g, normalized=TRUE)
    g <- abs(L)
  }
  #D <- matrix(data = rep(0, size), ncol=size, nrow=size)
  #diag(D) <- degree(g)
  #A <- 1/sqrt(degree(g)) * (D - Ad) * 1/sqrt(degree(g))
  g
}

