library(igraph)

network_smoothing <- function(net, mat_intensities, y, iter=1, alpha, network_type="laplacian", stop_val = 0)
{
    if (dim(mat_intensities)[1] == dim(net)[1])
        {
            if (stop_val != 0)
                {
                    delta = 0
                    i = 0
                    while(delta < stop_val)
                        {
                            t_mat_intensities <- mat_intensities
                            mat_intensities <- apply(t_mat_intensities, 2, function(x) as.vector(propagate_algebra(net, x, alpha, 1)))
                            delta <- mean((mat_intensities - t_mat_intensities)^2, na.rm = TRUE)
                            cat(paste(delta, 'delta \n'))
                            cat(paste(i, 'i \n'))
                            i = i + 1
                        }
                    smooth_net <- apply(mat_intensities, 2 , I)
                    smooth_net <- as.matrix(smooth_net)
                    colnames(smooth_net) <- y
                    smooth_net
                }
            else
                {
                    smooth_net <- apply(mat_intensities, 2, function(x) as.vector(propagate_algebra(net, x, alpha, iter)))
                    smooth_net <- apply(smooth_net, 2 , I)
                    smooth_net <- as.matrix(smooth_net)
                    colnames(smooth_net) <- y
                    smooth_net
                }
        }
    else
        {
            cat("Please ensure correspondance between network and intensity matrix.")
        }
    
}



transformation <- function(data)
{
  data <- data[c(2:dim(data)[2])]
  data[is.na(data)] <- 0.0000000000000000001
  data <- data.matrix(data, rownames.force = NA)
  data <- normalize.quantiles(data) 
  data <- log2(data)
  return(data)
}

randomize_repeat_numer <- function(x)
{
  for (i in 1:nrow(x)) 
  {
    for (j in 1:ncol(x))
    { 
      if (x[i,j]==0) { # assume i have too many 1 in my_column
        x[i,j] = as.numeric(runif(1,min=0.001,max=0.1)) # replace 1 with nearby values randomly
      }
    }
  }
  x
}

construct_network <- function(edges, type="laplacian")
{
  nodes <- unique(sort(c(edges[,1], edges[,2])))
  size <- length(nodes)
  cat(paste("Network with:", size, "nodes\n"))
  #Construct interaction matrix
  
  indices_s <- cbind(match(edges[,1], nodes), match(edges[,2], nodes))
  indices_g <- cbind(match(edges[,2], nodes), match(edges[,1], nodes))
  
  indices <- rbind(indices_s, indices_g)
  
  mat <- Matrix(0, nrow = size, ncol= size, sparse = TRUE)
  mat[indices] <- 1
  diag(mat) <- 0
  
  
  if (type == "laplacian")
      {
          g <- graph.adjacency(as.matrix(mat), mode="undirected")
          L <- graph.laplacian(g, normalized=TRUE)
          G <- abs(L)
          return(G)
      }
  if (type == "adjacency")
      {
          g <- graph.adjacency(as.matrix(mat), mode="undirected")
          #L <- graph.laplacian(g, normalized=TRUE)
          G <- abs(g)
          return(G)
      }
   
  #D <- matrix(data = rep(0, size), ncol=size, nrow=size)
  #diag(D) <- degree(g)
  #A <- 1/sqrt(degree(g)) * (D - Ad) * 1/sqrt(degree(g))
  
}

extract_exp <- function(x, val)
    {
        values_1 <- unlist(lapply(1:val, function(x) {paste("value_1.", x, sep="")}))
        values_2 <- unlist(lapply(1:val, function(x) {paste("value_2.", x, sep="")}))
        sample_1 <- unlist(lapply(1:val, function(x) {paste("sample_1.", x, sep="")}))
        sample_2 <- unlist(lapply(1:val, function(x) {paste("sample_2.", x, sep="")}))
        exp <- c("value_1", values_1, "value_2", values_2)
        sample <- c("sample_1", sample_1, "sample_2", sample_2)
        exp <- x[exp]
        sample <- x[sample][1,]
        colnames(exp) <- sample
        exp
    }

network_mapping <- function(network, exp_list, type = "laplacian", merge.by = 'gene_id', global = TRUE)
    {
        #Contruct a matrix representation of a network from a list of interactions.
        interaction_list <- cbind(as.vector(network[,1]), as.vector(network[,2]))
        ppi <- sort(unique(c(as.vector(network[,1]), as.vector(network[,2]))))
        ppi <- as.data.frame(ppi)
        colnames(ppi) <- merge.by      
        
        expr_start <- merge(ppi, exp_list, by=merge.by, all.x = global)
        expr_start <- expr_start[!duplicated(expr_start[,merge.by]),]       
                                        #interact <- interaction_list
        select_prot <- as.matrix(sort(expr_start[,c(merge.by)]))
        all_prot <- as.matrix(ppi)
        
                                        #select_int <- ifelse((interact[,1] %in% select_prot) | (interact[,2] %in% select_prot), TRUE, FALSE)
        g <- construct_network(interaction_list, type)
        S <- g[all_prot %in% select_prot, all_prot %in% select_prot]
        cat(paste("\n", "network defined by", dim(S), "matrix"))
        return(list('G' = S, 'mat_intensities' = expr_start))
    }
