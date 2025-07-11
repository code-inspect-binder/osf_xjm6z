# !!!!! 

# Note, that this script is no longer up to date as all of the auxiliary functions are now within
# the R package easybgm. 
# For an example script of easybgm see the annotated script here: https://osf.io/9gfrz
# The most up to date version of easybgm can be found here: https://github.com/KarolineHuth/easybgm
# Author: k.huth@uva.nl
# Update: 14.03.2025

# !!!!!!










# 1. Turns vector into matrix
vector2matrix <- function(vec, p, diag = F, bycolumn = F) {
  m <- matrix(0, p, p)
  
  if(bycolumn == F){
    m[lower.tri(m, diag = diag)] <- vec
    m <- t(m)
    m[lower.tri(m)] <- t(m)[lower.tri(m)]
  } else {
    m[upper.tri(m, diag = diag)] <- vec
    m <- t(m)
    m[upper.tri(m)] <- t(m)[upper.tri(m)]
  }
  return(m)
}

# 2. Transform precision into partial correlations for interpretation
pr2pc <- function(K) {
  D.Prec = diag(diag(K)^(-.5))
  R <- diag(2,dim(K)[1])-D.Prec%*%K%*%D.Prec
  colnames(R) <- colnames(K)
  rownames(R) <- rownames(K)
  return(R)
}

# 3. BDgraph stores graphs as byte strings for efficiency
string2graph <- function(Gchar, p) {
  Gvec = rep(0, p*(p-1)/2)
  edges <- which(unlist(strsplit(as.character(Gchar), "")) == 1)
  Gvec[edges] = 1
  G <- matrix(0, p, p)
  G[upper.tri(G)] <- Gvec
  G = G + t(G)
  return(G) 
}

# 4. BDgraph extract posterior distribution for estimates
extractposterior <- function(fit, data, method = c("ggm", "gcgm"), not.cont){
  m <- length(fit$all_graphs)
  k <- 30000 
  n <- nrow(data)
  p <- ncol(data)
  j <- 1
  densities <- rep(0, k)
  #Rs = array(0, dim=c(k, p, p))
  Rs = matrix(0, nrow = k, ncol = (p*(p-1))/2)
  if(method == "gcgm") {
    S <- get_S_n_p(data, method = method, n = n, not.cont = not.cont)$S
  } else {
    S <- t(data) %*% data
  }
  for (i in seq(1, m, length.out=k)) {
    graph_ix <- fit$all_graphs[i]
    G <- string2graph(fit$sample_graphs[graph_ix], p)
    K <- BDgraph::rgwish(n=1, adj=G, b=3+n, D=diag(p) + S)
    
    #Rs[j,,] <- pr2pc(K)
    Rs[j,] <- as.vector(pr2pc(K)[upper.tri(pr2pc(K))])
    densities[j] <- sum(sum(G)) / (p*(p-1))
    j <- j + 1
  }
  return(list(Rs, densities))
}

# 5. Samples from the G-wishart distribution
gwish_samples <- function(G, S, nsamples=1000) {
  p <- nrow(S)
  #Rs <- array(0, dim=c(nsamples, p, p))
  Rs = matrix(0, nrow = nsamples, ncol = (p*(p-1))/2)
  
  for (i in 1:nsamples) {
    K <- BDgraph::rgwish(n=1, adj=G, b=3+n, D=diag(p) + S)*(G + diag(p))
    Rs[i,] <- as.vector(pr2pc(K)[upper.tri(pr2pc(K))])
    #Rs[i,,] <- .pr2pc(K)
  }
  return(Rs)
}


# 6. Centrality of weighted graphs

# Strength centrality only ## FASTER CODE
centrality_strength <- function(res){
  Nsamples <- nrow(res$samples_posterior)
  p <- nrow(res$estimates_bma)
  strength_samples <- matrix(0, nrow = Nsamples, ncol = p)
  for(i in 1:Nsamples){
    strength_samples[i, ] <- rowSums(abs(vector2matrix(res$samples_posterior[i,], p, bycolumn = T)))
  }
  strength_mean <- colMeans(strength_samples)
  strength_median <- apply(strength_samples,2,median)
  return(list(centrality_strength_samples = strength_samples, centrality_strength_mean = strength_mean, 
              centrality_strength_median = strength_median))
}

# Strength, betweenness and closeness centrality ## SLOWER CODE
centrality <- function(res, include = c("Strength", "Closeness", "Betweenness", "ExpectedInfluence")){
  Nsamples <- nrow(res$samples_posterior)
  p <- as.numeric(nrow(res$estimates_bma))
  #degree_samples <- betweenness_samples <- closeness_samples <-influence_samples <- matrix(0, nrow = Nsamples, ncol = p)
  for(i in 1:Nsamples){
    graph <- centralityPlot(vector2matrix(res$samples_posterior[i,], p, bycolumn = T), 
                            include = c("Strength", "Closeness", "Betweenness", "ExpectedInfluence"),
                            verbose = F, print = F)
    if(i > 1){
      centrality_output[, i+2] <- graph$data[, 5]
    } else {
      centrality_output <- graph$data[, 3:5]
    }
  }
  return(centrality_output)
}

# 7. Centrality of unweighted graphs 
centrality_graph <- function(fit, include = c("degree", "closeness", "betweenness") ){
  # amount of visited structures
  len <- length(fit$sample_graphs)
  
  # objects to store graph centrality measures
  degree <- matrix(0, nrow = len, ncol = p)
  betweenness <- matrix(0, nrow = len, ncol = p)
  closeness <- matrix(0, nrow = len, ncol = p)
  
  # Obtain centrality measures for each graph
  for (i in 1:len){
    graph_matrix <- vector2matrix(as.numeric(unlist(strsplit(fit$sample_graphs[1], ""))), p , bycolumn = T)
    graph_graph <- igraph::as.undirected(graph.adjacency(graph_matrix, weighted = T))
    
    degree[i, ] <- igraph::degree(graph_graph)
    betweenness[i, ] <- igraph::betweenness(graph_graph)
    closeness[i, ] <- igraph::closeness(graph_graph)
  }
  # save centrality measures of interest
  centrality_graph <- list()
  if("degree" %in% include){
    degree_samples <- degree[rep(1:nrow(degree), fit$graph_weights),]
    centrality_graph[["degree_mean"]] <- colMeans(degree_samples)
    centrality_graph[["degree_samples"]] <- degree_samples
  }  
  if("betweenness" %in% include) {
    betweenness_samples <- betweenness[rep(1:nrow(betweenness), fit$graph_weights),]
    centrality_graph[["betweenness_mean"]] <- colMeans(betweenness_samples)
    centrality_graph[["betweenness_samples"]] <- betweenness_samples
  }  
  if  ("closeness" %in% include){
    closeness_samples <- closeness[rep(1:nrow(closeness), fit$graph_weights),]
    centrality_graph[["closeness_mean"]] <- colMeans(closeness_samples)
    centrality_graph[["closeness_samples"]] <- closeness_samples
  }
  return(centrality_graph)
}

# 8. turn list into matrix 
list2matrix <- function(obj, p) {
  nlist <- length(obj)/(p*p)
  m <- obj[, , 1]
  nest <- sum(lower.tri(m))
  res <- matrix(0, nrow = nlist, ncol = nest)
  for(i in 1:nlist){
    m <- obj[, , i]
    res[i, ] <- as.vector(m[lower.tri(m)])
  }
  return(res)
}
