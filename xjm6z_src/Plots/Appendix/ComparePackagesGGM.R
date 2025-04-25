
pacman::p_load(BDgraph, BGGM, tidyverse, mice, MASS, jtools)


pr2pc <- function(K) {
  D.Prec = diag(diag(K)^(-.5))
  R <- diag(2,dim(K)[1])-D.Prec%*%K%*%D.Prec
  colnames(R) <- colnames(K)
  rownames(R) <- rownames(K)
  return(R)
}


#----------------------------------------------------------------------
#----------------------------------------------------------------------
# Function to run single iteration of simulation
#----------------------------------------------------------------------
#----------------------------------------------------------------------

ComparisonFunctionGGM <- function(i, p, n, prob_interaction = 0.5){
  
  #----------------------------------------------------------------------
  # Simulate Data
  #----------------------------------------------------------------------
  #data.sim <- bdgraph.sim( n = ni, p = p, type = "Gaussian", graph = "fixed", K = K_estim)
  data.sim <- bdgraph.sim( n = n, p = p, type = "Gaussian", graph = "random", prob = prob_interaction)
  data <- data.sim$data
  res <- list()
  
  res$true_graph <- data.sim$G[upper.tri(data.sim$G)]
  res$true_sigma <- data.sim$sigma[upper.tri(data.sim$sigma)]
  
  #----------------------------------------------------------------------
  # Estimation
  #----------------------------------------------------------------------
  
  # obtain BGGM fit
  start_BGGM <- Sys.time()
  fit_BGGM <- BGGM::explore(Y = data, type = "continuous", 
                            iter = 5e4, impute = F, seed = NULL)
  select_BGGM <- BGGM::select(fit_BGGM)
  end_BGGM <- Sys.time()
  out_bggm <- summary(select_BGGM, alternative = "exhaustive")
  res$BGGM_parameter <- out_bggm$summary$Post.mean
  res$BGGM_incprob <- out_bggm$summary$Pr.H1
  res$BGGM_structure <- (out_bggm$summary$Pr.H1 > 0.5)*1
  
  
  # obtain BDgraph fit prior 0.5
  start_bdgraph <- Sys.time()
  fit_bdgraph <- BDgraph::bdgraph(data = data, method = "ggm", g.prior = 0.5, iter = 5e4)
  end_bdgraph <- Sys.time()
  pr <- pr2pc(fit_bdgraph$K_hat)
  res$bdgraph_parameter <- pr[upper.tri(pr)]
  incprobs <- plinks(fit_bdgraph)
  res$bdgraph_incprobs <- incprobs[upper.tri(incprobs)]
  structure <- BDgraph::get_graph(fit_bdgraph)
  res$bdgraph_structure <- structure[upper.tri(structure)]
  
  #----------------------------------------------------------------------
  # Comparison
  #----------------------------------------------------------------------
  
  # difference between parameter estimates
  diff_pars_BGGM <- sum((res$BGGM_parameter - res$true_sigma)^2) 
  diff_pars_bdgraph <- sum((res$bdgraph_parameter - res$true_sigma)^2) 
  
  
  # difference in inc_probs for included edges
  diff_probsinc_BGGM <- sum(1- res$BGGM_incprob[res$true_graph == 1])/sum(res$true_graph) 
  diff_probsinc_bdgraph <- sum(1-res$bdgraph_incprob[res$true_graph == 1])/sum(res$true_graph) 
  
  
  # difference in inc_probs for excluded edges
  diff_probsexc_BGGM <- sum(res$BGGM_incprob[res$true_graph == 0])/sum(res$true_graph==0) 
  diff_probsexc_bdgraph <- sum(res$bdgraph_incprob[res$true_graph == 0])/sum(res$true_graph==0) 
  
  # sensitivity 
  sens_BGGM <- sum(res$BGGM_structure*res$true_graph)/sum(res$true_graph)
  sens_bdgraph <- sum(res$bdgraph_structure*res$true_graph)/sum(res$true_graph)
  
  # Specificity
  spec_BGGM <- sum((res$BGGM_structure==0)*(res$true_graph==0))/sum(res$true_graph==0)
  spec_bdgraph<- sum((res$bdgraph_structure==0)*(res$true_graph==0))/sum(res$true_graph==0)
  
  # Correlation matrix
  cor_BGGM <- cor(res$BGGM_parameter, res$true_sigma)
  cor_bdgraph <- cor(res$bdgraph_parameter, res$true_sigma)
  
  # Time for running
  time_BGGM <-  end_BGGM - start_BGGM
  time_bdgraph <- end_bdgraph - start_bdgraph
  
  perf <- list(n = n, p = p, prob_interaction = prob_interaction, 
               diff_pars_BGGM = diff_pars_BGGM, diff_pars_bdgraph = diff_pars_bdgraph,
               diff_probsinc_BGGM = diff_probsinc_BGGM, diff_probsinc_bdgraph = diff_probsinc_bdgraph, 
               diff_probsexc_BGGM = diff_probsexc_BGGM, diff_probsexc_bdgraph = diff_probsexc_bdgraph,
               sens_BGGM = sens_BGGM , sens_bdgraph = sens_bdgraph, 
               spec_BGGM = spec_BGGM , spec_bdgraph = spec_bdgraph,
               cor_BGGM = cor_BGGM, cor_bdgraph = cor_bdgraph, 
               time_BGGM = time_BGGM , time_bdgraph = time_bdgraph 
  )
  return(perf)
  
}
ComparisonFunctionGGM(1, 5, 100)

#----------------------------------------------------------------------
#----------------------------------------------------------------------
# Parallel loop to run through all simulation setups
#----------------------------------------------------------------------
#----------------------------------------------------------------------

rep_iter <- 100
prob_interaction <- .5
p <- c(10, 15)
n <- c(100, 500, 1000, 2500)
numCores <- detectCores()
res_Comparison_GGM <- list()

for (pi in 1:length(p)) {
  for (ni in 1:length(n)){
    
    pit <- p[pi]
    nit <- n[ni]
    res <- mclapply(1:rep_iter, ComparisonFunctionGGM, p = pit, n = nit,
                    mc.cores = numCores)
    res_Comparison_GGM <- append(res_Comparison_GGM, res)
    save(res_Comparison_GGM, file = paste0("Performance_GGM_Simulation_", nit,"_", pit, ".RData"))
    
    print(nit)
  }
  print(pit)
}


#----------------------------------------------------------------------
#----------------------------------------------------------------------
# Extract and visualize results
#----------------------------------------------------------------------
#----------------------------------------------------------------------

res_unlist <- as.numeric(as.vector(unlist(res_Comparison_GGM)))
res <- as.data.frame(matrix(res_unlist, ncol = 17, byrow = T))
colnames(res) <- c("n","p", "prob_interaction" , 
                   "diff_pars_BGGM", "diff_pars_bdgraph" ,
                   "diff_probsinc_BGGM", "diff_probsinc_bdgraph",
                   "diff_probsexc_BGGM", "diff_probsexc_bdgraph", 
                   "sens_BGGM", "sens_bdgraph",
                   "spec_BGGM", "spec_bdgraph", 
                   "cor_BGGM", "cor_bdgraph",
                   "time_BGGM", "time_bdgraph")

# ------------------------------------------------------------------------------
write.csv(res, file = "simulationggm.csv")

perf <- res %>% 
  gather(estimate, value, diff_pars_BGGM:time_bdgraph) %>% 
  separate(estimate, into = c(paste0("name", 1:3)), sep = "_")

perf_diff <- perf %>% 
  filter(name1 == "diff") %>% 
  dplyr::select(-name1)%>% 
  rename(outcome = name2, package = name3)
perf_other <- perf %>% 
  filter(name1 !="diff") %>% 
  dplyr::select(-name3) %>% 
  rename(outcome = name1, package = name2)

perf_ggm <- rbind(perf_diff, perf_other) %>% 
  mutate(n = as.factor(n)) %>% 
  mutate(across('package', str_replace, 'bdgraph', 'BDgraph'), 
         across('outcome', str_replace, 'sens', 'Sensitivity'), 
         across('outcome', str_replace, 'spec', 'Specificity'), 
         across('outcome', str_replace, 'time', 'Time'), 
         across('outcome', str_replace, 'cor', 'Correlation'), 
         across('outcome', str_replace, 'probsexc', 'Error PIP exc.'), 
         across('outcome', str_replace, 'probsinc', 'Error PIP inc.'))

pdf("Plot/Performance_GGM.pdf", width=15, height=10)
perf_ggm %>% 
  filter(outcome != "pars") %>% 
  filter(outcome != "Correlation") %>% 
  filter(outcome != "Time") %>% 
  ggplot(aes(x = n, y = value, col = package)) +
  # geom_violin(alpha = 0.1) +
  # geom_jitter(alpha = 0.3)+
  geom_boxplot(alpha=0.4, size = 1.1) + 
  scale_discrete_manual(values=c("#E69F00", "#175d8d"), aesthetics = "colour", name = "Package") +
  facet_wrap(p~outcome, scales = "free", nrow = 2)+
  theme_minimal()+
  coord_cartesian(ylim = c(0, 1)) +
  theme_apa(remove.x.gridlines = F, remove.y.gridlines = F, legend.font.size = 13.5,
            x.font.size = 16, y.font.size = 16, facet.title.size = 16) +
  theme(legend.position = "right", axis.text=element_text(size=16), 
        legend.background = element_rect(fill = NULL), panel.border = element_blank(), 
        axis.line = element_line(colour = "black", linewidth = 1.1), axis.ticks.length=unit(.25, "cm"),
        axis.ticks = element_line(linewidth= 1), panel.spacing = unit(3, "lines"), 
        strip.text.y = element_text(angle = 0), 
        axis.title=element_text(size=16,face="bold")) +
  ylab("") +
  xlab("Number of Cases")
dev.off()

perf_ising %>% 
  filter(outcome == "Time") %>% 
  dplyr::select(value) %>% 
  max()

