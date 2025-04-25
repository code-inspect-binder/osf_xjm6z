
pacman::p_load(BDgraph, BGGM, tidyverse, IsingSampler, MASS, parallel, bgms)

#----------------------------------------------------------------------
#----------------------------------------------------------------------
# Function to run single iteration of simulation
#----------------------------------------------------------------------
#----------------------------------------------------------------------

ComparisonFunction <- function(i, p, n, prob_interaction = 0.5){
  
  #----------------------------------------------------------------------
  # Simulate Data
  #----------------------------------------------------------------------
  
  # Ising parameters:
  Graph <- matrix(data = 0, nrow = p, ncol = p)
  Graph[lower.tri(Graph)] <- rbinom(n = p * (p - 1) / 2,
                                    size = 1,
                                    prob = prob_interaction)
  Sigma <- Graph * runif(n = p ^ 2, min = 0.3, max = 1.5)
  true_graph <- Graph + t(Graph)
  true_sigma <- Sigma + t(Sigma)
  Thresholds <- -rowSums(true_sigma) / 2
  # Simulate:
  Data <- IsingSampler(n = n, graph = true_sigma, thresholds = Thresholds)
  
  res <- list()
  res$true_graph <- true_graph[upper.tri(true_graph)]
  res$true_sigma <- true_sigma[upper.tri(true_sigma)]
  
  #----------------------------------------------------------------------
  # Estimation
  #----------------------------------------------------------------------
  
  # bgms
  start_bgms <- Sys.time()
  fit_bgm <- bgm(as.matrix(Data))
  end_bgms <- Sys.time()
  res$bgms_parameter <- fit_bgm$interactions[upper.tri(fit_bgm$interactions)]
  res$bgms_incprob <- fit_bgm$gamma[upper.tri(fit_bgm$gamma)]
  res$bgms_structure <- (res$bgms_incprob > 0.5)*1
  
  # BGGM
  start_BGGM <- Sys.time()
  fit_BGGM <- explore(Data,  type = "binary", iter = 1e4, impute = F, seed = NULL)
  out_bggm <- summary(BGGM::select(fit_BGGM))
  end_BGGM <- Sys.time()
  res$BGGM_parameter <- out_bggm$summary$Post.mean
  res$BGGM_incprob <- out_bggm$summary$Pr.H1
  res$BGGM_structure <- (out_bggm$summary$Pr.H1 > 0.5)*1
  
  # BD graph
  start_bdgraph <- Sys.time()
  fit_bdgraph <- bdgraph.mpl(Data, method = "dgm-binary", iter = 1e4, save = T)
  end_bdgraph <- Sys.time()
  incprobs <- plinks(fit_bdgraph)
  res$bdgraph_incprobs <- incprobs[upper.tri(incprobs)]
  structure <- BDgraph::get_graph(fit_bdgraph)
  res$bdgraph_structure <- structure[upper.tri(structure)]
  
  #----------------------------------------------------------------------
  # Comparison
  #----------------------------------------------------------------------
  
  # difference between parameter estimates
  diff_pars_bgms <- sum((res$bgms_parameter - res$true_sigma)^2) 
  diff_pars_BGGM <- sum((res$BGGM_parameter - res$true_sigma)^2) 
  
  
  # difference in inc_probs for included edges
  diff_probsinc_bgms <- sum(1- res$bgms_incprob[res$true_graph == 1])/sum(res$true_graph) 
  diff_probsinc_BGGM <- sum(1- res$BGGM_incprob[res$true_graph == 1])/sum(res$true_graph) 
  diff_probsinc_bdgraph <- sum(1-res$bdgraph_incprob[res$true_graph == 1])/sum(res$true_graph) 
  
  
  # difference in inc_probs for excluded edges
  diff_probsexc_bgms <- sum(res$bgms_incprob[res$true_graph == 0])/sum(res$true_graph==0) 
  diff_probsexc_BGGM <- sum(res$BGGM_incprob[res$true_graph == 0])/sum(res$true_graph==0) 
  diff_probsexc_bdgraph <- sum(res$bdgraph_incprob[res$true_graph == 0])/sum(res$true_graph==0) 
  
  # sensitivity 
  sens_bgms <- sum(res$bgms_structure*res$true_graph)/sum(res$true_graph)
  sens_BGGM <- sum(res$BGGM_structure*res$true_graph)/sum(res$true_graph)
  sens_bdgraph <- sum(res$bdgraph_structure*res$true_graph)/sum(res$true_graph)
  
  
  # Specificity
  spec_bgms <- sum((res$bgms_structure==0)*(res$true_graph==0))/sum(res$true_graph==0)
  spec_BGGM <- sum((res$BGGM_structure==0)*(res$true_graph==0))/sum(res$true_graph==0)
  spec_bdgraph <- sum((res$bdgraph_structure==0)*(res$true_graph==0))/sum(res$true_graph==0)
  
  
  # Correlation matrix
  cor_bgms <- cor(res$bgms_parameter, res$true_sigma)
  cor_BGGM <- cor(res$BGGM_parameter, res$true_sigma)
  
  # Time for running
  time_bgms <-  end_bgms - start_bgms
  time_BGGM <-  end_BGGM - start_BGGM
  time_bdgraph <- end_bdgraph - start_bdgraph
  
  perf <- list(n = n, p = p, prob_interaction = prob_interaction, 
              diff_pars_bgms = diff_pars_bgms, diff_pars_BGGM = diff_pars_BGGM, 
              diff_probsinc_bgms = diff_probsinc_bgms, diff_probsinc_BGGM = diff_probsinc_BGGM, diff_probsinc_bdgraph = diff_probsinc_bdgraph, 
              diff_probsexc_bgms = diff_probsexc_bgms, diff_probsexc_BGGM = diff_probsexc_BGGM, diff_probsexc_bdgraph = diff_probsexc_bdgraph,
              sens_bgms = sens_bgms, sens_BGGM = sens_BGGM , sens_bdgraph = sens_bdgraph, 
              spec_bgms = spec_bgms, spec_BGGM = spec_BGGM , spec_bdgraph = spec_bdgraph,
              cor_bgms = cor_bgms, cor_BGGM = cor_BGGM, 
              time_bgms = time_bgms, time_BGGM = time_BGGM , time_bdgraph = time_bdgraph 
              )
  return(perf)
}

# ComparisonFunction(1, 5, 100)

#----------------------------------------------------------------------
#----------------------------------------------------------------------
# Parallel loop to run through all simulation setups
#----------------------------------------------------------------------
#----------------------------------------------------------------------

rep_iter <- 100
prob_interaction <- .5
p <- c(15) # c(10, 15)
n <- c(100, 500, 1000, 2500)
numCores <- detectCores()
#res_Comparison <- list()
for (pi in 1:length(p)) {
  for (ni in 1:length(n)){

    pit <- p[pi]
    nit <- n[ni]
    res <- mclapply(1:rep_iter, ComparisonFunction, p = pit, n = nit,
                               mc.cores = numCores)
    res_Comparison <- append(res_Comparison, res)
    save(res_Comparison, file = paste0("Performance_Ising_Simulation_", nit,"_", pit, ".RData"))
    
    print(nit)
  }
  print(pit)
}

#----------------------------------------------------------------------
#----------------------------------------------------------------------
# Extract and visualize results
#----------------------------------------------------------------------
#----------------------------------------------------------------------

res_unlist <- as.numeric(as.vector(unlist(res_Comparison)))
res <- as.data.frame(matrix(na.omit(res_unlist), ncol = 22, byrow = T))
colnames(res) <- c("n","p", "prob_interaction" , 
                   "diff_pars_bgms" , "diff_pars_BGGM", 
                   "diff_probsinc_bgms", "diff_probsinc_BGGM", "diff_probsinc_bdgraph",
                   "diff_probsexc_bgms", "diff_probsexc_BGGM", "diff_probsexc_bdgraph", 
                   "sens_bgms", "sens_BGGM", "sens_bdgraph",
                   "spec_bgms", "spec_BGGM", "spec_bdgraph", 
                   "cor_bgms", "cor_BGGM", 
                   "time_bgms", "time_BGGM", "time_bdgraph")

# ------------------------------------------------------------------------------
write.csv(res, file = "simulationising.csv")
res <- read.csv("simulationising.csv")
perf <- res %>% 
  gather(estimate, value, diff_pars_bgms:time_bdgraph) %>% 
  separate(estimate, into = c(paste0("name", 1:3)), sep = "_")

perf_diff_ising <- perf %>% 
  filter(name1 == "diff") %>% 
  dplyr::select(-name1)%>% 
  rename(outcome = name2, package = name3)
perf_other_ising <- perf %>% 
  filter(name1 !="diff") %>% 
  dplyr::select(-name3) %>% 
  rename(outcome = name1, package = name2)

perf_ising <- rbind(perf_diff_ising, perf_other_ising) %>% 
  mutate(n = as.factor(n)) %>% 
  mutate(across('package', str_replace, 'bdgraph', 'BDgraph'), 
         across('outcome', str_replace, 'sens', 'Sensitivity'), 
         across('outcome', str_replace, 'spec', 'Specificity'), 
         across('outcome', str_replace, 'time', 'Time'), 
         across('outcome', str_replace, 'cor', 'Correlation'), 
         across('outcome', str_replace, 'probsexc', 'Error PIP exc.'), 
         across('outcome', str_replace, 'probsinc', 'Error PIP inc.'))

pdf("Plot/Performance_Ising.pdf", width=15, height=10)
perf_ising %>% 
  filter(outcome != "pars") %>% 
  filter(outcome != "Correlation") %>% 
  filter(outcome != "Time") %>% 
  ggplot(aes(x = n, y = value, col = package)) +
  # geom_violin(alpha = 0.1) +
  # geom_jitter(alpha = 0.3)+
  geom_boxplot(alpha=0.4, size = 1.1) +
  scale_discrete_manual(values=c( "#E69F00", "#175d8d","#a64646"), aesthetics = "colour", name = "Package") +
  facet_wrap(p~outcome, scales = "free", nrow = 2)+
  theme_minimal() +
  coord_cartesian(ylim = c(0, 1))  +
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

