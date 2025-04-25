# !!!!! 

# Note, that this script is no longer up to date as all of the functionalities can be achieved with 
# the R package easybgm. 
# For an example script of easybgm see the annotated script here: https://osf.io/9gfrz
# The most up to date version of easybgm can be found here: https://github.com/KarolineHuth/easybgm
# Author: k.huth@uva.nl
# Update: 14.03.2025

# !!!!!!













# ============================================================================================== # 
#  Tutorial for Bayesian cross-sectional network analysis: BDgraph, BGGM, bgms                   #
# ============================================================================================== # 


# ---------------------------------------------------------------------------------------------- #
# ---------- TABLE OF CONTENTS ----------------------------------------------------------------- #
# ---------------------------------------------------------------------------------------------- #
# ---------- 01. Essential Procedures ---------------------------------------------------------- #
# ---------- 02. Data Loading ------------------------------------------------------------------ # 
# ---------- 03. BDgraph ----------------------------------------------------------------------- #
# ---------- 04. BGGM -------------------------------------------------------------------------- #
# ---------- 05. rbinnet ----------------------------------------------------------------------- #
# ---------- 06. Data Visualization ------------------------------------------------------------ #
# ---------- 07. ... --------------------------------------------------------------------------- #

# For details on the background of the methods used here, please refer to the 
# introduction paper.

# -----------------------------------------------------------------------------------------------------------------

# =========================
#  1. ESSENTIAL PROCEDURES
# =========================

# 1. Set Working Directory
setwd("~/project") 
# Specify your working directory here

# 2. Install/Load the necessary packages
# install.packages("pacman")
devtools::install_github("MaartenMarsman/rbinnet")
pacman::p_load(BDgraph, BGGM, rbinnet, qgraph, tidyverse, igraph, HDInterval) 
source("R/AuxiliaryFunctions.R") # Download and store the file from the OSF page

# -----------------------------------------------------------------------------------------------------------------

# =========================
#  2. LOAD & CLEAN DATASETS
# =========================

# Load here the dataset you want to analyse. 
# As an exemplary dataset you can download the Answers to the Depression Anxiety Stress Scales data here:
# https://openpsychometrics.org/_rawdata/

data_dass <- read_tsv("data/data.csv")
data <- data_dass[, c((0:41*3)+1, 132:141, 158:162, 165:172)]
data <- data %>% 
  mutate(stress = rowSums(.[, c(1, 6, 8, 11, 12, 14, 22, 27, 29, 32, 33, 35, 39)]), 
         anxiety = rowSums(.[, c(2, 4, 7, 9, 15, 19, 20, 23, 25, 28, 30, 36, 40, 41)]), 
         depression = rowSums(.[, c(3, 5, 10, 13, 16, 17, 21, 24, 26, 31, 34, 37, 38, 42)]), 
         zstress = scale(stress), 
         zanxiety = scale(anxiety), 
         zdepression = scale(depression)) %>% 
  dplyr::select(zstress, zanxiety, zdepression, TIPI1:TIPI10) %>%
  rename(depression = zdepression, stress = zstress, anxiety = zanxiety, 
         extraverted = TIPI1, quarrelsome = TIPI2, dependable = TIPI3, anxious = TIPI4, 
         open = TIPI5, reserved = TIPI6, sympathetic = TIPI7, disordganized = TIPI8, 
         calm = TIPI9, uncreative = TIPI10) %>% 
  sample_n(2000) %>% 
  as.matrix()

n <- nrow(data)
p <- ncol(data)
Nsamples <- 5e3 # Choose at least 5e4 for stable results; for demonstrative purposes lower numbers are used
# -----------------------------------------------------------------------------------------------------------------

# =========================
#  3. BDGRAPH
# =========================

# 1. Fit the model
bdgraph_fit <- BDgraph::bdgraph(data=data,               #(M) n*p matrix of responses
                                method="gcgm",           #(M) type of data      
                                algorithm="rjmcmc",      #(O) type of sampling algorithm
                                iter=Nsamples,           #(O) no. iterations sampler
                                save=TRUE,               #(O) Should samples be stored
                                burnin=Nsamples/10,            #(O) no. burnin iterations sampler
                                g.start = "empty",       #(O) starting point of graph
                                not.cont = c(rep(0, 3), rep(1, 10)), #(O) Specifies not continuous variables if method is ggm
                                df.prior = 3,            #(M) prior: degree of freedom for G-Wishart distribution
                                g.prior = 0.5)           #(M) prior: inclusion probability for edges

# 2. Extract results
bdgraph_res <- list()
#Bayesian model-averaged estimates
bdgraph_res$estimates_bma <- pr2pc(bdgraph_fit$K_hat)
diag(bdgraph_res$estimates_bma) <- 0
bdgraph_res$inc_probs <- as.matrix(BDgraph::plinks(bdgraph_fit))
bdgraph_res$inc_probs  <- bdgraph_res$inc_probs + t(bdgraph_res$inc_probs)
bdgraph_res$BF <- bdgraph_res$inc_probs / (1 - bdgraph_res$inc_probs)
bdgraph_res$structure_bma <- 1*(bdgraph_res$inc_probs > 0.5)
bdgraph_res$samples_posterior <- extractposterior(bdgraph_fit, data, method = "gcgm", not.cont = c(rep(0, 3), rep(1, 10)))[[1]]
bdgraph_res$graph_weights <- bdgraph_fit$graph_weights
bdgraph_res$sample_graph <- bdgraph_fit$sample_graphs

#Maximum aposteriori estimates
S <- t(data) %*% data
bdgraph_res$structure_map <- BDgraph::pgraph(bdgraph_fit, number.g=1)$selected_g[[1]]
bdgraph_res$structure_map <- bdgraph_res$structure_map + t(bdgraph_res$structure_map)
posterior_map_samples <- gwish_samples(bdgraph_res$structure_map, S, nsamples=Nsamples)
bdgraph_res$estimates_map <- vector2matrix(colMeans(posterior_map_samples), p, bycolumn = T)

# Centrality indices
bdgraph_res$centrality_strength <- centrality_strength(bdgraph_res)
bdgraph_res$centrality <- centrality(bdgraph_res)
# -----------------------------------------------------------------------------------------------------------------

# =========================
#  4. BGGM
# =========================

# 1. Fit the model
bggm_fit <- BGGM::explore(data,                        #(M) n*p matrix of responses
                          type = "continuous",         #(O) type of data
                          mixed_type = NULL,           #(O) which data should be treated as ranks
                          prior_sd = 0.25,             #(M) prior distribution standard deviation
                          iter = Nsamples,             #(O) no. iterations sampler
                          progress = FALSE,            #(O) Should a progress bar be plotted?
                          impute = FALSE,              #(O) Should missings be imputed?
                          seed = 1)                    #(O) Integer for random seed

# 2. Extract results
bggm_out <- summary(BGGM::select(bggm_fit, BF_cut = 3))
bggm_res <- list()
bggm_res$estimates_bma <- vector2matrix(bggm_out$summary$Post.mean, p, bycolumn = T)
bggm_res$inc_probs <- vector2matrix(bggm_out$summary$Pr.H1, p, bycolumn = T)
bggm_res$BF <- bggm_res$inc_probs/(1-bggm_res$inc_probs)
bggm_res$structure_bma <- 1*(bggm_res$inc_probs > 0.5)
bggm_res$samples_posterior <- list2matrix(bggm_fit$post_samp$pcors, p)
bggm_res$centrality_strength <- centrality_strength(bggm_res)

# -----------------------------------------------------------------------------------------------------------------

# =========================
#  5. RBINNET
# =========================
# Please note, this package only allows for binary data estimation
# Will not work with the example data provided

# 1. Fit the model
rbinnet_sigma.map <- fit_pseudoposterior(x = as.matrix(data), prior_var = 1)$sigma
rbinnet_fit <- select_structure(x = as.matrix(data),                  #(M) n * p matrix of binary responses
                                sigma = sigma.map,              #(O) p * p matrix of Ising parameters
                                theta = 0.5,                    #(O) prior inclusion probability
                                prior_var_intercepts = 1,       #(O) prior threshold parameter
                                output_samples = T,             #(O) if TRUE, outputs posterior draws
                                number_iterations = 1e5,        #(O) no. iterations Gibbs sampler
                                number_burnin_iterations = 0,   #(O) no. burnin iterations Gibbs sampler
                                precision = .975)               #(O) prior precision for edge selection

# 2. Extract results
rbinnet_res <- list()
rbinnet_res$estimates_bma <- vector2matrix(colMeans(rbinnet_fit$sigma_samples), p, diag = T)
diag(rbinnet_res$estimates_bma) <- 0
rbinnet_res$inc_probs <- vector2matrix(rbinnet_fit$posterior_probability %*% rbinnet_fit$structures, p) 
rbinnet_res$BF <- rbinnet_res$inc_probs/(1 - rbinnet_res$inc_probs)
rbinnet_res$structure_bma <- 1*(rbinnet_res$inc_probs > 0.5)
rbinnet_res$samples_posterior <- rbinnet_fit$sigma_samples

# Centrality indices
rbinnet_res$centrality_strength <- centrality_strength(rbinnet_res)
rbinnet_res$centrality <- centrality(rbinnet_res)
# -----------------------------------------------------------------------------------------------------------------

# =========================
#  6. DATA VISUALISATION
# =========================
# Note: Exchange [Package]_res with the respective package that you used
res <- bdgraph_res

# 0. Plot posterior structure estimates (NOTE: not usable with BGGM)

# 0.a Posterior structure probability
sorted_structure_prob <- as.data.frame(sort(res$graph_weights/sum(res$graph_weights), decreasing = T))
colnames(sorted_structure_prob) <- "posterior_prob"
ggplot(sorted_structure_prob, aes(x = 1:nrow(sorted_structure_prob), y = posterior_prob)) +
  geom_point()+
  ylab("Posterior Structure Probability") +
  xlab("Structure")  +
  theme_minimal()+
  theme(legend.position = c(.85, 0.25), axis.text=element_text(size=20), 
        legend.background = element_rect(fill = NULL), panel.border = element_blank(), 
        axis.line = element_line(colour = "black", size = 1.1), axis.ticks.length=unit(.2, "cm"),
        axis.ticks = element_line(size= .8), legend.text = element_text( size=14), 
        axis.title.x = element_text(size=18,face="bold"), 
        axis.title.y = element_text(size=18,face="bold"),
        text=element_text(  family="Times New Roman"), 
        panel.grid.major = element_blank()
        #, axis.text.y = element_blank(), axis.ticks.y = element_blank()
  )

# 0.b. Posterior complexity probability
complexity <- c()
for(i in 1:length(res$sample_graph)){
  complexity[i] <- sum(as.numeric(unlist(strsplit(res$sample_graph[i], ""))))
}

data_complexity <- tibble(complexity, weights = res$graph_weights)  %>%
  group_by(complexity) %>% 
  summarise(complexity_weight = sum(weights)) %>% 
  mutate(complexity_weight = complexity_weight/sum(complexity_weight))

ggplot(data_complexity, aes(x = complexity, y = complexity_weight))+
  geom_point() + 
  ylab("Posterior Complexity Probability") +
  xlab("Complexity")  +
  theme_minimal()+
  theme(legend.position = c(.85, 0.25), axis.text=element_text(size=20), 
        legend.background = element_rect(fill = NULL), panel.border = element_blank(), 
        axis.line = element_line(colour = "black", size = 1.1), axis.ticks.length=unit(.2, "cm"),
        axis.ticks = element_line(size= .8), legend.text = element_text( size=14), 
        axis.title.x = element_text(size=18,face="bold"), 
        axis.title.y = element_text(size=18,face="bold"),
        text=element_text(  family="Times New Roman"), 
        panel.grid.major = element_blank()
  )


# 1. Plot structure
qgraph(res$structure_bma, layoutScale = c(.8,1), palette = "R",
       theme = "TeamFortress", vsize = 6, edge.width = .3, layout = "spring")

# 2. Plot network model (also referred to as median probability model)
qgraph(res$estimates_bma*res$structure_bma, layout = "spring", 
       layoutScale = c(.8,1), palette = "R",
       theme = "TeamFortress", vsize = 6)

# 3. Plot evidence plot
graph_color <- res$BF
graph_color <-  ifelse(res$BF < 10 & res$BF > .1, graph_color <- "#bfbfbf", graph_color <- "#36648b")
graph_color[res$BF < .1] <- "#990000"

qgraph(matrix(1, ncol = p, nrow = p), edge.color = graph_color, edge.width = 2)

# 4. Plot density of posterior parameter estimates
i <- 10 # Choose respective density that should be ploted
plot(density(res$samples_posterior[, i]), bty = "n", 
     main = paste0("Density parameter ", i))

# 5. Plot strength centrality estimate + HDI 
centrality_means <- res$centrality_strength$centrality_strength_mean
centrality_hdi_intervals <- apply(res$centrality_strength$centrality_strength_samples, MARGIN = 2, FUN = hdi, allowSplit = F)

centrality <- as.data.frame(t(rbind(centrality_means, centrality_hdi_intervals)))

centrality %>%
  mutate(centrality_means = as.numeric(centrality_means)) %>%
  ggplot(aes(x = 1:nrow(centrality), y = centrality_means)) +
  geom_line(size = 2) + 
  geom_point(size = 5) +
  geom_errorbar(aes(ymin = as.numeric(lower), ymax = as.numeric(upper)), size = 1, width = .2) +
  theme_minimal() +
  labs(y = "Centrality", x = "Node") +
  theme(text = element_text(size = 22, family = "Times New Roman"), axis.title.x = element_text(size=22), 
        axis.title.y = element_text(size=22), panel.border = element_blank(),
        axis.line = element_line(colour = "black", size = 1.1), axis.ticks.length=unit(.3, "cm"), 
        legend.position = c(.85, 0.9), legend.text = element_text(size = 20)) +
  scale_x_continuous(breaks = 1:nrow(centrality))


# 6. Plot all centrality measures + HDI 
centrality_means <- apply(res$centrality[, 3:(Nsamples+2)], MARGIN = 1, mean)
centrality_hdi_intervals <- apply(res$centrality[, 3:(Nsamples+2)], MARGIN = 1, FUN = hdi, allowSplit = F)

centrality_summary <- cbind(res$centrality[, 1:2], centrality_means, t(centrality_hdi_intervals))
centrality_summary %>% 
  ggplot(aes(x = node, y=centrality_means, group = measure))+
  geom_line()+
  geom_point()+
  geom_errorbar(aes(y= centrality_means, ymin =lower, ymax = upper), size = .5, width = 0.4)+
  facet_wrap(~ measure, ncol = 4) +
  coord_flip() +
  ylab("") +
  xlab("")

