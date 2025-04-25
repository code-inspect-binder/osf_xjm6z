# !!!!! 

# Note, that this script is no longer up to date and some of the functionalities will not work anymore with the current version of easybgm. 
# We keep this script on OSF as it belongs to the original publication. 
# For an example script of easybgm see the annotated script here: https://osf.io/9gfrz
# The most up to date version of easybgm can be found here: https://github.com/KarolineHuth/easybgm
# Author: k.huth@uva.nl
# Update: 14.03.2025

# !!!!!!









# ============================================================================================== # 
#  Tutorial for Bayesian cross-sectional network analysis: BDgraph, BGGM, rbinnet                #
# ============================================================================================== # 


# ---------------------------------------------------------------------------------------------- #
# ---------- TABLE OF CONTENTS ----------------------------------------------------------------- #
# ---------------------------------------------------------------------------------------------- #
# ---------- 01. Essential Procedures ---------------------------------------------------------- #
# ---------- 02. Data Loading ------------------------------------------------------------------ # 
# ---------- 03. BDgraph ----------------------------------------------------------------------- #
# ---------- 04. Data Visualization ------------------------------------------------------------ #
# ---------- 05. ... --------------------------------------------------------------------------- #


# -----------------------------------------------------------------------------------------------------------------

# =========================
#  1. ESSENTIAL PROCEDURES
# =========================


# 1. Set Working Directory
setwd("~/Library/CloudStorage/OneDrive-UvA/01_Projects/M003_BayesianNetworkIntro") # Specify your working directory here

# 2. Install/Load the necessary packages
library(easybgm)
library(bgms)
library(BDgraph)
library(BGGM)
library(qgraph)
library(igraph)
library(HDInterval)
library(tidyverse)
library(jtools)

# -----------------------------------------------------------------------------------------------------------------

# =========================
#  2. LOAD & CLEAN DATASETS
# =========================

# For an exemplary dataset we used the Answers to the Depression Anxiety Stress Scales from:
# https://openpsychometrics.org/_rawdata/
set.seed(123)
data_dass <- read_tsv("data/data.csv")
data <- data_dass[, c((0:41*3)+1, 132:141, 158:162, 165:172)] %>% sample_n(3000)
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
  as.matrix()
# write_csv(as.data.frame(data), "data/data_cleaned_new.csv")

data <- as.matrix(read.csv("data/data_cleaned_new.csv"))
colnames(data) <- 1:13
n <- nrow(data)
p <- ncol(data)
Nsamples <- 5e4

# -----------------------------------------------------------------------------------------------------------------


# ==========================
#  3. Model fit with BDGRAPH
# ==========================

# 1. Fit the model
bdgraph_fit <- easybgm::easybgm(data = data, type = "mixed", not_cont = c(rep(0, 3), rep(1, 10)), 
                                iter = Nsamples, save = TRUE, centrality = TRUE)

# save(bdgraph_fit, file = "R/BDgraph_example.RData")
load("R/BDgraph_example.RData")

# -----------------------------------------------------------------------------------------------------------------

# =========================
#  4. DATA VISUALISATION
# =========================

# 0. Plot posterior structure probability
pdf('Plot/Posterior_Structure_Probability.pdf', width=9, height=6.5)
par(cex.main = 2.5, mar = c(5, 8, 4, 5) + 0.1, mgp = c(3.5, 1, 0), cex.lab = 2.5 , font.lab = 2, cex.axis = 1.9, bty = "n", las = 1)
plot(sort(bdgraph_fit$graph_weights/sum(bdgraph_fit$graph_weights), decreasing = T), bty = "n", xlab = "", 
     ylab = " ", ylim = c(0, 0.06), xlim = c(0, 600), cex = 2, 
     pch = 1, axes = FALSE)
abline(h = 1/2^78, col = "grey", lty = 2, lwd = 4)
axis(1)
axis(2) 
par(las = 0)
mtext("Structure Index", side = 1, line = 3.5, cex = 2.5)
mtext("Posterior Structure Probability", side = 2, line = 5, cex = 2.5)
dev.off()

# 0.a) Plot posterior structure complexity probability
complexity <- c()
for(i in 1:length(bdgraph_fit$sample_graph)){
  complexity[i] <- sum(as.numeric(unlist(strsplit(bdgraph_fit$sample_graph[i], ""))))
}

data_complexity <- cbind(complexity, bdgraph_fit$graph_weights)  %>%
  as_tibble() %>% 
  group_by(complexity) %>% 
  summarise(complexity_weight = sum(V2))

prior_complexity <- c()
for (i in 1:78){
  prior_complexity[i] <- (1/2^78)*choose(78, i) 
}

# Plot posterior mass
pdf('Plot/Posterior_Complexity_Probability.pdf', width=9, height=6.5)
par(cex.main = 2.5, mar = c(5, 8, 4, 5) + 0.1, mgp = c(3.5, 1, 0), cex.lab = 2.5 , font.lab = 2, cex.axis = 1.9, bty = "n", las = 1)
plot(data_complexity$complexity, data_complexity$complexity_weight/sum(data_complexity$complexity_weight), bty = "n", 
     xlab = " ", 
     ylab = "", cex = 2, xlim = c(30, 60), ylim = c(0, 0.5),
     pch = 1, axes = FALSE)
points(1:78, prior_complexity, col = "grey", cex = 2, pch = 4)
axis(1)
axis(2) 
par(las = 0)
mtext("Structure Complexity", side = 1, line = 3.5, cex = 2.5)
mtext("Posterior Complexity Probability", side = 2, line = 5.5, cex = 2.5)
dev.off()


Layout <- averageLayout(bdgraph_fit$parameters*bdgraph_fit$structure)
Labels = c("stress", "anxiety", "depression", "extraverted", "critical", 
           "dependable", "anxious", "open, complex", 
           "reserved", "sympathetic", "disorganized", "calm", 
           "conventional")
Groups = c(rep("Mental Health", 3), rep("Personality", 10))

# 1. Plot structure
easybgm::plot_structure(bdgraph_fit, layout = Layout, legend = T, nodeNames = Labels, 
                        groups = Groups, layoutScale = c(.8,1), palette = "R", color= c("#fbb20a","#AF601A" ),
                        theme = "TeamFortress", vsize = 6)

# 2. Plot median probability model
pdf('Plot/Network_Median_Prob.pdf', width=30, height=17)
easybgm::plot_network(bdgraph_fit, layout = Layout, legend = T, nodeNames = Labels,
                      groups = Groups,  layoutScale = c(.9,1), palette = "R", color= c("#fbb20a","#E59866" ),
                      theme = "TeamFortress", vsize = 7, legend.cex = 2.5, layoutOffset = c(-.2, 0), label.cex = 1.1)
dev.off()


pdf('Plot/Network_Median_Prob_nolegend.pdf', width=30, height=20)
easybgm::plot_network(bdgraph_fit, layout = Layout, legend = F, nodeNames = Labels,
                      groups = Groups,  layoutScale = c(.8,1), palette = "R", color= c("#fbb20a","#E59866" ),
                      theme = "TeamFortress", vsize = 9, legend.cex = 1, layoutOffset = c(-.3, 0), label.cex = 1.3)
dev.off()

# 3. Plot evidence plot
colnames(m) <- rownames(m) <- colnames(data)
pdf('Plot/Network_Evidence_Plot.pdf', width=20, height=20)
easybgm::plot_edgeevidence(bdgraph_fit, split = T, edge.width = 3, 
                           layout = Layout, legend = F, nodeNames = Labels,
                           groups = Groups, palette = "R", color= c("#fbb20a","#E59866"),
                           theme = "TeamFortress", vsize = 9, legend.cex = 2.3,  label.cex = 1.3, 
                           repulsion = 0.5)
dev.off()


# 4. Plot density of posterior parameter estimates
p <- ncol(data)
Prior_samples = matrix(0, nrow = 10000, ncol = (p*(p-1))/2)
for(i in 1:10000){
  adj <- graph.sim(p= 13, graph = "random")
  samples <- rgwish(n = 1, b = 3, adj = bdgraph_fit$structure)
  wishart_samples <- wi2net(samples)
  Prior_samples[i,] <- as.vector(wishart_samples[upper.tri(wishart_samples)])
}


pdf('Plot/Density_Included.pdf', width=9, height=6.5)
par(cex.main = 2.5, mar = c(5, 5, 4, 2) + 0.1, mgp = c(3.5, 1, 0), cex.lab = 2.5 , font.lab = 2, cex.axis = 1.9, bty = "n", las = 1)
plot(density(bdgraph_fit$samples_posterior[, 7]), bty = "n", 
     xlab = " ", main = "",
     ylab = "", cex = 2, lwd = 2,
     pch = 1, axes = FALSE, xlim = c(0, 0.2))
#lines(density(Prior_samples[, 3]), lty = 2)
axis(1)
axis(2, labels = FALSE, lwd.ticks = 0)
par(las = 0)
mtext(expression(theta[ij]), side = 1, line = 3.5, cex = 2.5)
mtext("Density", side = 2, line = 2, cex = 2.5, las = 0)
dev.off()

pdf('Plot/Density_Inconclusive.pdf', width=9, height=6.5)
par(cex.main = 2.5, mar = c(5, 5, 4, 2) + 0.1, mgp = c(3.5, 1, 0), cex.lab = 2.5 , 
    cex.axis = 1.9, bty = "n", las = 1)
plot(density(bdgraph_fit$samples_posterior[, 11]), bty = "n", 
     xlab = " ", main = "",
     ylab = "", cex = 2, lwd = 2,
     pch = 1, axes = FALSE, xlim = c(0, 0.2))
axis(1)
axis(2, labels = FALSE, lwd.ticks = 0)
par(las = 0)
mtext(expression(theta[ij]), side = 1, line = 3.5, cex = 2.5)
mtext("Density", side = 2, line = 2, cex = 2.5,  las = 0)
dev.off()

# 5. Centrality plot
pdf('Plot/Centrality.pdf', width=11, height=6.5)
plot_centrality(bdgraph_fit) +
  geom_point(size = 5, col = c(rep("#fbb20a", 3), rep("#E59866", 10))) +
  geom_errorbar(aes(ymin = as.numeric(lower), ymax = as.numeric(upper)), 
                size = 1, width = .2, col = c(rep("#fbb20a", 3), rep("#E59866", 10))) +
  theme_apa() +
  labs(y = "Strength Centrality", x = "Node") +
  theme(text = element_text(size = 24), axis.title.x = element_text(size=22), 
        axis.title.y = element_text(size=22), panel.border = element_blank(),
        axis.line = element_line(colour = "black", size = 1.1), axis.ticks.length=unit(.3, "cm"), 
        legend.position = c(.85, 0.9), legend.text = element_text(size = 20), 
        axis.text.x = element_text(size=20, angle=75, vjust=0.5)) +
  scale_color_manual(values = c( "#ebac00", "#142f69")) 
dev.off()

