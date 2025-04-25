
# ============================================================================================== # 
#  Tutorial for Bayesian cross-sectional network analysis: BDgraph, BGGM, rbinnet                #
# ============================================================================================== # 

# ---------------------------------------------------------------------------------------------- #
# ---------- TABLE OF CONTENTS ----------------------------------------------------------------- #
# ---------------------------------------------------------------------------------------------- #
# ---------- 01. Essential Procedures ---------------------------------------------------------- #
# ---------- 02. Data Loading ------------------------------------------------------------------ # 
# ---------- 03. BDgraph ----------------------------------------------------------------------- #
# ---------- 04. BGGM -------------------------------------------------------------------------- #
# ---------- 05. bgms ----------------------------------------------------------------------- #
# ---------- 06. Data Visualization ------------------------------------------------------------ #
# ---------- 07. ... --------------------------------------------------------------------------- #

# For details on the background of the methods, please refer to the 
# tutorial paper (see preprint on OSF https://osf.io/xjm6z/?view_only=ca0ed443d4674c4dbe42a09a17c8d463).
# Author: k.huth@uva.nl
# Update: 14.03.2025
# For the most up-to-date functionalities of easybgm, see https://github.com/KarolineHuth/easybgm

# -----------------------------------------------------------------------------------------------------------------

# =========================
#  1. ESSENTIAL PROCEDURES
# =========================

# 1. Set Working Directory
setwd("~/project") # Specify your working directory here

# 2. Install/Load the necessary packages
# install.packages("easybgm")
library(easybgm)
library(qgraph)
library(tidyverse)
library(readr)

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

bdgraph_res <- easybgm::easybgm(data = data, type = "mixed", 
                           package = "BDgraph", not_cont = c(rep(0, 3), rep(1, 10)), 
                           centrality = TRUE, save = TRUE, iter = Nsamples)

# -----------------------------------------------------------------------------------------------------------------

# =========================
#  4. BGGM
# =========================

bggm_res <- easybgm::easybgm(data = data, type = "mixed", 
                       package = "BGGM", not_cont = c(rep(0, 3), rep(1, 10)), 
                       centrality = TRUE, save = TRUE, iter = Nsamples)

# -----------------------------------------------------------------------------------------------------------------

# =========================
#  5. bgms
# =========================
# Please note, this package only allows for ordinal data estimation

bgms_res <- easybgm::easybgm(data = data[, 4:13], type = "ordinal", 
                    package = "bgms", 
                    centrality = TRUE, save = TRUE, iter = Nsamples)

# -----------------------------------------------------------------------------------------------------------------

# =========================
#  6. DATA VISUALISATION
# =========================
# Note: Exchange [Package]_res with the respective package that you used
res <- bgms_res

# 0. Plot posterior structure estimates 

# 0.a Posterior structure probability
plot_structure_probabilities(res, as_BF = F) #(NOTE: not usable with BGGM)

# 0.b. Posterior complexity probability
plot_complexity_probabilities(res) #(NOTE: not usable with BGGM)

# 1. Plot structure
plot_structure(res, layoutScale = c(.8,1), palette = "R",
               theme = "TeamFortress", vsize = 6, edge.width = .3, layout = "spring")

# 2. Plot network model (also referred to as median probability model)
plot_network(res, layout = "spring", 
             layoutScale = c(.8,1), palette = "R",
             theme = "TeamFortress", vsize = 6)

# 3. Plot evidence plot
plot_edgeevidence(res, edge.width = 2, split = F, legend = F)

# 4. Plot parameter forest plot  (Note: posterior_samples in bgm_extract needs to be set to TRUE)
plot_parameterHDI(res)

# 5. Plot strength centrality estimate + HDI 
plot_centrality(res) 



