# Name:     02_fnct_con_div.R
# Updated:  20230123
# User:     Andrew Kleist

# packages, working directory
setwd("~/Desktop/Volkman_Lab/BABU_VOLKMAN/CHEMOKINE_PRJ/Data/chemokine_gpcr_ms/")
library(tidyverse)
library(ggrepel)
  
##### PART 1: CHEMOKINE FUNCTIONAL DISTANCE MATRIX #############################

data <- read_csv("06_network/output/chemokine_gpcr_interactome_validated.csv")
data$n <- c(1)
data <- data %>% pivot_wider(names_from = "ckr", values_from = n)
data[is.na(data)] <- 0
data <- as.matrix(data)
rownames(data) <- data[,1]
data <- data[, -1]

dist <- dist(data, method = "euclidean")
dist <- as.data.frame(as.matrix(dist))
dist$a <- rownames(dist)
dist <- dist %>% pivot_longer(cols = 1:(ncol(dist)-1), names_to = "b", values_to = "dist"  )
dist <- dist %>% select(a:dist)

dist <- dist %>% mutate(dist_nl = dist / max(dist))
dist <- dist %>% mutate(sim_nl = 1-dist_nl)
# temp <- dist
# colnames(temp) <- c("b", "a", "dist", "dist_nl", "sim_nl")
# dist <- rbind(dist, temp)
# rm(temp)
# dist <- dist %>% filter(ck1 != ck2)

write_csv(dist, "06_network/output/network_similarity_scores_ck.csv")


##### PART 2: RECEPTOR FUNCTIONAL DISTANCE MATRIX #############################

data <- read_csv("06_network/output/chemokine_gpcr_interactome_validated.csv")
data$n <- c(1)
data <- data %>% pivot_wider(names_from = "ckr", values_from = n)
data[is.na(data)] <- 0
data <- as.matrix(data)
rownames(data) <- data[,1]
data <- data[, -1]
data <- t(data)

dist <- dist(data, method = "euclidean")
dist <- as.data.frame(as.matrix(dist))
dist$a <- rownames(dist)
dist <- dist %>% pivot_longer(cols = 1:(ncol(dist)-1), names_to = "b", values_to = "dist"  )

dist <- dist %>% mutate(dist_nl = dist / max(dist))
dist <- dist %>% mutate(sim_nl = 1-dist_nl)
# dist <- dist %>% filter(a != b)

write_csv(dist, "06_network/output/network_similarity_scores_ckr.csv")


