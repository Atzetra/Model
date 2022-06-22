rm(list = ls())
library(tidyverse)
library(ggpubr)
library("RColorBrewer")
library("svglite")
theme_set(theme_pubr())


###############################################
# Palette
###############################################

# #000063, #B07312, #633F05, #A5B012, #5D6305

###############################################
# Steadystate
###############################################

# Import data for steadystate
df_pea_levels <- readRDS("Data/SteadyState/df_pea_levels.rdata")
df_steadymeasure <- readRDS("Data/SteadyState/df_steadymeasure.rdata")

load("figures/ggplots/f.rdata")
load("figures/ggplots/plot_all.rdata")


figure <- ggarrange(f, plot_all,
labels = c("A", "B"),
ncol = 2)
figure