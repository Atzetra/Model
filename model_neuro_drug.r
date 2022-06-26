rm(list = ls())
graphics.off()
library(deSolve)
library(tidyverse)
library(rootSolve)
library(scales)
library(ggpubr)

#################
# Model parameters
ks <- 0.2054
# synthesis rate of NAPE
kdegnape <- 0.001
# degredation rate of NAPE
kdeg <- 0.001
# degradation rate of PA
kpa <- 0.4
# basal synthesis of PA
kpeasyn <- 0.155
# basal synthesis of PEA
kpeadeg <- 0.048
# basal degredation of PEA
kpeaextra <- 0.08
# extra pathway of PEA
kass <- 0.001
# assosiation constant of PPAR
kdis <- 0.001
# dissasociation of PPARa
Vmpld <- 97.9 * 0.938
# vmax for NAPE to PEA conversion by PLD enzyme
Kmpld <- 3300
# Km for NAPE to PEA conversion by PLD enzyme
Vmgde <- 3 * 1.236
# vmax for NAPE to PEA conversion by GDE1_4 enzyme
Kmgde <- 12000000
# Km for NAPE to PEA conversion by GDE1_4 enzyme
Vmfa <- 2.6 * 620 * 1.039
# vmax for PEA to PA conversion by FAAH enzyme
Kmfa <- 5000
# Km for PEA to PA conversion by FAAH enzyme
Vmna <- 12 * 620 * 0.910
# vmax for PEA to PA conversion by NAAA enzyme
Kmna <- 97000
# Km for PEA to PA conversion by NAAA enzyme

para_steady <- unlist(c(data.frame(
    ks,
    kdegnape,
    kdeg,
    kass,
    kdis,
    kpa,
    kpeasyn,
    kpeaextra,
    kpeadeg,
    Vmpld, # vmax for NAPE to PEA conversion by PLD enzyme
    Kmpld, # Km for NAPE to PEA conversion by PLD enzyme
    Vmgde, # vmax for NAPE to PEA conversion by GDE1_4 enzyme
    Kmgde, # Km for NAPE to PEA conversion by GDE1_4 enzyme
    Vmfa, # vmax for PEA to PA conversion by FAAH enzyme
    Kmfa, # Km for PEA to PA conversion by FAAH enzyme
    Vmna, # vmax for PEA to PA conversion by NAAA enzyme
    Kmna # Km for PEA to PA conversion by NAAA enzyme
)))

# Use discrete functions in order to capture end of effect with NAAA and FAAH

PEA_model_steady_no_drug <- function(t, x, parms) {
    with(as.list(c(parms, x)), {
        dNAPE <- ks - Vmpld * NAPE / (Kmpld + NAPE) - Vmgde * NAPE / (Kmgde + NAPE) - kdegnape * NAPE
        dPEA <- kpeasyn + Vmpld * NAPE / (Kmpld + NAPE) + Vmgde * NAPE / (Kmgde + NAPE) - Vmfa * PEA / (Kmfa + PEA) - Vmna * PEA / (Kmna + PEA) - kass * PEA + kdis * PPAR - kpeadeg * PEA - kpeaextra * PEA
        dPA <- kpa + Vmfa * PEA / (Kmfa + PEA) + Vmna * PEA / (Kmna + PEA) - kdeg * PA
        dPPAR <- kass * PEA - kdis * PPAR

        res <- c(dNAPE, dPEA, dPA, dPPAR)
        list(res)
    })
}

t <- seq(0, 72, 1)

xstart <- c(
    NAPE = 6.710947,
    PEA = 0.6798838,
    PA = 666.3042,
    PPAR = 0.6798813,
    DRUG = 1
)

xstart_no_drug <- c(
    NAPE = 6.710947,
    PEA = 0.6798838,
    PA = 666.3042,
    PPAR = 0.6798813
)

#############################################
#               COMPARING STEADY            #
#############################################

# Get neuro steady state levels
rs_steady_neuro <- runsteady(y = xstart_no_drug, func = PEA_model_steady_no_drug, parms = para_steady, times = c(0, 1e18))
neuro_steady <- as.data.frame(rs_steady_neuro$y)
neuro_steady <- tibble::rownames_to_column(neuro_steady)
names(neuro_steady) <- c("name", "level")
neuro_steady <- neuro_steady %>%
    pivot_wider(
    names_from = name,
    values_from = level
    )
neuro_steady$model <- "Neuro2a"

# Load RAW264.7 steady states
load(file = "Data/SteadyState/RAW.rdata")
raw_steady <- as.data.frame(rs_steady_raw$y)
raw_steady <- tibble::rownames_to_column(raw_steady)
names(raw_steady) <- c("name", "level")
raw_steady <- raw_steady %>%
    pivot_wider(
    names_from = name,
    values_from = level
    )
raw_steady$model <- "RAW264.7"

# Combine RAW and Neuro data frames
raw_compare <- rbind(neuro_steady, raw_steady)
raw_compare

#############################################
#             RUNNING NEURO ODE             #
#############################################

PEA_model_no_both <- function(t, x, parms) {
    with(as.list(c(parms, x)), {
        if (t < 24 & DRUG >= 0) {
            dDRUG <- -0.12
        } else if (t > 24 & DRUG <= 1) {
            dDRUG <- 0.017
        } else {
            dDRUG <- 0
        }
        dNAPE <- ks - Vmpld * NAPE / (Kmpld + NAPE) - Vmgde * NAPE / (Kmgde + NAPE) - kdegnape * NAPE
        dPEA <- kpeasyn + Vmpld * NAPE / (Kmpld + NAPE) + Vmgde * NAPE / (Kmgde + NAPE) - DRUG * Vmfa * PEA / (Kmfa + PEA) - DRUG * Vmna * PEA / (Kmna + PEA) - kpeadeg * PEA - kpeaextra * PEA - kass * PEA + kdis * PPAR
        dPA <- kpa + DRUG * Vmfa * PEA / (Kmfa + PEA) + DRUG * Vmna * PEA / (Kmna + PEA) - kdeg * PA
        dPPAR <- kass * PEA - kdis * PPAR

        res <- c(dNAPE, dPEA, dPA, dPPAR, dDRUG)
        list(res)
    })
}

para_no_both <- unlist(c(data.frame(
    ks,
    kdegnape,
    kdeg,
    kass,
    kdis,
    kpa,
    kpeasyn,
    kpeadeg,
    Vmpld, # vmax for NAPE to PEA conversion by PLD enzyme
    Kmpld, # Km for NAPE to PEA conversion by PLD enzyme
    Vmgde, # vmax for NAPE to PEA conversion by GDE1_4 enzyme
    Kmgde # Km for NAPE to PEA conversion by GDE1_4 enzyme
)))

outNeuro <- ode(y = xstart, func = PEA_model_no_both, parms = para_no_both, times = t, method = "ode45")

outNeuro <- as.data.frame(outNeuro)
outNeuro <- outNeuro %>%
    add_column(model = "Neuro2a")
outNeuro

#############################################
#                   GGPLOT                  #
#############################################

#############################################
#                 Column plot               #
#############################################

steady_compare_NAPE <- ggplot(raw_compare, aes(model, NAPE, fill = model)) +
    geom_col() +
    scale_fill_manual(values = c("#000063", "#B07312")) +
    theme(aspect.ratio = 20 / 9,
    legend.position = "none",
    axis.text.x = element_blank(),
    axis.title.x = element_blank()) +
    labs(y = "NAPE (fmol)", fill = "Model")

steady_compare_PEA <- ggplot(raw_compare, aes(model, PEA, fill = model)) +
    geom_col() +
    scale_fill_manual(values = c("#000063", "#B07312")) +
    theme(aspect.ratio = 20 / 9,
    legend.position = "none",
    axis.text.x = element_blank(),
    axis.title.x = element_blank()) +
    labs(y = "PEA (fmol)", fill = "Model")

steady_compare_PPAR <- ggplot(raw_compare, aes(model, PPAR, fill = model)) +
    geom_col() +
    scale_fill_manual(values = c("#000063", "#B07312")) +
    theme(aspect.ratio = 20 / 9,
    legend.position = "none",
    axis.text.x = element_blank(),
    axis.title.x = element_blank()) +
    labs(y = "PPAR binding", fill = "Model")

steady_compare_PA <- ggplot(raw_compare, aes(model, PA, fill = model)) +
    geom_col() +
    scale_fill_manual(values = c("#000063", "#B07312")) +
    theme(aspect.ratio = 20 / 9,
    legend.position = "none",
    axis.text.x = element_blank(),
    axis.title.x = element_blank()) +
    labs(y = "PA (fmol)", fill = "Model")

plot_col_total <- ggarrange(steady_compare_NAPE, steady_compare_PEA,
steady_compare_PPAR, steady_compare_PA,
ncol = 2,
nrow = 2,
common.legend = TRUE)

#############################################
#                Line plot                  #
#############################################

# Read RAW264 data from model_drug.r
load("Data/SteadyState/RAWout2.rdata")
raw1 <- out2
raw1 <- raw1 %>%
    add_column(model = "RAW264.7")

# Combine RAW and Neuro data
comb_line_data <- rbind(outNeuro, raw1)
comb_line_data

# Make plots

plot_NAPE <- ggplot(comb_line_data, aes(time, NAPE, color = model)) +
    geom_line() +
    scale_y_continuous(
        labels = label_number_auto()
    ) +
    labs(x = "Time (h)", y = "NAPE (fmol)", color = "Model") +
    scale_color_manual(values = c("#000063", "#B07312"))

plot_PEA <- ggplot(comb_line_data, aes(time, PEA, color = model)) +
    geom_line() +
    scale_y_continuous(
        labels = label_number_auto()
    ) +
    labs(x = "Time (h)", y = "PEA (fmol)", color = "Model") +
    scale_color_manual(values = c("#000063", "#B07312"))

plot_PPAR <- ggplot(comb_line_data, aes(time, PPAR, color = model)) +
    geom_line() +
    scale_y_continuous(
        labels = label_number_auto()
    ) +
    labs(x = "Time (h)", y = "PPAR binding", color = "Model") +
    scale_color_manual(values = c("#000063", "#B07312"))

plot_PA <- ggplot(comb_line_data, aes(time, PA, color = model)) +
    geom_line() +
    scale_y_continuous(
        labels = label_number_auto()
    ) +
    labs(x = "Time (h)", y = "PA (fmol)", color = "Model") +
    scale_color_manual(values = c("#000063", "#B07312"))

plot_all <- ggarrange(plot_NAPE, plot_PEA, plot_PPAR, plot_PA,
    ncol = 2,
    nrow = 2,
    common.legend = TRUE
)

#############################################
#                Combine plot               #
#############################################
plot_complete <- ggarrange(
plot_col_total, plot_all,
ncol = 2,
nrow = 1,
labels = c("A", "B")
)

plot_complete