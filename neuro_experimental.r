rm(list = ls())
graphics.off()
library(deSolve)
library(tidyverse)
library(scales)
library(rootSolve)
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

#############################################
#                   MODELS                  #
#############################################

t <- seq(0, 96, 1)

xstart <- c(
    NAPE = 6.710947,
    PEA = 0.6798838,
    PA = 676.4779672,
    PPAR = 0.6798813,
    DRUG = 1
)

# No inhibition model

PEA_model_steady <- function(t, x, parms) {
    with(as.list(c(parms, x)), {
        if (t < 24 & DRUG >= 0) {
            dDRUG <- -0.8
        } else if (t > 24 & DRUG <= 1) {
            dDRUG <- 0.05
        } else {
            dDRUG <- 0
        }
        dNAPE <- ks - Vmpld * NAPE / (Kmpld + NAPE) - Vmgde * NAPE / (Kmgde + NAPE) - kdegnape * NAPE
        dPEA <- kpeasyn + Vmpld * NAPE / (Kmpld + NAPE) + Vmgde * NAPE / (Kmgde + NAPE) - Vmfa * PEA / (Kmfa + PEA) - Vmna * PEA / (Kmna + PEA) - kass * PEA + kdis * PPAR - kpeadeg * PEA - kpeaextra * PEA
        dPA <- kpa + Vmfa * PEA / (Kmfa + PEA) + Vmna * PEA / (Kmna + PEA) + kpeaextra * PEA - kdeg * PA
        dPPAR <- kass * PEA - kdis * PPAR

        res <- c(dNAPE, dPEA, dPA, dPPAR, dDRUG)
        list(res)
    })
}

# FAAH + NAAA Inhibition

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

PEA_model_no_faah <- function(t, x, parms) {
    with(as.list(c(parms, x)), {
        if (t < 24 & DRUG >= 0) {
            dDRUG <- -0.12
        } else if (t > 24 & DRUG <= 1) {
            dDRUG <- 0.017
        } else {
            dDRUG <- 0
        }
        dNAPE <- ks - Vmpld * NAPE / (Kmpld + NAPE) - Vmgde * NAPE / (Kmgde + NAPE) - kdegnape * NAPE
        dPEA <- kpeasyn + Vmpld * NAPE / (Kmpld + NAPE) + Vmgde * NAPE / (Kmgde + NAPE) - DRUG * Vmfa * PEA / (Kmfa + PEA) - Vmna * PEA / (Kmna + PEA) - kpeadeg * PEA - kpeaextra * PEA - kass * PEA + kdis * PPAR
        dPA <- kpa + DRUG * Vmfa * PEA / (Kmfa + PEA) +Vmna * PEA / (Kmna + PEA) - kdeg * PA
        dPPAR <- kass * PEA - kdis * PPAR

        res <- c(dNAPE, dPEA, dPA, dPPAR, dDRUG)
        list(res)
    })
}

PEA_model_no_naaa <- function(t, x, parms) {
    with(as.list(c(parms, x)), {
        if (t < 24 & DRUG >= 0) {
            dDRUG <- -0.12
        } else if (t > 24 & DRUG <= 1) {
            dDRUG <- 0.017
        } else {
            dDRUG <- 0
        }
        dNAPE <- ks - Vmpld * NAPE / (Kmpld + NAPE) - Vmgde * NAPE / (Kmgde + NAPE) - kdegnape * NAPE
        dPEA <- kpeasyn + Vmpld * NAPE / (Kmpld + NAPE) + Vmgde * NAPE / (Kmgde + NAPE) - Vmfa * PEA / (Kmfa + PEA) - DRUG * Vmna * PEA / (Kmna + PEA) - kpeadeg * PEA - kpeaextra * PEA - kass * PEA + kdis * PPAR
        dPA <- kpa + Vmfa * PEA / (Kmfa + PEA) + DRUG * Vmna * PEA / (Kmna + PEA) - kdeg * PA
        dPPAR <- kass * PEA - kdis * PPAR

        res <- c(dNAPE, dPEA, dPA, dPPAR, dDRUG)
        list(res)
    })
}

#############################################
#                  RUN ODEs                 #
#############################################

out_none <- ode(y = xstart, func = PEA_model_steady, parms = para_steady, times = t, method = "ode45")
out_no_faah <- ode(y = xstart, func = PEA_model_no_faah, parms = para_steady, times = t, method = "ode45")
out_no_naaa <- ode(y = xstart, func = PEA_model_no_naaa, parms = para_steady, times = t, method = "ode45")
out_no_both <- ode(y = xstart, func = PEA_model_no_both, parms = para_steady, times = t, method = "ode45")

# Calculate relative difference compared to no_inhibition
# And round numbers to prevent strange plotting

out_none <- as.data.frame(out_none)
out_none <- out_none %>% mutate(model = "None",
relative_pea = PEA / out_none$PEA,
relative_ppar = PPAR / out_none$PPAR,
relative_pa = PA / out_none$PA,
relative_nape = NAPE / out_none$NAPE,
across(starts_with("relative"), round, 4))

out_none

out_no_faah <- as.data.frame(out_no_faah)
out_no_faah <- out_no_faah %>% mutate(model = "FAAH",
relative_pea = PEA / out_none$PEA,
relative_ppar = PPAR / out_none$PPAR,
relative_pa = PA / out_none$PA,
relative_nape = NAPE / out_none$NAPE,
across(starts_with("relative"), round, 4))

out_no_naaa <- as.data.frame(out_no_naaa)
out_no_naaa <- out_no_naaa %>% mutate(model = "NAAA",
relative_pea = PEA / out_none$PEA,
relative_ppar = PPAR / out_none$PPAR,
relative_pa = PA / out_none$PA,
relative_nape = NAPE / out_none$NAPE,
across(starts_with("relative"), round, 4))

out_no_both <- as.data.frame(out_no_both)
out_no_both <- out_no_both %>% mutate(model = "FAAH + NAAA",
relative_pea = PEA / out_none$PEA,
relative_ppar = PPAR / out_none$PPAR,
relative_pa = PA / out_none$PA,
relative_nape = NAPE / out_none$NAPE,
across(starts_with("relative"), round, 4))

# Combine data frames

combined_df <- rbind(out_none, out_no_faah, out_no_naaa, out_no_both)

combined_df


#############################################
#                   GGPLOT                  #
#############################################

#############################################
#                  Line plot                #
#############################################

plot_PEA <- ggplot(combined_df, aes(time, relative_pea, color = model)) +
    geom_line() +
    scale_y_continuous(
        labels = label_number_auto()
    ) +
    labs(x = "Time (h)", y = "PEA Change", color = "Inhibited pathway") +
    scale_color_manual(values = c("#000063", "#B07312", "#633F05", "#A5B012"))

plot_NAPE <- ggplot(combined_df, aes(time, relative_nape, color = model)) +
    geom_line() +
    scale_y_continuous(
        labels = label_number(accuracy = 0.02)
    ) +
    labs(x = "Time (h)", y = "NAPE Change", color = "Inhibited pathway") +
    scale_color_manual(values = c("#000063", "#B07312", "#633F05", "#A5B012"))

plot_PPAR <- ggplot(combined_df, aes(time, relative_ppar, color = model)) +
    geom_line() +
    scale_y_continuous(
        labels = label_number_auto()
    ) +
    labs(x = "Time (h)", y = "PPAR Change", color = "Inhibited pathway") +
    scale_color_manual(values = c("#000063", "#B07312", "#633F05", "#A5B012"))

plot_PA <- ggplot(combined_df, aes(time, relative_pa, color = model)) +
    geom_line() +
    scale_y_continuous(
        labels = label_number_auto()
    ) +
    labs(x = "Time (h)", y = "PA Change", color = "Inhibited pathway") +
    scale_color_manual(values = c("#000063", "#B07312", "#633F05", "#A5B012"))

plot_all <- ggarrange(plot_NAPE, plot_PEA, plot_PPAR, plot_PA,
    ncol = 2,
    nrow = 2,
    common.legend = TRUE
)
plot_all