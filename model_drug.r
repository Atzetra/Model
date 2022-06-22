rm(list = ls())
graphics.off()
library(deSolve)
library(tidyverse)
library(scales)

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
Vmpld <- 97.9
# vmax for NAPE to PEA conversion by PLD enzyme
Kmpld <- 3300
# Km for NAPE to PEA conversion by PLD enzyme
Vmgde <- 3
# vmax for NAPE to PEA conversion by GDE1_4 enzyme
Kmgde <- 12000000
# Km for NAPE to PEA conversion by GDE1_4 enzyme
Vmfa <- 2.6 * 620
# vmax for PEA to PA conversion by FAAH enzyme
Kmfa <- 5000
# Km for PEA to PA conversion by FAAH enzyme
Vmna <- 12 * 620
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

# Use decrete functions in order to capture end of effect with NAAA and FAAH

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
        dPA <- kpa + Vmfa * PEA / (Kmfa + PEA) + Vmna * PEA / (Kmna + PEA) - kdeg * PA
        dPPAR <- kass * PEA - kdis * PPAR

        res <- c(dNAPE, dPEA, dPA, dPPAR, dDRUG)
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


PEA_model_no_both <- function(t, x, parms) {
    with(as.list(c(parms, x)), {
        if (t < 24 & DRUG >= 0) {
            dDRUG <- -0.12
        } else if (t > 24 & DRUG <= 0.2) {
            dDRUG <- 0.002
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

out1 <- ode(y = xstart, func = PEA_model_steady, parms = para_steady, times = t, method = "ode45")
out2 <- ode(y = xstart, func = PEA_model_no_both, parms = para_no_both, times = t, method = "ode45")

out1 <- as.data.frame(out1)
out1 <- out1 %>%
    add_column(model = "None")
out1
out2 <- as.data.frame(out2)
out2 <- out2 %>%
    add_column(model = "FAAH + NAAH")

combined_df <- rbind(out1, out2)

combined_df

#############################################
#                EXPERIMENTAL               #
#############################################

model <- c("None", "None", "None", "None", "FAAH + NAAH", "FAAH + NAAH", "FAAH + NAAH", "FAAH + NAAH")
PEA <- c(0.537, 0.514, 0.706, 0.432, 0.518, 1.003, 2.523, 1.349)
sd <- c(0.074, 0.065, 0.047, 0.050, 0.046, 0.114, 0.142, 1.175)
time <- c(0, 6, 24, 48, 0, 6, 24, 48)

df_experimental <- data.frame(model, PEA, sd, time)

df_steady.PEA <- c(0.671, 2.523)
df_steady.sd <- c(0.015, 0.142)
df_steady.model <- c("None", "FAAH + NAAA")
df_steady <- data.frame(df_steady.model, df_steady.PEA, df_steady.sd)
names(df_steady) <- c("model", "PEA", "sd")

#############################################
#                   GGPLOT                  #
#############################################


#############################################
#                 Column plot               #
#############################################

line_plot <- ggplot(df_steady) +
    geom_col(aes(model, PEA, fill = model)) +
    geom_errorbar(aes(model, ymin = PEA - sd, ymax = PEA + sd)) +
    scale_fill_manual(values = c("#000063", "#B07312")) +
    theme(aspect.ratio = 20 / 9,
    legend.position = "none")

#############################################
#                  Line plot                #
#############################################

plot_PEA <- ggplot(combined_df, aes(time, PEA, color = model)) +
    geom_line(size = 1.2) +
    scale_y_continuous(
        labels = label_number_auto()
    ) +
    geom_point(data = df_experimental, aes(time, PEA, color = model)) +
    geom_errorbar(data = df_experimental,
    aes(x = time, ymin = PEA - sd, ymax = PEA + sd, color = model)) +
    labs(x = "Time", y = "PEA Levels", color = "Inhibited pathway") +
    scale_color_manual(values = c("#000063", "#B07312"))

plot_NAPE <- ggplot(combined_df, aes(time, NAPE, color = model)) +
    geom_line(size = 1.2) +
    scale_y_continuous(
        labels = label_number(accuracy = 0.02)
    ) +
    labs(x = "Time", y = "NAPE Levels", color = "Inhibited pathway") +
    scale_color_manual(values = c("#000063", "#B07312"))

plot_PPAR <- ggplot(combined_df, aes(time, PPAR, color = model)) +
    geom_line(size = 1.2) +
    scale_y_continuous(
        labels = label_number_auto()
    ) +
    labs(x = "Time", y = "PPAR Levels", color = "Inhibited pathway") +
    scale_color_manual(values = c("#000063", "#B07312"))

plot_PA <- ggplot(combined_df, aes(time, PA, color = model)) +
    geom_line(size = 1.2) +
    scale_y_continuous(
        labels = label_number_auto()
    ) +
    labs(x = "Time", y = "PA Levels", color = "Inhibited pathway") +
    scale_color_manual(values = c("#000063", "#B07312"))

plot_all <- ggarrange(plot_NAPE, plot_PEA, plot_PPAR, plot_PA,
    ncol = 2,
    nrow = 2,
    common.legend = TRUE
)

plot_steady <- ggarrange(line_plot, plot_all,
ncol = 2,
labels = c("A", "B"),
widths = c(1,2))
plot_steady

ggsave("figures/basemodel.pdf")