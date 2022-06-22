rm(list = ls())
library(deSolve)
library(rootSolve)
library(tidyverse)
library("RColorBrewer")
library("svglite")
library("scales")

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
para_steady

# Use decrete functions in order to capture end of effect with NAAA and FAAH

PEA_model_steady <- function(t, x, parms) {
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
    PPAR = 0.6798813
)


rs_steadystate <- runsteady(y = xstart, func = PEA_model_steady, parms = para_steady, times = c(0, 1e18))

(df_pea_levels <- rbind(rs_steadystate$y))

PEA_model_no_both <- function(t, x, parms) {
    with(as.list(c(parms, x)), {
        dNAPE <- ks - Vmpld * NAPE / (Kmpld + NAPE) - Vmgde * NAPE / (Kmgde + NAPE) - kdegnape * NAPE
        dPEA <- kpeasyn + Vmpld * NAPE / (Kmpld + NAPE) + Vmgde * NAPE / (Kmgde + NAPE) - kpeadeg * PEA - kpeaextra * PEA - kass * PEA + kdis * PPAR
        dPA <- kpa - kdeg * PA
        dPPAR <- kass * PEA - kdis * PPAR

        res <- c(dNAPE, dPEA, dPA, dPPAR)
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

rs_no_both <- runsteady(y = xstart, func = PEA_model_no_both, parms = para_no_both, times = c(0, 1e18))
df_pea_levels <- rbind(df_pea_levels, rs_no_both$y)

df_pea_levels <- as.data.frame(df_pea_levels)
df_pea_levels$model <- c("steadystate", "no_both")
df_pea_levels


out1 <- ode(y = xstart, func = PEA_model_steady, parms = para_steady, times = t)
out2 <- ode(y = xstart, func = PEA_model_no_both, parms = para_no_both, times = t)


out_comb <- out2 / out1


out_comb[, "time"] <- out1[, "time"]

out_comb

df_pea_levels
saveRDS(df_pea_levels, file = "Data/SteadyState/df_pea_levels.rdata")

df_steadymeasure <- read.csv2("Experimental.csv")
df_steadymeasure
saveRDS(df_steadymeasure, file = "Data/SteadyState/df_steadymeasure.rdata")


f <- ggplot() +
    geom_col(data = df_pea_levels, aes(model, PEA, fill = model)) +
    geom_errorbar(data = df_steadymeasure, aes(ï..treatment, SD, ymin = level - SD, ymax = level + SD)) +
    theme(axis.text.x = element_text(angle = 60, hjust = 1), aspect.ratio = 4/3,
legend.position = "none") +
scale_x_discrete(labels = c("FAAH + NAAA", "None"))
f

save(f, file = "figures/ggplots/f.rdata")


# Multiply the Vmax by the fold change of the neuro cell lines

# i <- ggplot(as.data.frame(out_comb), aes(time, PEA))
# i + geom_line()
