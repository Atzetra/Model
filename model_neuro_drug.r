rm(list = ls())
library(deSolve)
library(rootSolve)
library(tidyverse)

#################
# Model parameters
ks <- 0.2054
# synthesis rate of NAPE
kdegnape <- 0.001
# degredation rate of NAPE
kdeg <- 0.001
# degradation rate of PA
kpa <- 0.43
# basal synthesis of PA
kpeasyn <- 0.126
# basal synthesis of PEA
kpeadeg <- 0.13
# basal degredation of PEA
kass <- 0.001
# assosiation constant of PPAR
kdis <- 0.005
# dissasociation of PPARa
Vmpld <- 97.9 * 0.94
# vmax for NAPE to PEA conversion by PLD enzyme
Kmpld <- 3300
# Km for NAPE to PEA conversion by PLD enzyme
Vmgde <- 3 * 1.24
# vmax for NAPE to PEA conversion by GDE1_4 enzyme
Kmgde <- 12000000
# Km for NAPE to PEA conversion by GDE1_4 enzyme
Vmfa <- 2.6 * 550 * 1.04
# vmax for PEA to PA conversion by FAAH enzyme
Kmfa <- 5000
# Km for PEA to PA conversion by FAAH enzyme
Vmna <- 12 * 550 * 0.91
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
        dPEA <- kpeasyn + Vmpld * NAPE / (Kmpld + NAPE) + Vmgde * NAPE / (Kmgde + NAPE) - Vmfa * PEA / (Kmfa + PEA) - Vmna * PEA / (Kmna + PEA) - kass * PEA + kdis * PPRA - kpeadeg * PEA
        dPA <- kpa + Vmfa * PEA / (Kmfa + PEA) + Vmna * PEA / (Kmna + PEA) - kdeg * PA
        dPPRA <- kass * PEA - kdis * PPRA

        res <- c(dNAPE, dPEA, dPA, dPPRA)
        list(res)
    })
}

t <- seq(0, 72, 1)

xstart <- c(
    NAPE = 6.7113176,
    PEA = 0.6714116,
    PA = 600.4321418,
    PPRA = 0.6714116,
)


rs_steadystate <- runsteady(y = xstart, func = PEA_model_steady, parms = para_steady, times = c(0, 1e18))

(df_pea_levels <- rbind(rs_steadystate$y))

PEA_model_no_faah <- function(t, x, parms) {
    with(as.list(c(parms, x)), {
        dNAPE <- ks - Vmpld * NAPE / (Kmpld + NAPE) - Vmgde * NAPE / (Kmgde + NAPE) - kdegnape * NAPE
        dPEA <- kpeasyn + Vmpld * NAPE / (Kmpld + NAPE) + Vmgde * NAPE / (Kmgde + NAPE) - Vmna * PEA / (Kmna + PEA) - kass * PEA + kdis * PPRA - kpeadeg * PEA
        dPA <- kpa + Vmna * PEA / (Kmna + PEA) + Vmfa * PEA / (Kmfa + PEA) - kdeg * PA
        dPPRA <- kass * PEA - kdis * PPRA

        res <- c(dNAPE, dPEA, dPA, dPPRA)
        list(res)
    })
}

para_no_faah <- unlist(c(data.frame(
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
    Kmgde, # Km for NAPE to PEA conversion by GDE1_4 enzyme
    Vmna, # vmax for PEA to PA conversion by NAAA enzyme
    Kmna,
    Vmfa,
    Kmfa
)))
para_no_faah

rs_no_faah <- runsteady(y = xstart, func = PEA_model_no_faah, parms = para_no_faah, times = c(0, 1e18))
(df_pea_levels <- rbind(df_pea_levels, rs_no_faah$y))

PEA_model_no_naaa <- function(t, x, parms) {
    with(as.list(c(parms, x)), {
        dNAPE <- ks - Vmpld * NAPE / (Kmpld + NAPE) - Vmgde * NAPE / (Kmgde + NAPE) - kdegnape * NAPE
        dPEA <- kpeasyn + Vmpld * NAPE / (Kmpld + NAPE) + Vmgde * NAPE / (Kmgde + NAPE) - Vmfa * PEA / (Kmfa + PEA) - kass * PEA + kdis * PPRA - kpeadeg * PEA
        dPA <- kpa + Vmfa * PEA / (Kmfa + PEA) - kdeg * PA
        dPPRA <- kass * PEA - kdis * PPRA

        res <- c(dNAPE, dPEA, dPA, dPPRA)
        list(res)
    })
}

para_no_naaa <- unlist(c(data.frame(
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
    Kmgde, # Km for NAPE to PEA conversion by GDE1_4 enzyme
    Vmfa, # vmax for PEA to PA conversion by FAAH enzyme
    Kmfa # Km for PEA to PA conversion by FAAH enzyme
)))
para_no_naaa

rs_no_naaa <- runsteady(y = xstart, func = PEA_model_no_naaa, parms = para_no_naaa, times = c(0, 1e18))
(df_pea_levels <- rbind(df_pea_levels, rs_no_naaa$y))

PEA_model_no_both <- function(t, x, parms) {
    with(as.list(c(parms, x)), {
        dNAPE <- ks - Vmpld * NAPE / (Kmpld + NAPE) - Vmgde * NAPE / (Kmgde + NAPE) - kdegnape * NAPE
        dPEA <- kpeasyn + Vmpld * NAPE / (Kmpld + NAPE) + Vmgde * NAPE / (Kmgde + NAPE) - kpeadeg * PEA
        dPA <- kpa - kdeg * PA
        dPPRA <- kass * PEA - kdis * PPRA

        res <- c(dNAPE, dPEA, dPA, dPPRA)
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
df_pea_levels$model <- c("steadystate", "no_faaa", "no_naaa", "no_both")
df_pea_levels


f <- ggplot(df_pea_levels, aes(model, PEA))
f + geom_col()

out1 <- ode(y = xstart, func = PEA_model_steady, parms = para_steady, times = t)
out2 <- ode(y = xstart, func = PEA_model_no_both, parms = para_no_both, times = t)


out_comb <- out2/out1


out_comb[,'time'] <- out1[,'time']

out_comb

plot(out_comb)
df_pea_levels

# Multiply the Vmax by the fold change of the neuro cell lines