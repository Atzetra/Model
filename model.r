rm(list = ls())
library(deSolve)
library(rootSolve)
library(tidyverse)

#################
# Model parameters
ks <- 0.1987
# synthesis rate of NAPE
kdeg <- 0.001
# degradation rate of PA
kpa <- 0.6
# basal synthesis of PA
kpea <- 0.2953
# basal degredation of PEA
kass <- 0.005
# assosiation constant of PPAR
kdis <- 0.005
# dissasociation of PPARa
Vmpld <- 97.9
# vmax for NAPE to PEA conversion by PLD enzyme
Kmpld <- 3300
# Km for NAPE to PEA conversion by PLD enzyme
Vmgde <- 3
# vmax for NAPE to PEA conversion by GDE1_4 enzyme
Kmgde <- 12000000
# Km for NAPE to PEA conversion by GDE1_4 enzyme
Vmfa <- 2.6
# vmax for PEA to PA conversion by FAAH enzyme
Kmfa <- 5000
# Km for PEA to PA conversion by FAAH enzyme
Vmna <- 12
# vmax for PEA to PA conversion by NAAA enzyme
Kmna <- 97000
# Km for PEA to PA conversion by NAAA enzyme

# Add basal synthesis of PEA

para_steady <- unlist(c(data.frame(
    ks,
    kdeg,
    kass,
    kdis,
    kpa,
    kpea,
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
        dNAPE <- ks - Vmpld * NAPE / (Kmpld + NAPE) - Vmgde * NAPE / (Kmgde + NAPE)
        dPEA <- Vmpld * NAPE / (Kmpld + NAPE) + Vmgde * NAPE / (Kmgde + NAPE) - Vmfa * PEA / (Kmfa + PEA) - Vmna * PEA / (Kmna + PEA) - kass * PEA + kdis * PPRA - kpea * PEA
        dPA <- kpa + Vmfa * PEA / (Kmfa + PEA) + Vmna * PEA / (Kmna + PEA) - kdeg * PA
        dPPRA <- kass * PEA - kdis * PPRA

        res <- c(dNAPE, dPEA, dPA, dPPRA)
        list(res)
    })
}

t <- seq(0, 48, 1)
xstart <- c(
    NAPE = 6.7113176,
    PEA = 0.6714116,
    PA = 600.4321418,
    PPRA = 0.6714116
)


rs_steadystate <- runsteady(y = xstart, func = PEA_model_steady, parms = para_steady, times = c(0, 1e18))

(df_pea_levels <- rbind(rs_steadystate$y))

PEA_model_no_faah <- function(t, x, parms) {
    with(as.list(c(parms, x)), {
        dNAPE <- ks - Vmpld * NAPE / (Kmpld + NAPE) - Vmgde * NAPE / (Kmgde + NAPE)
        dPEA <- Vmpld * NAPE / (Kmpld + NAPE) + Vmgde * NAPE / (Kmgde + NAPE) - Vmna * PEA / (Kmna + PEA) - kass * PEA + kdis * PPRA - kpea * PEA
        dPA <- kpa + Vmna * PEA / (Kmna + PEA) - kdeg * PA
        dPPRA <- kass * PEA - kdis * PPRA

        res <- c(dNAPE, dPEA, dPA, dPPRA)
        list(res)
    })
}

para_no_faah <- unlist(c(data.frame(
    ks,
    kdeg,
    kass,
    kdis,
    kpa,
    kpea,
    Vmpld, # vmax for NAPE to PEA conversion by PLD enzyme
    Kmpld, # Km for NAPE to PEA conversion by PLD enzyme
    Vmgde, # vmax for NAPE to PEA conversion by GDE1_4 enzyme
    Kmgde, # Km for NAPE to PEA conversion by GDE1_4 enzyme
    Vmna, # vmax for PEA to PA conversion by NAAA enzyme
    Kmna
)))
para_no_faah

rs_no_faah <- runsteady(y = xstart, func = PEA_model_no_faah, parms = para_no_faah, times = c(0, 1e18))
(df_pea_levels <- rbind(df_pea_levels, rs_no_faah$y))

PEA_model_no_naaa <- function(t, x, parms) {
    with(as.list(c(parms, x)), {
        dNAPE <- ks - Vmpld * NAPE / (Kmpld + NAPE) - Vmgde * NAPE / (Kmgde + NAPE)
        dPEA <- Vmpld * NAPE / (Kmpld + NAPE) + Vmgde * NAPE / (Kmgde + NAPE) - Vmfa * PEA / (Kmfa + PEA) - kass * PEA + kdis * PPRA - kpea * PEA
        dPA <- kpa + Vmfa * PEA / (Kmfa + PEA) - kdeg * PA
        dPPRA <- kass * PEA - kdis * PPRA

        res <- c(dNAPE, dPEA, dPA, dPPRA)
        list(res)
    })
}

para_no_naaa <- unlist(c(data.frame(
    ks,
    kdeg,
    kass,
    kdis,
    kpa,
    kpea,
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
df_pea_levels <- as.data.frame(df_pea_levels)
df_pea_levels$model <- c("steadystate", "no_faaa", "no_naaa")
df_pea_levels



f <- ggplot(df_pea_levels, aes(model, PEA))
f + geom_boxplot()

out1 <-  ode(y = xstart, func = PEA_model_steady, parms = para_steady, times = t)
plot(out1)
