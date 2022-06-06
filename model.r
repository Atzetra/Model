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



para <- unlist(c(data.frame(
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
para

PEA_model <- function(t, x, parms) {
    with(as.list(c(parms, x)), {
        dNAPE <- ks - Vmpld * NAPE / (Kmpld + NAPE) - Vmgde * NAPE / (Kmgde + NAPE)
        dPEA <- Vmpld * NAPE / (Kmpld + NAPE) + Vmgde * NAPE / (Kmgde + NAPE) - Vmfa * PEA / (Kmfa + PEA) - Vmna * PEA / (Kmna + PEA) - kass * PEA + kdis * PPRA - kpea * PEA
        dPA <- kpa + Vmfa * PEA / (Kmfa + PEA) + Vmna * PEA / (Kmna + PEA) - kdeg * PA
        dPPRA <- kass * PEA - kdis * PPRA

        res <- c(dNAPE, dPEA, dPA, dPPRA)
        list(res)
    })
}

t <- seq(0, 240, 1)
xstart <- c(
    NAPE = 0.3371102,
    PEA = 14.0134102,
    PA = 48.9999720,
    PPRA = 70.0670294
)

out <- ode(y = xstart, times = t, func = PEA_model, parms = para)
plot(out)

RS <- runsteady(y = xstart, func = PEA_model, parms = para, times = c(0, 1e18))

unlist(RS)
